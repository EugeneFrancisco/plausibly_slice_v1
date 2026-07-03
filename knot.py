"""
A ``Knot`` class for the Theorem 5.10 pipeline of Dunfield-Gong.

A knot is naturally more than a name: it carries an exterior, sometimes a
diagram, and -- as we run the pipeline -- *discovered* properties like its
Rasmussen s-invariants, its smooth-slice status, and the 0-friends and RBG
links relating it to other knots.  This module bundles all of that onto a
single object so those discoveries stick to the knot they describe.

The three stages of the paper map onto the class as follows:

    1. 0-friends: ``Knot.zero_surgery``, ``Knot.is_zero_friend``,
       ``Knot.find_zero_friends``.
    2. RBG links: ``Knot.search_for_rbg_link`` (live Section 5.11 search)
       and ``Knot.recorded_rbg_link`` (fast, from the paper's CSV), plus
       ``Knot.blue_green_knots`` to recover K_B, K_G from a link.
    3. s-invariants: ``Knot.rasmussen_s`` (shells out to KnotJob) and the
       static ``Knot.conclude_via_theorem_5_9``.

Computing the Rasmussen s-invariant uses KnotJob and needs Java >= 23; the
path to it is auto-detected (see ``_find_java``) and can be overridden with
the ``KNOTJOB_JAVA`` environment variable.

Everything here runs in the ``sage`` conda env (SnapPy + Sage).
"""

from __future__ import annotations

import os
import re
import csv
import sys
import ast
import shutil
import tempfile
import subprocess

import snappy

# The paper's code lives in ``code/``; make it importable before importing.
_HERE = os.path.dirname(os.path.abspath(__file__))
if os.path.join(_HERE, "code") not in sys.path:
    sys.path.insert(0, os.path.join(_HERE, "code"))

from rbg import (                                      # noqa: E402
    RedBlueGreenLink,
    BlueGreenExterior,
    blue_green_exteriors,
    blue_green_exteriors_alt,
    is_Z_homology_solid_torus,
    normalize_slope,
)
from sage.all import vector                             # noqa: E402
from snappy.exceptions import (                         # noqa: E402
    InsufficientPrecisionError,
    SnapPeaFatalError,
)
from find_0_friends import find_common_zero_surgery_via_words  # noqa: E402

DATA_CSV = os.path.join(_HERE, "data", "unknown_with_0-friend_final.csv")
KNOTJOB_JAR = os.path.join(_HERE, "tools", "knotjob", "KnotJob.jar")

# Numerical/geometric failures that SnapPy can raise on a single "bad"
# drilled curve during the Section 5.11 search.  None of these mean the
# pair has no RBG link -- they just mean *this* curve could not be
# evaluated at the working precision -- so the robust search skips the
# curve and moves on rather than aborting the whole (hours-long) run.
_SEARCH_ERRORS = (
    RuntimeError,
    InsufficientPrecisionError,
    SnapPeaFatalError,
    AssertionError,
    AttributeError,
    ValueError,
)


# ============================ KnotJob plumbing ============================

def _find_java() -> str:
    """
    Locate a Java >= 23 capable of running KnotJob.jar.

    Returns:
        Path to a suitable java executable.
    """
    if os.environ.get("KNOTJOB_JAVA"):
        return os.environ["KNOTJOB_JAVA"]
    # Prefer a conda env named jdk23 next to the active env.
    for base in [os.path.expanduser("~/miniforge3/envs/jdk23/bin/java"),
                 os.path.expanduser("~/miniconda3/envs/jdk23/bin/java"),
                 os.path.expanduser("~/anaconda3/envs/jdk23/bin/java")]:
        if os.path.exists(base):
            return base
    found = shutil.which("java")
    if found:
        return found
    raise RuntimeError("No java found; set KNOTJOB_JAVA to a Java >= 23.")


def _pd_to_knotjob(link: snappy.Link, name: str) -> str:
    """
    Format a SnapPy Link as a KnotJob planar-diagram line.

    Args:
        link: The SnapPy Link to format.
        name: Name to label the diagram with.

    Returns:
        A KnotJob 'PD' input line.
    """
    crossings = ", ".join(
        "X[" + ",".join(str(x) for x in c) + "]" for c in link.PD_code()
    )
    return f"{name} = PD {crossings}"


def _rasmussen_s_via_knotjob(
    link: snappy.Link,
    primes: tuple[int, ...] = (2, 3),
    name: str = "knot",
    max_heap: str = "8g",
) -> dict[int, int]:
    """
    Compute Rasmussen s-invariants of a knot diagram via KnotJob.

    Args:
        link: The knot as a SnapPy Link.
        primes: Fields F_p to compute over (p = 0 means the rationals).
        name: Name to label the diagram with.
        max_heap: Max Java heap size, e.g. '8g'.

    Returns:
        A dict {p: s_value} of s-invariants keyed by prime.
    """
    java = _find_java()
    if not os.path.exists(KNOTJOB_JAR):
        raise RuntimeError(f"KnotJob.jar not found at {KNOTJOB_JAR}")

    with tempfile.TemporaryDirectory() as tmp:
        txt = os.path.join(tmp, "knot.txt")
        with open(txt, "w") as f:
            f.write(_pd_to_knotjob(link, name) + "\n")

        cmd = [java, f"-Xmx{max_heap}", "-jar", KNOTJOB_JAR, txt]
        cmd += [f"-s{p}" for p in primes]
        cmd += ["-nf"]  # print to terminal, don't write a results file
        out = subprocess.run(cmd, capture_output=True, text=True).stdout

    # Lines look like:  "S-Invariant mod 3 : 2"
    result = {}
    for p, val in re.findall(r"S-Invariant mod (\d+)\s*:\s*(-?\d+)", out):
        result[int(p)] = int(val)
    return result


# ==================== Robust Section 5.11 search drivers ==================
# These mirror ``rbg.blue_green_exteriors`` /
# ``BlueGreenExterior.search_for_nice_dual_curves`` exactly, except a
# numerical failure on one drilled curve is caught and skipped instead of
# aborting the whole search.  They back ``Knot.search_for_rbg_link``.

def _super_special_from_bge(
    bge: BlueGreenExterior,
    max_segments: int = 12,
    radius: float = 6.0,
    verbose: bool = True,
):
    """
    Robust version of ``BlueGreenExterior.search_for_nice_dual_curves``.

    Drills the 2-cusped blue-green exterior a second time to find a red
    cusp making a super-special RBG link, yielding each hit.

    Args:
        bge: The blue-green exterior (paper's X = Y(eps_R, ., .)).
        max_segments: Segment cutoff for ``dual_curves``.
        radius: Filled-length cutoff on the drilled curves.
        verbose: If True, print each skipped curve.

    Yields:
        Super-special RedBlueGreenLink objects.
    """
    M = bge.manifold
    for curve in M.dual_curves(max_segments):
        if curve.filled_length.real() > radius:
            break
        try:
            E = M.drill(curve)
            # The drilled cusp becomes "red", so move it to the front.
            E._reindex_cusps([1, 2, 0])
            ans = bge._drilled_manifold_is_super_special(E)
        except _SEARCH_ERRORS as e:
            if verbose:
                print(f"    [skip curve {curve}: {type(e).__name__}]", flush=True)
            continue
        if ans is not None:
            yield ans


def _blue_green_exteriors_robust(
    blue_exterior: snappy.Manifold,
    blue_merid: tuple[int, int],
    green_exterior: snappy.Manifold,
    green_merid: tuple[int, int],
    max_segments: int = 12,
    verbose: bool = True,
):
    """
    Robust version of ``rbg.blue_green_exteriors``.

    Drills the blue exterior along curves in its dual 1-skeleton, keeps a
    drilling whose blue-longitude (0-)filling is a Z-homology solid torus
    isometric to the green exterior, packages the match as a
    ``BlueGreenExterior``, and hands it to ``_super_special_from_bge``.

    Args:
        blue_exterior: The exterior E_K we drill (the "blue" knot).
        blue_merid: Meridian slope of the blue exterior.
        green_exterior: The isometry target E_{K'} (the "green" knot).
        green_merid: Meridian slope of the green exterior.
        max_segments: Segment cutoff for ``dual_curves``.
        verbose: If True, print progress.

    Yields:
        Super-special RedBlueGreenLink objects.
    """
    blue_long = blue_exterior.homological_longitude()
    green_long = green_exterior.homological_longitude()
    green_vol = green_exterior.volume()
    for d in blue_exterior.dual_curves(max_segments):
        if verbose:
            print(f"  {d}", flush=True)
        try:
            E = blue_exterior.drill(d)
            E1 = E.copy()
            E1.dehn_fill(blue_long, 0)
            E1 = E1.filled_triangulation()
            if not (is_Z_homology_solid_torus(E1)
                    and abs(E1.volume() - green_vol) < 1e-8):
                continue
            if E1.solution_type().startswith("contains"):
                continue
            isos = green_exterior.is_isometric_to(E1, True)
            if not isos:
                continue
            A = isos[0].cusp_maps()[0]
            new_green_merid = normalize_slope(A * vector(green_merid))
            new_green_long = normalize_slope(A * vector(green_long))
            bge = BlueGreenExterior(
                snappy.Manifold(E), blue_merid, blue_long,
                new_green_merid, new_green_long,
            )
        except _SEARCH_ERRORS as e:
            if verbose:
                print(f"  [skip outer curve {d}: {type(e).__name__}]", flush=True)
            continue
        if verbose:
            print(f"  FOUND matching green exterior: volume "
                  f"{bge.manifold.volume()}", flush=True)
        yield from _super_special_from_bge(bge, verbose=verbose)


# ================================ The Knot ================================

class Knot:
    """
    A knot together with the properties we discover about it.

    A ``Knot`` can be built from a SnapPy name/isosig, a knot exterior
    (a one-cusped ``snappy.Manifold``), or a diagram (a ``snappy.Link``).
    Whatever you give, the others are computed lazily on demand.

    Discovered properties, filled in as methods run:

        s_invariants:    dict {p: s} of Rasmussen s-invariants.
        smoothly_slice:  None (unknown), or a bool once determined.
        zero_friends:    list of other Knots sharing this 0-surgery.
        rbg_link:        a RedBlueGreenLink relating this knot to a friend.

    Examples:

        >>> K = Knot('K11n34')
        >>> K.exterior().num_cusps()
        1
        >>> K.smoothly_slice is None
        True
    """

    def __init__(self, identifier=None, *, exterior=None, diagram=None,
                 name=None):
        """
        Args:
            identifier: A name/isosig (str), an exterior (snappy.Manifold),
                or a diagram (snappy.Link).  Optional if you pass the
                keyword forms instead.
            exterior: A one-cusped knot exterior as a Manifold.
            diagram: A knot diagram as a SnapPy Link.
            name: A human-readable name (defaults to the identifier if it
                was a string).
        """
        self.name = name
        self._exterior = None
        self._diagram = None

        if isinstance(identifier, str):
            self.name = self.name or identifier
        elif isinstance(identifier, snappy.Link):
            self._diagram = identifier
        elif isinstance(identifier, snappy.Manifold):
            self._exterior = identifier.copy()
        elif identifier is not None:
            raise TypeError(f"Cannot build a Knot from {identifier!r}")

        if exterior is not None:
            self._exterior = exterior.copy()
        if diagram is not None:
            self._diagram = diagram

        # Discovered properties.
        self.s_invariants: dict[int, int] = {}
        self.smoothly_slice: bool | None = None
        self.zero_friends: list[Knot] = []
        self.rbg_link: RedBlueGreenLink | None = None

    # ---- construction helpers ----------------------------------------

    @classmethod
    def coerce(cls, knot) -> "Knot":
        """
        Return ``knot`` if it is already a Knot, else wrap it in one.

        Args:
            knot: A Knot, or anything a Knot can be built from.

        Returns:
            A Knot instance.
        """
        return knot if isinstance(knot, cls) else cls(knot)

    @property
    def label(self) -> str:
        """A short human-readable label for printing."""
        return self.name or repr(self._exterior) or "knot"

    def __repr__(self) -> str:
        bits = [self.label]
        if self.s_invariants:
            bits.append(f"s={self.s_invariants}")
        if self.smoothly_slice is not None:
            bits.append("slice" if self.smoothly_slice else "not-slice")
        return f"<Knot {' '.join(bits)}>"

    # ---- geometry ----------------------------------------------------

    def exterior(self) -> snappy.Manifold:
        """
        The knot exterior as a one-cusped SnapPy Manifold (cached).

        Returns:
            The exterior Manifold.
        """
        if self._exterior is None:
            if self.name is not None:
                self._exterior = snappy.Manifold(self.name)
            elif self._diagram is not None:
                self._exterior = self._diagram.exterior()
            else:
                raise ValueError("Knot has no name, exterior, or diagram")
        return self._exterior

    def _search_exterior(self) -> snappy.Manifold:
        """A fresh one-cusped exterior copy for the RBG search."""
        M = self.exterior().copy()
        if M.num_cusps() != 1:
            raise ValueError(f"Expected a one-cusped knot exterior, got {M}")
        return M

    def diagram(self, simplify: bool = True) -> snappy.Link:
        """
        A diagram of the knot as a SnapPy Link (cached).

        Args:
            simplify: If True, return a simplified copy of the diagram.

        Returns:
            The knot as a SnapPy Link.
        """
        if self._diagram is None:
            self._diagram = self.exterior().exterior_to_link()
        link = self._diagram
        if simplify:
            link = link.copy()
            link.simplify("global")
        return link

    def zero_surgery(self) -> snappy.Manifold:
        """
        The 0-surgery of the knot as a closed Manifold.

        Returns:
            The 0-surgery as a closed SnapPy Manifold.
        """
        M = self.exterior().copy()
        M.dehn_fill((0, 1))
        return M

    # ---- Part 1: 0-friends -------------------------------------------

    def is_zero_friend(self, other, tries: int = 5) -> bool:
        """
        Certify whether ``self`` and ``other`` are 0-friends.

        Two knots are 0-friends when their 0-surgeries are the same
        closed manifold.  On success both knots record each other in
        their ``zero_friends`` list.

        Args:
            other: Another Knot (or anything Knot.coerce accepts).
            tries: Number of randomized retries for the isometry check.
                Each retry quadruples the triangulation size.

        Returns:
            True iff the two 0-surgeries are certified isometric.
        """
        other = Knot.coerce(other)
        A, B = self.zero_surgery(), other.zero_surgery()
        friends = False
        for _ in range(tries):
            try:
                friends = bool(A.is_isometric_to(B))
                break
            except RuntimeError:
                A.randomize(); B.randomize()
        if friends:
            self._add_zero_friend(other)
            other._add_zero_friend(self)
        return friends

    def find_zero_friends(
        self, max_len: float = 3.0
    ) -> list[tuple[str, complex, float, str]] | None:
        """
        Search for 0-friends via short closed geodesics.

        Thin wrapper around the paper's main search routine: it drills the
        short geodesics of the hyperbolic 0-surgery and keeps those whose
        result is again a knot exterior.  Each hit is recorded as a Knot
        in ``self.zero_friends``.

        Args:
            max_len: Maximum geodesic length to search; larger is slower.

        Returns:
            A list of 4-tuples (word, complex_length, volume, isosig), or
            None.
        """
        hits = find_common_zero_surgery_via_words(self.exterior(), max_len)
        for hit in hits or []:
            isosig = hit[-1]
            self._add_zero_friend(Knot(isosig))
        return hits

    def _add_zero_friend(self, other: "Knot") -> None:
        """Record ``other`` as a 0-friend of this knot (dedup by label)."""
        if other is self:
            return
        if all(f.label != other.label for f in self.zero_friends):
            self.zero_friends.append(other)

    # ---- Part 2: RBG links -------------------------------------------

    @staticmethod
    def search_for_rbg_link(
        knot_K,
        knot_Kprime,
        meridian: tuple[int, int] = (1, 0),
        use_words: bool = False,
        radius: float = 4.0,
        try_both_orders: bool = True,
        robust: bool = True,
        verbose: bool = True,
    ) -> RedBlueGreenLink | None:
        """
        Search for a super-special RBG link realizing a 0-friend pair.

        This is the live version of the paper's Stage D (Dunfield-Gong
        Section 5.11).  It drills the blue exterior E_K along curves in
        its dual 1-skeleton, looks for a drilling whose blue-longitude
        filling is isometric to the green exterior E_{K'}, and from each
        match drills once more searching for a red meridian that makes the
        whole thing a super-special RBG link.

        The section 5.11 search is NOT symmetric in K and K': the roles of
        "blue" (the knot we drill) and "green" (the isometry target) are
        not interchangeable, so by default we retry with the roles swapped
        if the first order finds nothing.

        On success the RBG link is recorded on *both* knots (``rbg_link``)
        and they are marked as 0-friends of each other.

        Args:
            knot_K: One knot of the pair (a Knot, name/isosig, or Manifold).
            knot_Kprime: The other knot of the pair.
            meridian: Meridian slope of each exterior; (1, 0) is standard.
            use_words: If True, drill closed geodesics by word (the weaker
                Section 5.12 variant); if False, drill dual curves.
            radius: Geodesic-length cutoff used only when ``use_words``.
            try_both_orders: If True, also try the (K', K) role-swap.
            robust: If True, use the exception-tolerant drivers that skip
                numerically-bad curves; if False, use the exact
                ``rbg.blue_green_exteriors`` generator.  Ignored when
                ``use_words`` is True.
            verbose: If True, print progress.

        Returns:
            The first super-special RedBlueGreenLink found, or None.
        """
        knot_K = Knot.coerce(knot_K)
        knot_Kprime = Knot.coerce(knot_Kprime)
        E_K = knot_K._search_exterior()
        E_Kprime = knot_Kprime._search_exterior()

        orders = [(E_K, E_Kprime, "K -> K'")]
        if try_both_orders:
            orders.append((E_Kprime, E_K, "K' -> K"))

        for blue, green, lbl in orders:
            if verbose:
                print(f"\n[search] blue={blue}  green={green}   ({lbl})",
                      flush=True)
            if use_words:
                gen = blue_green_exteriors_alt(
                    blue, meridian, green, meridian, radius=radius
                )
            elif robust:
                gen = _blue_green_exteriors_robust(
                    blue, meridian, green, meridian, verbose=verbose
                )
            else:
                gen = blue_green_exteriors(blue, meridian, green, meridian)

            for rbg in gen:
                # The generators only yield links that already passed the
                # super-special check; re-assert it as a cheap guard.
                if rbg.is_super_special():
                    if verbose:
                        print(f"[search] SUPER-SPECIAL: {rbg}", flush=True)
                    knot_K.rbg_link = rbg
                    knot_Kprime.rbg_link = rbg
                    knot_K._add_zero_friend(knot_Kprime)
                    knot_Kprime._add_zero_friend(knot_K)
                    return rbg

        if verbose:
            print("[search] no super-special RBG link found", flush=True)
        return None

    @staticmethod
    def recorded_rbg_link(base_knot: str) -> RedBlueGreenLink:
        """
        Load the precomputed super-special RBG link for a base knot.

        Reads the paper's CSV (``data/unknown_with_0-friend_final.csv``)
        for the DT code and framing recorded against ``base_knot`` and
        rebuilds the RedBlueGreenLink.  This is the fast alternative to
        ``search_for_rbg_link``.

        Args:
            base_knot: The base knot name to look up.

        Returns:
            The recorded RedBlueGreenLink.
        """
        row = _load_csv_row(base_knot)
        framing = ast.literal_eval(row["framing"])
        link = snappy.Link(row["RBG_DT"])
        return RedBlueGreenLink(link, framing)

    @staticmethod
    def blue_green_knots(
        rbg: RedBlueGreenLink, simplify: bool = True
    ) -> tuple["Knot", "Knot"]:
        """
        Recover the two knots K_B, K_G encoded by an RBG link.

        Each returned Knot carries both its diagram (from the link) and
        its exterior (from the link's blue/green fillings).

        Args:
            rbg: The RedBlueGreenLink to extract from.
            simplify: If True, simplify each diagram first.

        Returns:
            A tuple (K_B, K_G) of Knots.
        """
        blue, green = rbg.blue_knot(), rbg.green_knot()
        if simplify:
            blue, green = blue.copy(), green.copy()
            blue.simplify("global")
            green.simplify("global")
        k_b = Knot(diagram=blue, exterior=rbg.blue_exterior, name="K_B")
        k_g = Knot(diagram=green, exterior=rbg.green_exterior, name="K_G")
        k_b.rbg_link = rbg
        k_g.rbg_link = rbg
        return k_b, k_g

    # ---- Part 3: s-invariants ----------------------------------------

    def rasmussen_s(
        self,
        primes: tuple[int, ...] = (2, 3),
        max_heap: str = "8g",
    ) -> dict[int, int]:
        """
        Compute this knot's Rasmussen s-invariants via KnotJob.

        The result is stored in ``self.s_invariants`` (merged with any
        previously computed values) and returned.

        Args:
            primes: Fields F_p to compute over (p = 0 means the rationals).
            max_heap: Max Java heap size, e.g. '8g'.

        Returns:
            A dict {p: s_value} keyed by prime.
        """
        result = _rasmussen_s_via_knotjob(
            self.diagram(), primes=primes, name=self.label, max_heap=max_heap
        )
        self.s_invariants.update(result)
        return result

    @staticmethod
    def conclude_via_theorem_5_9(
        knot_b: "Knot", knot_g: "Knot"
    ) -> tuple[bool, str | None]:
        """
        Apply Theorem 5.9 to a pair of 0-friends via their s-invariants.

        If either knot has a nonzero s-invariant (over some field), then
        by Theorem 5.9 *both* knots fail to be smoothly slice; this sets
        ``smoothly_slice = False`` on both.

        Args:
            knot_b: One knot of the pair (K_B), with ``s_invariants`` set.
            knot_g: The other knot (K_G), with ``s_invariants`` set.

        Returns:
            A tuple (is_not_slice, witness_string_or_None).
        """
        for knot in (knot_b, knot_g):
            for p, val in knot.s_invariants.items():
                if val != 0:
                    field = "Q" if p == 0 else f"F_{p}"
                    knot_b.smoothly_slice = False
                    knot_g.smoothly_slice = False
                    return True, f"s_{{{field}}}({knot.label}) = {val} != 0"
        return False, None

    # ---- End-to-end driver (Parts 1-3) -------------------------------

    def conclude_slice_status_via_rbg(
        self,
        partner,
        primes: tuple[int, ...] = (2, 3),
        verbose: bool = True,
    ) -> bool:
        """
        Settle the smooth-slice status of this knot and a 0-friend via RBG.

        Runs the full Theorem 5.10 pipeline for the pair (``self``,
        ``partner``) end to end, combining all three stages of the paper:

            1. Search live (Section 5.11) for a super-special RBG link
               realizing the pair as 0-friends (``search_for_rbg_link``).
            2. Recover the two simpler knots K_B, K_G the link encodes
               (``blue_green_knots``).
            3. Compute their Rasmussen s-invariants (``rasmussen_s``) and
               apply Theorem 5.9 (``conclude_via_theorem_5_9``).

        By Theorem 5.9, a non-zero s-invariant on *either* K_B or K_G shows
        that *both* fail to be smoothly slice.  Because K_B and K_G are just
        simpler representatives of ``self`` and ``partner``, that verdict
        transfers to the pair: on success we set ``smoothly_slice = False``
        on both ``self`` and ``partner`` (and on K_B, K_G).

        Args:
            partner: The 0-friend to pair with (a Knot, name/isosig, or
                Manifold).
            primes: Fields F_p to compute s over (p = 0 means the rationals).
            verbose: If True, print progress mirroring the demonstration.

        Returns:
            True if the slice status was concluded (both knots shown not
            smoothly slice); False if no super-special RBG link was found or
            every computed s-invariant vanished.
        """
        partner = Knot.coerce(partner)

        # ---- Parts 1 & 2: find the super-special RBG link --------------
        if verbose:
            print(f"[1/3] Searching for a super-special RBG link for "
                  f"({self.label}, {partner.label})...", flush=True)
        rbg = Knot.search_for_rbg_link(self, partner, verbose=verbose)
        if rbg is None:
            if verbose:
                print("      No super-special RBG link found; cannot "
                      "conclude slice status.", flush=True)
            return False
        if verbose:
            print(f"      Found: {rbg}", flush=True)

        # ---- Recover the two knots the link encodes -------------------
        k_b, k_g = Knot.blue_green_knots(rbg)
        if verbose:
            print(f"[2/3] Encoded knots  K_B ({len(k_b.diagram().PD_code())} "
                  f"crossings),  K_G ({len(k_g.diagram().PD_code())} "
                  f"crossings)", flush=True)

        # ---- Part 3: s-invariants and Theorem 5.9 ---------------------
        if verbose:
            print("[3/3] Computing Rasmussen s-invariants...", flush=True)
        k_b.rasmussen_s(primes=primes)
        k_g.rasmussen_s(primes=primes)
        if verbose:
            print(f"      s(K_B) = {k_b.s_invariants}   "
                  f"s(K_G) = {k_g.s_invariants}", flush=True)

        not_slice, witness = Knot.conclude_via_theorem_5_9(k_b, k_g)
        if not not_slice:
            if verbose:
                print("      All computed s-invariants vanish; slice status "
                      "not concluded.", flush=True)
            return False

        # Transfer the verdict to the pair (same knots, larger diagrams).
        self.smoothly_slice = False
        partner.smoothly_slice = False
        self._add_zero_friend(partner)
        partner._add_zero_friend(self)
        if verbose:
            print(f"      {witness}", flush=True)
            print(f"      => by Theorem 5.9, {self.label} and "
                  f"{partner.label} are NOT smoothly slice.", flush=True)
        return True


def _load_csv_row(base_knot: str) -> dict[str, str]:
    """
    Fetch the RBG record for a base knot from the paper's CSV.

    Args:
        base_knot: The base knot name to look up.

    Returns:
        The matching CSV row as a dict.
    """
    with open(DATA_CSV) as f:
        for row in csv.DictReader(f):
            if row["base_knot"] == base_knot:
                return row
    raise KeyError(f"{base_knot} not in {DATA_CSV}")


if __name__ == "__main__":
    import doctest
    print(doctest.testmod())
