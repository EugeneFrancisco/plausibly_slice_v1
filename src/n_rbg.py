"""Search for n-special RBG links associated to a pair of n-friends."""

from __future__ import annotations

import snappy
from sage.all import matrix, vector

from find_n_friends import n_surgery_slope
from rbg import (
    invert_perm,
    is_Z_homology_solid_torus,
    is_hopf_link,
    is_hopf_link_exterior,
    is_three_sphere,
    is_unlink,
    isometry_preserving_curves,
    normalize_slope,
    reorder_link_components,
)


def _cyclic_homology_order(manifold) -> int:
    """Return |H_1|, using zero for infinite cyclic homology."""
    divisors = manifold.homology().elementary_divisors()
    if divisors == [0]:
        return 0
    if not divisors:
        return 1
    if len(divisors) == 1:
        return int(divisors[0])
    return -1


class NRedBlueGreenLink:
    """An |n|-RBG link, with components ordered red, blue, green."""

    def __init__(self, link, n: int, framings):
        if n < 0:
            raise ValueError("n must be non-negative")
        self.n = n
        self.link = link.copy()
        self.framings = [normalize_slope(s) for s in framings]
        self.exterior = self.link.exterior()
        r, b, g = self.framings

        E = self.exterior.copy()
        E.dehn_fill([r, (0, 0), g])
        self.blue_exterior = E.filled_triangulation()
        E = self.exterior.copy()
        E.dehn_fill([r, b, (0, 0)])
        self.green_exterior = E.filled_triangulation()
        self._verify()

    def _verify(self):
        E = self.exterior.copy()
        E.dehn_fill(self.framings)
        if _cyclic_homology_order(E) != self.n:
            raise ValueError("filled link does not have H_1 = Z/n")

    def __repr__(self):
        r, b, g = self.framings
        return (f"<{self.n}-RBG link of {len(self.link.crossings)} crossings "
                f"with r={r}, b={b}, g={g}>")

    def is_n_special(self) -> bool:
        """Check the stronger, diagrammatically testable n-special condition."""
        r, b, g = self.framings
        if b != (0, 1) or g != (0, 1) or r[1] != 1:
            return False
        L = self.link
        if any(not is_unlink(L.sublink([i])) for i in range(3)):
            return False
        if not (is_hopf_link(L.sublink([0, 1]))
                and is_hopf_link(L.sublink([0, 2]))):
            return False

        linking = matrix(L.linking_matrix())
        linking[0, 0] = r[0]
        linking[1, 1] = linking[2, 2] = 0
        return abs(int(linking.det())) == self.n

    # Backwards-compatible name used on henris-branch.
    is_n_super_special = is_n_special


class NBlueGreenExterior:
    """A two-cusped exterior encoding two knots with a common n-surgery."""

    def __init__(self, manifold, n, blue_merid, blue_long,
                 green_merid, green_long):
        self.manifold = snappy.Manifold(manifold)
        self.n = n
        self.blue_merid = normalize_slope(blue_merid)
        self.blue_long = normalize_slope(blue_long)
        self.green_merid = normalize_slope(green_merid)
        self.green_long = normalize_slope(green_long)
        self.blue_n = self._integer_slope(self.blue_merid, self.blue_long, n)
        self.green_n = self._integer_slope(self.green_merid, self.green_long, n)
        self._verify()

    @staticmethod
    def _integer_slope(meridian, longitude, n):
        return normalize_slope((n * meridian[0] + longitude[0],
                                n * meridian[1] + longitude[1]))

    def _verify(self):
        M = self.manifold.copy()
        M.dehn_fill([self.blue_n, self.green_n])
        if _cyclic_homology_order(M) != self.n:
            raise ValueError("mutual n-surgery has wrong homology")

        M = self.manifold.copy()
        M.dehn_fill([self.blue_merid, self.green_n])
        if not is_three_sphere(M):
            raise ValueError("blue meridian filling is not S^3")
        M = self.manifold.copy()
        M.dehn_fill([self.blue_n, self.green_merid])
        if not is_three_sphere(M):
            raise ValueError("green meridian filling is not S^3")

    def search_for_nice_dual_curves(self, max_segments=12, radius=6.0):
        """Yield n-special links found by drilling the red component."""
        for curve in self.manifold.dual_curves(max_segments):
            if curve.filled_length.real() > radius:
                break
            E = self.manifold.drill(curve)
            E._reindex_cusps([1, 2, 0])
            yield from self._drilled_manifold_is_n_special(E)

    def _drilled_manifold_is_n_special(self, manifold):
        E = manifold.copy()
        E1 = E.copy()
        E1.dehn_fill(self.blue_merid, 1)
        if not is_hopf_link_exterior(E1.filled_triangulation()):
            return None
        E2 = E.copy()
        E2.dehn_fill(self.green_merid, 2)
        if not is_hopf_link_exterior(E2.filled_triangulation()):
            return None
        if E.solution_type().startswith("contains"):
            return None

        for red_merid in E.short_slopes(12, first_cusps=[0])[0]:
            meridians = [red_merid, self.blue_merid, self.green_merid]
            S3 = E.copy()
            S3.dehn_fill(meridians)
            if not is_three_sphere(S3):
                continue
            answer = self._link_from_exterior(E, S3, meridians)
            if answer is not None:
                yield answer

    def _link_from_exterior(self, exterior, filled, meridians):
        L = filled.exterior_to_link(check_answer=False)
        L.simplify("global")
        L = L.mirror()
        iso = isometry_preserving_curves(
            exterior, meridians, L.exterior(), 3 * [(1, 0)])
        if iso is None:
            return None
        L = reorder_link_components(L, invert_perm(iso.cusp_images()))
        iso = isometry_preserving_curves(
            exterior, meridians, L.exterior(), 3 * [(1, 0)])
        if iso is None or iso.cusp_images() != [0, 1, 2]:
            return None
        maps = iso.cusp_maps()
        red_frame = maps[0] * vector((1, 0))
        try:
            answer = NRedBlueGreenLink(
                L, self.n, [red_frame, (0, 1), (0, 1)])
        except ValueError:
            return None
        return answer if answer.is_n_special() else None


def n_blue_green_exteriors(blue_exterior, green_exterior, n,
                           blue_merid=(1, 0), green_merid=(1, 0),
                           radius=4.0, max_segments=12):
    """Yield n-special RBG links whose associated knots are the inputs."""
    if n <= 0:
        raise ValueError("use rbg.blue_green_exteriors_alt when n = 0")
    blue = snappy.Manifold(blue_exterior)
    green = snappy.Manifold(green_exterior)
    blue_long = blue.homological_longitude()
    green_long = green.homological_longitude()
    green_volume = green.volume()

    curves = blue.length_spectrum(radius, include_words=True, grouped=False)
    for curve in curves:
        E = blue.drill_word(curve.word, bits_prec=100)
        candidate = E.copy()
        candidate.dehn_fill(n_surgery_slope(blue, n), 0)
        candidate = candidate.filled_triangulation()
        if not is_Z_homology_solid_torus(candidate):
            continue
        if abs(candidate.volume() - green_volume) >= 1e-8:
            continue
        isometries = green.is_isometric_to(candidate, True)
        for iso in isometries:
            cusp_map = iso.cusp_maps()[0]
            mapped_merid = normalize_slope(cusp_map * vector(green_merid))
            mapped_long = normalize_slope(cusp_map * vector(green_long))
            try:
                bge = NBlueGreenExterior(
                    E, n, blue_merid, blue_long, mapped_merid, mapped_long)
            except ValueError:
                continue
            for answer in bge.search_for_nice_dual_curves(max_segments):
                if (answer.blue_exterior.is_isometric_to(blue)
                        and answer.green_exterior.is_isometric_to(green)):
                    yield answer


def find_n_special_rbg_link(blue_exterior, green_exterior, n, **kwargs):
    """Return the first n-special RBG link found, or None."""
    return next(n_blue_green_exteriors(
        blue_exterior, green_exterior, n, **kwargs), None)
