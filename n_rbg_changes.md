# How `n_rbg.py` generalizes `rbg.py`

This note compares the original 0-surgery RBG search in `code/rbg.py` with
the positive-`n` search in `code/n_rbg.py`. Line numbers refer to the current
versions of the files and may move as the code changes.

The overall search strategy is unchanged:

1. Start with two knot exteriors.
2. Drill a curve from the first exterior to obtain a two-cusped blue-green
   exterior.
3. Fill one cusp and test whether the result is the second knot exterior.
4. Drill another curve to obtain a three-cusped candidate RBG exterior.
5. Recover a link diagram and test its special structure and framings.

## 1. RBG link representation and homology

### 0-surgery version

[`RedBlueGreenLink.__init__`](code/rbg.py#L215) constructs the blue and green
exteriors by leaving the corresponding component unfilled
([`code/rbg.py:215-226`](code/rbg.py#L215)). Its `_verify` method requires the
fully filled link to have

$$
H_1 \cong \mathbb Z,
$$

which appears as the elementary divisors `[0]` in
[`code/rbg.py:230-238`](code/rbg.py#L230).

### Generalized version

[`NRedBlueGreenLink.__init__`](code/n_rbg.py#L37) follows the same filling
construction but also stores `n` ([`code/n_rbg.py:37-52`](code/n_rbg.py#L37)).
Its `_verify` method requires the order of the filled link's first homology to
equal `n` ([`code/n_rbg.py:54-58`](code/n_rbg.py#L54)). The helper
[`_cyclic_homology_order`](code/n_rbg.py#L22) treats the three cases uniformly:

- `n = 0`: infinite cyclic homology, represented by `[0]`;
- `n = 1`: trivial homology, represented by no elementary divisors;
- `n > 1`: finite cyclic homology of order `n`.

### Relation to Henri's branch

Henri introduced `NRedBlueGreenLink`, stored `n`, and identified
`H_1 \cong \mathbb Z/n` as the replacement for the original `H_1 \cong
\mathbb Z` condition. The current code makes that verification active,
handles `n = 0` and `n = 1` explicitly, and fixes constructor/import issues
present in the branch implementation.

## 2. The special-link condition

### 0-surgery version

[`RedBlueGreenLink.is_super_special`](code/rbg.py#L264) checks that:

- the blue and green framings are `(0, 1)`
  ([`code/rbg.py:281-283`](code/rbg.py#L281));
- all three individual components are unknots
  ([`code/rbg.py:284-287`](code/rbg.py#L284));
- the red-blue and red-green sublinks are Hopf links
  ([`code/rbg.py:289`](code/rbg.py#L289)).

### Generalized version

[`NRedBlueGreenLink.is_n_special`](code/n_rbg.py#L65) retains those
diagrammatic checks in [`code/n_rbg.py:67-75`](code/n_rbg.py#L67). It then
adds Qin's signed framed linking-matrix condition

$$
n=-\det M_L.
$$

Here the framed linking matrix is

$$
M_L=
\begin{pmatrix}
r & \operatorname{lk}(R,B) & \operatorname{lk}(R,G) \\
\operatorname{lk}(R,B) & 0 & \operatorname{lk}(B,G) \\
\operatorname{lk}(R,G) & \operatorname{lk}(B,G) & 0
\end{pmatrix},
$$

where the diagonal entries are the red, blue, and green framings and the
off-diagonal entries are the pairwise linking numbers.

The code inserts the red framing and the two zero framings on the diagonal of
the linking matrix before taking its determinant
([`code/n_rbg.py:77-80`](code/n_rbg.py#L77)).

### Relation to Henri's branch

Henri added an `is_n_super_special` method and recognized the determinant
condition. The current implementation additionally puts the red surgery
coefficient into the matrix; `link.linking_matrix()` alone does not encode
that framing. Absolute determinant is sufficient for the homology condition
of an `|n|`-RBG link, but Qin's `n`-special condition retains the sign. The
name `is_n_special` follows Qin's terminology, while
[`code/n_rbg.py:82-83`](code/n_rbg.py#L82) preserves Henri's method name as an
alias.

## 3. Blue-green exteriors and their filling slopes

### 0-surgery version

[`BlueGreenExterior`](code/rbg.py#L292) stores a meridian and longitude on
each cusp ([`code/rbg.py:332-338`](code/rbg.py#L332)). Since the common
surgery coefficient is zero, `_verify` fills both cusps along their
longitudes and checks for `H_1 \cong \mathbb Z`
([`code/rbg.py:340-347`](code/rbg.py#L340)). It then checks the two ways of
recovering `S^3`:

- blue meridian plus green longitude
  ([`code/rbg.py:349-352`](code/rbg.py#L349));
- blue longitude plus green meridian
  ([`code/rbg.py:354-357`](code/rbg.py#L354)).

### Generalized version

[`NBlueGreenExterior`](code/n_rbg.py#L86) also stores `n`. Its constructor
computes the two integer-surgery slopes

$$
n\mu_B+\lambda_B, \qquad n\mu_G+\lambda_G
$$

in [`code/n_rbg.py:89-104`](code/n_rbg.py#L89). Its `_verify` method replaces
each longitude filling from the 0-surgery code with the corresponding
`n`-slope:

- both `n`-slopes must give homology of order `n`
  ([`code/n_rbg.py:106-110`](code/n_rbg.py#L106));
- blue meridian plus green `n`-slope must give `S^3`
  ([`code/n_rbg.py:112-115`](code/n_rbg.py#L112));
- blue `n`-slope plus green meridian must give `S^3`
  ([`code/n_rbg.py:116-119`](code/n_rbg.py#L116)).

### Relation to Henri's branch

Henri introduced `NBlueGreenExterior` and computed positive-`n` filling
coefficients. The current implementation isolates that calculation in
`_integer_slope`, enables the knot-exterior checks, and uses the recorded
meridian and longitude coordinates rather than assuming every peripheral
basis is literally `((1, 0), (0, 1))`.

## 4. The outer search for the second knot exterior

### 0-surgery version

The word-based search is
[`blue_green_exteriors_alt`](code/rbg.py#L606). After drilling a word from the
blue exterior, it fills the original cusp along `blue_long`
([`code/rbg.py:621-630`](code/rbg.py#L621)). It then applies a homology/volume
filter, looks for an isometry to the green exterior, and transports the green
meridian and longitude through the first isometry
([`code/rbg.py:631-642`](code/rbg.py#L631)).

The older dual-curve variant makes the same longitude filling in
[`blue_green_exteriors`, `code/rbg.py:579-603`](code/rbg.py#L579).

### Generalized version

[`n_blue_green_exteriors`](code/n_rbg.py#L176) performs the same search, but
the key filling is now

```python
candidate.dehn_fill(n_surgery_slope(blue, n), 0)
```

at [`code/n_rbg.py:188-193`](code/n_rbg.py#L188). The imported
`n_surgery_slope` uses the actual homological longitude, so this is
`n*meridian + longitude` even when the peripheral coordinates are not in the
standard basis.

After the homology and volume filters
([`code/n_rbg.py:194-197`](code/n_rbg.py#L194)), the generalized code tries
every returned isometry and transports both green peripheral curves
([`code/n_rbg.py:198-207`](code/n_rbg.py#L198)).

### Relation to Henri's branch

Henri made the central conceptual change from longitude filling to an
`n`-surgery filling and used increased precision when drilling. The current
code uses the shared `n_surgery_slope` helper instead of assuming the slope is
always `(n, 1)`, and it tries all peripheral maps rather than only the first
isometry.

## 5. Searching for and reconstructing the red component

### 0-surgery version

[`BlueGreenExterior.search_for_nice_dual_curves`](code/rbg.py#L369) drills
short dual curves, moves the drilled cusp into the red position, and calls
`_drilled_manifold_is_super_special`
([`code/rbg.py:369-382`](code/rbg.py#L369)).

That method:

- checks the red-blue and red-green Hopf-link fillings
  ([`code/rbg.py:517-528`](code/rbg.py#L517));
- enumerates possible red meridians and looks for an `S^3` filling
  ([`code/rbg.py:534-540`](code/rbg.py#L534));
- checks the blue and green longitudes
  ([`code/rbg.py:541-556`](code/rbg.py#L541));
- recovers, mirrors, and reorders the link diagram, then transports its
  framings ([`code/rbg.py:558-576`](code/rbg.py#L558)).

### Generalized version

[`NBlueGreenExterior.search_for_nice_dual_curves`](code/n_rbg.py#L121) keeps
the same drill-and-reindex structure
([`code/n_rbg.py:121-128`](code/n_rbg.py#L121)). The positive-`n` candidate
test performs the Hopf-link checks in
[`code/n_rbg.py:130-141`](code/n_rbg.py#L130), then enumerates red meridians
whose simultaneous filling gives `S^3` in
[`code/n_rbg.py:143-151`](code/n_rbg.py#L143).

Link recovery is separated into
[`NBlueGreenExterior._link_from_exterior`](code/n_rbg.py#L153). It performs
the diagram conversion and tries both the recovered diagram and its mirror.
For each orientation it performs component reordering in
[`code/n_rbg.py:153-165`](code/n_rbg.py#L153), then transports the red filling
slope and installs the diagrammatic zero framings in
[`code/n_rbg.py:166-173`](code/n_rbg.py#L166).

Trying both orientations is necessary only when the surgery coefficient is
signed. Mirroring negates the three-component linking matrix, so it changes
an `n`-special link into a `(-n)`-special link. The original zero-surgery code
could mirror unconditionally because this distinction disappears for `n=0`.

### Relation to Henri's branch

Henri adapted the same drilling pipeline and passed `n` into the resulting
link object. The current version factors reconstruction into a smaller
method, reuses the component-ordering and peripheral-isometry helpers from
`rbg.py`, enumerates every successful peripheral marking, and treats
unsuccessful reconstructions as rejected candidates rather than aborting the
whole search.

## 6. Determining the red framing

### 0-surgery version

In the original code, the red framing is the transported red longitude:

```python
new_red_long = maps[0] * vector((1, 0))
```

See [`code/rbg.py:568-574`](code/rbg.py#L568). This is appropriate for the
0-framed red component in the super-special 0-surgery setting.

### Generalized version

For positive `n`, the slope `(1, 0)` on the newly drilled red cusp is the
slope that refills the drilled surgery core and recovers the blue-green
exterior. It is therefore the red framing, not the red meridian. After the
link components have been ordered, the code transports this slope through
the peripheral isometry:

```python
red_frame = maps[0] * vector((1, 0))
```

This is implemented in [`code/n_rbg.py:166-170`](code/n_rbg.py#L166). The
blue and green framings are set to the diagrammatic zero slope `(0, 1)`, as
required for an `n`-special link. The resulting link is then checked by
`is_n_special`, including the condition `det(M_L) == -n`, in
[`code/n_rbg.py:171-173`](code/n_rbg.py#L171).

### Relation to Henri's branch

This interpretation corrects the earlier positive-`n` implementation, which
treated the drilled cusp's filling slope as a meridian and then searched for
an integral red framing using only the determinant. Transporting the actual
filling slope retains the peripheral information needed to distinguish the
correct red framing when more than one integral value has determinant of
absolute value `n`.

## 7. Public search API

### 0-surgery version

The original public generators are
[`blue_green_exteriors`](code/rbg.py#L579) and
[`blue_green_exteriors_alt`](code/rbg.py#L606). Callers iterate these
generators or use `next(...)` to obtain a link.

### Generalized version

[`n_blue_green_exteriors`](code/n_rbg.py#L176) is the positive-`n` generator.
[`find_n_special_rbg_link`](code/n_rbg.py#L214) is a convenience wrapper that
returns the first result or `None` rather than allowing `StopIteration` to
escape ([`code/n_rbg.py:214-217`](code/n_rbg.py#L214)).

### Relation to Henri's branch

Henri exposed a similar operation as `forms_super_special_NRBG_link`, but it
called `next(generator)` directly. Returning `None` makes the no-result case
explicit and is easier for tests and higher-level pipelines to handle.

## 8. Enumerating peripheral markings and recovering the requested pair

A three-cusped exterior can admit several `n`-special peripheral markings.
For the Example 3.2 exterior, the first valid marking found by the search
produces `6_1` and `K13n3663`, while a later marking produces the requested
pair `6_2` and `K13n3596`. Thus finding the correct unmarked exterior is not
enough: the search must retain the peripheral data and check the associated
blue and green knots.

The generalized search now implements this in two places:

1. [`code/n_rbg.py:128`](code/n_rbg.py#L128) uses `yield from` so that each
   drilled red curve can return every valid marking. The red-meridian loop at
   [`code/n_rbg.py:143-151`](code/n_rbg.py#L143) yields each successfully
   reconstructed `n`-special link instead of stopping at the first one.
2. The outer search checks that the reconstructed blue and green exteriors
   are isometric to the requested input exteriors before yielding the link
   ([`code/n_rbg.py:208-211`](code/n_rbg.py#L208)). This prevents an
   `n`-special marking of the right unmarked exterior but the wrong knot pair
   from becoming a false positive.

The old longitude comparison was also removed from the red-meridian loop.
The blue and green longitude coordinates on the intermediate blue-green
exterior change under the nonzero red surgery, so comparing them directly to
the zero longitudes of the reconstructed link incorrectly rejects the known
`n=1` example. The reconstruction now tests the diagrammatic zero framings
directly in [`code/n_rbg.py:166-173`](code/n_rbg.py#L166).

Together these changes make the Example 3.2 search an end-to-end passing
test: it returns a 1-special link whose blue knot is `6_2`, whose green knot
is `K13n3596`, and whose exterior is isometric to the saved diagram exterior.
Because the same knot pair has more than one 1-special exterior within the
search bounds, the example test supplies the saved exterior as an optional
target instead of assuming the first valid result is Qin's diagram.
