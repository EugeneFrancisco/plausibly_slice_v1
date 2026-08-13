# Obstructed n-friends in `data/n_friends(master version).csv`

**Prepared:** 2026-08-09 · **Source data:** `data/n_friends(master version).csv` (3697 rows)
**Companion machine-readable file:** [Data/results/obstructed_pairs.csv](obstructed_pairs.csv) (42 rows, includes both PD codes)

---

## 1. Headline

Scanning all 3697 rows against the obstruction criterion below:

| | count |
|---|---|
| Rows where a recorded invariant obstructs | **42** |
| Distinct census knots those 42 rows sit over | **12** |
| …of which the census knot's sliceness was **already known** (`slice = -1`) | **12 (all of them)** |
| …of which the census knot's sliceness is **open** (`slice = 0`) | **0** |

**So: the current data contains zero new sliceness results, but 42 fully-validated
test cases for the n-special RBG construction.** Every obstructing friend sits over a
census knot that is already known to be non-slice — and, worse for novelty, over a knot
whose *own* τ and s already obstruct it directly (see §4).

That is not a dead end. It is exactly the right input for the next step: these 12 knots
are the cases where you can **build the n-special RBG link and check that the machinery
reproduces a known answer** before trusting it on an open knot. §6 then lists where new
results are most likely to come from.

I recomputed the criterion from scratch rather than trusting the `obstructs` column; the
two agree exactly on all 42 rows. Two rows in the file have `verification = False`
(`id_num` 1442, 3455); neither is among the 42.

I also checked `Data/n_friends(Henri's version).csv` at `origin/main`, which carries **172
rows not yet merged into the master version** (only 17 of them have any invariant
computed). **None of the 172 obstructs.** So merging Henri's file will not change the
count in this report — but it should still be merged before the next invariant run.

---

## 2. Why an invariant on the friend obstructs the original knot

Reference: **Qianhe Qin, "An RBG construction of integral surgery homeomorphisms",
[arXiv:2308.04681](https://arxiv.org/abs/2308.04681)**, generalizing Manolescu–Piccirillo's
0-surgery RBG construction.

The chain of reasoning, as it applies to a row of this CSV:

1. **The row records a surgery homeomorphism.** The row's knot `K'` (stored in
   `knot_PD_code`) and the census knot `K` (`n_friend_PD_code`) satisfy
   `S³_n(K) ≅ S³_n(K')` for the row's `n`. This is what the search verified
   (`verification = True`).

2. **Theorem 1.2** says any such ⟨n⟩-surgery homeomorphism is realized by an `|n|`-RBG
   link `L_φ` producing the framed knots `{(K, n), (K', n)}`. If that link can be put in
   **n-special** form (Def. 1.3: `b = g = 0`, `n = −det(M_L)`, and the isotopies
   `R ∪ B ≅ R ∪ μ_R ≅ R ∪ G`), then `K = K_B` and `K' = K_G` are the two knots of an
   n-special RBG link. **This step is the work being handed off — it is not in the CSV.**

3. **Theorem 1.4(a)** then gives, for a `k`-special RBG link `L = {(R,r),(B,0),(G,0)}`:
   > if `R` is `r`-slice in some `#^m ℂℙ²` and `K_B` is `k`-slice in some `#^l ℂℙ̄²`,
   > then `s(K_G) ≤ k − √k`.

   (Def. 5.1: `K` is *n-slice* in `X°` if it bounds a properly embedded disk of
   self-intersection `−n`.)

4. **Sliceness feeds the hypothesis.** If `K_B` is smoothly slice in `B⁴`, tubing its
   slice disk with `k` exceptional spheres makes it `k`-slice in `#^k ℂℙ̄²`. So
   *`K_B` slice ⟹ `s(K_G) ≤ k − √k`*, and contrapositively

   > **`s(K_G) > k − √k` ⟹ `K_B` is not smoothly slice.**

   The pipeline applies the same shape of bound with `2τ` and `2ν` in place of `s`
   (see the caveat on ν in §5).

5. **Negative n is handled by mirroring.** `S³_{−m}(K) ≅ S³_{−m}(K')` is the same as
   `S³_m(mK) ≅ S³_m(mK')`, and `mK` is slice iff `K` is. So a row with `n < 0` is the
   `|n|` statement applied to the mirrors.

### The operative criterion

Let `m = |n|` and `thr = m − √m`. A row obstructs iff any of:

```
 s(K'_mirrored)  > thr      s_2(K'_mirrored) > thr      s_3(K'_mirrored) > thr
2·τ(K'_mirrored) > thr     2·ν(K'_mirrored) > thr
```

where "mirrored" means: take the mirror of the friend when `n < 0`, leave it alone when
`n > 0`. Thresholds:

| \|n\| | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| `thr` | 0 | 0.586 | 1.268 | 2 | 2.764 |

Since `s` and `2τ` are always even integers for knots, this collapses to a very simple rule:

* **|n| ≤ 3:** need `s ≥ 2` or `τ ≥ 1`.
* **|n| ≥ 4:** need `s ≥ 4` or `τ ≥ 2`.

This is why `|n| ∈ {4, 5}` has produced nothing: the bar jumps from `τ ≥ 1` to `τ ≥ 2`,
and `τ ≥ 2` simply does not occur anywhere in the current data.

---

## 3. ⚠️ Two data-layout traps — read before touching the CSV

### Trap 1: the column names are inverted relative to what you'd guess

`n_friend_name`, `n_friend_index`, `n_friend_PD_code`, and `slice` all describe the
**original census knot** from `data/plausibly_unknown.csv` (`n_friend_index` is its
1-based row number there). The row's *own* `knot_PD_code`, `num_crossings`, `volume`,
`s`, `s_2`, `s_3`, `nu`, `tau` describe the **found friend**.

Read `n_friend_name` as *"this row's knot is an n-friend **of** …"*.
Verified against `code/n_friends_search.py:216-224` and the PD codes read by
`code/compute_s_invariant.py:224` / `code/compute_floer_invariants.py:81`.

### Trap 2: `tau`/`nu` and `s`/`s_2`/`s_3` are stored in *different* sign conventions

This is the one that will silently produce fake results.

* `code/compute_floer_invariants.py:87-90` **mirrors the knot before computing HFK**
  when `n < 0`. So for `n < 0` rows, the `tau` and `nu` columns already hold the
  invariants **of the mirror**. Test them directly: `2·tau > thr`. **Do not flip again.**
* `code/compute_s_invariant.py` does **not** mirror. The `s`, `s_2`, `s_3` columns hold
  the un-mirrored friend's invariants, and the sign flip is applied at check time
  (`s → −s` when `n < 0`).

Confirmed independently from the data: among the 13 rows with `n < 0` where both `s` and
`τ` are present, **all 13** have `s = +2, τ = −1`. That is only consistent if `τ` is the
mirror's (`s(mK) = −2 = 2τ(mK)`); under a common convention it would read `s ≠ 2τ` in
13/13 rows, which is not credible.

**Consequence if you get this wrong:** flipping `τ` for `n < 0` yields 28 spurious
"new obstructions" at `slice = 0`, all of the form `n = −1 or −2, τ = −1`. They are
artifacts of double-mirroring. Under the correct convention `2τ = −2 < thr` and nothing
fires.

---

## 4. The 42 obstructing pairs

Grouped by the census knot to be obstructed. Every row's `n` and `id_num` is given so you
can pull the PD codes straight from
[obstructed_pairs.csv](obstructed_pairs.csv) or from the master CSV.

| Census knot `K` (to obstruct) | index in `plausibly_unknown.csv` | cr | own τ | own s₀ | # friends | rows (`id_num`) |
|---|---|---|---|---|---|---|
| `18nh_00182176` | 820 | 18 | 1 | 2 | 7 | 705(n=1), 706(n=1), 707(n=1), 708(n=1), 709(n=1), 710(n=2), 711(n=3) |
| `17nh_0024597` | 226 | 17 | 1 | 2 | 6 | 164(n=1), 165(n=1), 166(n=2), 167(n=2), 168(n=3), 169(n=3) |
| `18nh_00588827` | 933 | 18 | −1 | −2 | 5 | 3027(n=−2), 3028(n=−2), 3029(n=−1), 3030(n=−1), 3031(n=−1) |
| `18nh_00778169` | 963 | 18 | 1 | 2 | 4 | 911(n=1), 912(n=1), 913(n=1), 914(n=1) |
| `18nh_00244364` | 834 | 18 | 1 | 2 | 4 | 726(n=1), 727(n=2), 728(n=2), 729(n=3) |
| `18nh_00061316` | 757 | 18 | 1 | 2 | 4 | 558(n=1), 560(n=2), 561(n=2), 562(n=3) |
| `17nh_0466454` | 328 | 17 | 1 | 2 | 3 | 293(n=1), 294(n=2), 295(n=3) |
| `K15n77799` | 31 | 15 | 1 | 2 | 2 | 20(n=1), 21(n=3) |
| `18nh_00067017` | 763 | 18 | 1 | 2 | 2 | 573(n=1), 574(n=3) |
| `18nh_00004122` | 725 | 18 | −1 | −2 | 2 | 2649(n=−3), 2650(n=−1) |
| `17nh_0920835` | 390 | 17 | −1 | −2 | 2 | 2521(n=−3), 2522(n=−2) |
| `17nh_0013937` | 222 | 17 | 1 | 2 | 1 | 160(n=3) |

Notes on the table:

* **All 12 have `slice = -1`** in `plausibly_unknown.csv` (already known not smoothly slice).
* **All 12 are obstructed by their own invariants** (`τ = ±1`, `s₀ = ±2`), so even the
  *method* is not producing an independent proof here — the friend just inherits the same
  τ and s as the knot it is a friend of. In every one of the 42 rows,
  `τ(friend, mirrored) = |τ(K)| = 1`. Concordance-invariant agreement across an n-surgery
  homeomorphism is itself a sanity check that the surgery data is right.
* Only 4 of the 42 rows have `s` computed at all (`id_num` 20, 21, 166, 168); the other
  38 fire on `2τ` and `2ν` only.
* The sign pattern is consistent: `n > 0` rows have `τ(friend) = +1`; `n < 0` rows have
  `τ(friend, already mirrored) = +1` sitting over a census knot with `τ(K) = −1`.

### What to do with these

Use them as **the validation suite** for the n-special RBG step. For each pair,
attempt to realize the recorded surgery homeomorphism as an n-special RBG link
(`code/n_rbg.py` is the generalization of Dunfield–Gong's `rbg.py` for this). Success
means Theorem 1.4 reproduces a known `slice = -1`; failure isolates whether the blocker is
the n-special condition itself or the `R` is `r`-slice hypothesis. The best starting
points are the smallest ones:

* **`K15n77799` / `id_num` 20, n = 1** — census knot only 15 crossings, friend 39 crossings,
  both `s` and `τ` computed and obstructing. Smallest case in the set.
* **`id_num` 168, n = 3** (over `17nh_0024597`) — friend is only 33 crossings, full
  invariant set present.
* **`id_num` 166, n = 2** (same census knot) — 35-crossing friend, full invariant set.

---

## 5. Caveats on the criterion itself

* **The ν criterion is unverified.** `code/compute_floer_invariants.py:124` carries the
  comment *"double check that this is correct"*. Qin's Theorem 1.4 is stated for `s`;
  the `2τ` and `2ν` forms are the pipeline's extrapolation. Note that ν is **not**
  anti-symmetric under mirroring (`ν(mK) ≠ −ν(K)` in general) — which is precisely why
  the Floer code mirrors the knot rather than negating the output. In the current data ν
  never fires alone: every ν hit is accompanied by a τ hit (across the whole file,
  `(τ,ν)` only takes the values `(−1,0)`, `(0,0)`, `(1,1)`), so nothing in §4 depends on it.
  Still worth confirming with the advisor before it is used as a sole obstruction.
* **Theorem 1.4 has a hypothesis on `R` that the CSV cannot see** — that `R` is `r`-slice
  in some `#^m ℂℙ²`. It only becomes checkable once the RBG link exists. Treat every pair
  in §4 as a *candidate*, not a completed obstruction.
* **Theorem 1.4(b)** gives the alternative bound `s(K_G) ≤ k+1 − √(k+1)` under a different
  hypothesis on `R`. That is a strictly *weaker* obstruction (higher threshold) and is not
  used anywhere in the current pipeline. If the (a)-hypothesis fails for some link, check
  whether (b) applies with the shifted threshold: at `|n| = 1` it needs `s ≥ 2` still,
  but at `|n| = 3` it needs `s ≥ 4`.

---

## 6. Where new results would actually come from

The obstruction search is not limited by the criterion — it is limited by **missing
invariant computations**.

| | rows |
|---|---|
| `slice = 0` rows (census knot's sliceness open) | 3494 (over 937 distinct knots) |
| …with `τ`, `ν` computed | 916 |
| …with `s` computed | 66 |
| …with **no invariants at all** | **2578** |

Three concrete queues, in descending order of expected value:

### Queue A — the |n| = 5 near-miss cluster (highest value, cheapest)

**46 rows over 38 distinct `slice = 0` census knots have a friend with `τ = 1`, and every
single one of them is at `|n| = 5`, where `thr = 2.764 > 2 = 2τ`.** They miss by one step
of the threshold table.

These 38 knots are the most promising in the whole dataset: they demonstrably admit
friends with non-vanishing τ. The action is **not** to compute more invariants on the
5-friends — it is to **re-run the n-friend search on those 38 census knots at
|n| ≤ 3**, where `τ ≥ 1` suffices. Several already have `|n| ≤ 3` friends recorded with
`τ = 0`, but coverage is patchy (e.g. `17nh_0090600` and `18nh_01237232` have *no*
`|n| ≤ 3` friend at all; `18nh_00257425` and `18nh_00547552` have small-`n` friends with
τ never computed). Full list of the 38, with their per-`n` coverage, is reproducible from
the master CSV by filtering `slice == 0 and tau >= 1`.

Alternatively, compute `s` on the 46 five-friends: if any has `s = 4` while `τ = 1`, it
fires at `|n| = 5`. Lower probability, but the computation is already scripted.

### Queue B — `τ` known and 0, `s` unknown, `|n| ≤ 3`

693 rows. `s ≠ 2τ` does happen in this data (4 rows have `s = ∓2` with `τ = 0`), so an
`s` computation on a `τ = 0` friend can still fire. 43 of these are ≤ 45 crossings —
tractable for KnotJob. Start there.

### Queue C — 2578 rows with nothing computed

2024 of these are at `|n| ≤ 3` (where the bar is lowest), but they are overwhelmingly
large: only 10 are ≤ 60 crossings, and 742 exceed 200 crossings. The bottleneck is
diagram simplification, not the invariant computation — `code/simplify_master.py` and
`super_simplify_table` in `code/n_friends_search.py` are the relevant tools. Worth running
a bulk simplification pass before any more invariant computation is scheduled.

---

## 7. Bugs found while auditing (worth fixing, none affect §4)

1. **`check_obstruction()` in `code/compute_floer_invariants.py:134-177` is a no-op.**
   It does `row = data.iloc[i]` (a *copy*), mutates `row["obstructs"]`, then writes
   `data` back out — so the `obstructs` column is never actually updated by this
   function. Verified: `df.iloc[0]['a'] = 99` does not propagate to `df`. Fix with
   `data.at[i, "obstructs"] = obstructs`.

   The column happens to be correct anyway, because the per-row `compute_s_invariants`
   and `compute_floer_homology` set `obstructs` at compute time using the same logic.
   But anyone who re-runs `check_obstruction()` expecting a refresh will get silence.

2. **`check_obstruction()` also encodes the sign convention inconsistently with itself.**
   It flips `s`, `s_2`, `s_3` for `n < 0` but not `τ`/`ν` — which is *correct* given how
   the columns are stored (§3, Trap 2), but only by accident of the reader matching the
   writer. Add a comment there before someone "fixes" it in the wrong direction.

3. **`obstructs` is stale for any row whose invariants were filled in later.** It is only
   ever set to `True`, never back to `False`, and only at compute time. Recompute from
   the raw columns rather than trusting it. (For the current file it happens to agree
   exactly — 42/42.)

4. **Two rows have `verification = False`**: `id_num` 1442 (`n = 1`, over
   `18nh_09603071`) and 3455 (`n = −2`, over `18nh_04477958`). Neither has invariants
   computed. `rerun_all_unverified()` in `code/n_friends_search.py` exists to clear these.

---

## 8. Reproducing this analysis

```python
import csv, math
rows = list(csv.DictReader(open('data/n_friends(master version).csv')))

def obstructs(r):
    n = int(r['n']); m = abs(n); thr = m - math.sqrt(m); sgn = 1 if n > 0 else -1
    def f(x):
        try:
            v = float(x); return None if math.isnan(v) else v
        except ValueError:
            return None
    hits = []
    for k in ['s', 's_2', 's_3']:          # stored un-mirrored -> flip for n<0
        v = f(r[k])
        if v is not None and sgn * v > thr:
            hits.append(f'{k}={sgn*v:g}')
    for k in ['tau', 'nu']:                # stored ALREADY mirrored -> do NOT flip
        v = f(r[k])
        if v is not None and 2 * v > thr:
            hits.append(f'2*{k}={2*v:g}')
    return hits

hits = [(r, obstructs(r)) for r in rows]
hits = [(r, h) for r, h in hits if h]
print(len(hits), sum(1 for r, _ in hits if r['slice'] == '0'))   # -> 42 0
```
