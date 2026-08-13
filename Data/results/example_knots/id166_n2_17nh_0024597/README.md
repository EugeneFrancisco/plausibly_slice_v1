# Example pair: id_num 166 (n = 2)

- Source CSV: `data/n_friends(master version).csv` row `id_num = 166`
- Companion row: `Data/results/obstructed_pairs.csv`
- Census knot K: `17nh_0024597` (index 226 in `plausibly_unknown.csv`)
- Original K crossings: 17; K τ = 1; K s₀ = 2
- n-friend K' crossings: 35
- K' invariants: τ = 1.0, ν = 1.0, s = 2.0, s₂ = 2.0, s₃ = 2.0
- Obstruction: s=2; s_2=2; s_3=2; 2*tau=2; 2*nu=2
- RBG framings (r, b, g): [(0, 1), (0, 1), (0, 1)]

## Files
- `K_census.pd.txt` — PD code of K (the census knot).
- `K_friend.pd.txt` — PD code of K' (its n-friend).
- `RBG_link.pd.txt` — PD code of the n-special RBG link found.
- `.lnk` projections were not written (headless / plink unavailable). Re-open the PD codes with `snappy.Link(pd).view().save_file(...)` in an interactive session to generate them.
