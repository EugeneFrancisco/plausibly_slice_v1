# Preliminary-result diagrams

The `friends` directory contains four knot diagrams and their SnapPea
triangulations:

| File stem | Result | Diagram crossings |
|---|---|---:|
| `6_2_1_friend_K13n3596` | 1-friend of `6_2`, identified as `K13n3596` | 13 |
| `6_2_3_friend_K14n10164_mirror` | 3-friend of `6_2`, the mirror of `K14n10164` | 14 |
| `conway_1_friend_A` | first 1-friend of Conway (`K11n34`) | 34 |
| `conway_1_friend_B` | second 1-friend of Conway (`K11n34`) | 40 |

The `rbg_links` directory contains two 1-special RBG diagrams, with components
ordered red, blue, green and slopes
`[(-1, 1), (0, 1), (0, 1)]`:

| File stem | Associated pair | Diagram crossings |
|---|---|---:|
| `6_2_K13n3596_1_RBG` | `6_2` and `K13n3596` | 11 |
| `conway_friend_A_1_RBG` | Conway and `conway_1_friend_A` | 18 |

The `.lnk` files are SnapPy Link Editor projections. Open them from the Link
Editor, or load them without the GUI as follows:

```python
import snappy
from plink import LinkManager

manager = LinkManager()
with open("preliminary_results/friends/conway_1_friend_A.lnk") as stream:
    manager._from_string(stream.read())
link = snappy.Link(manager.PD_code())
```

The matching `.tri` files load directly:

```python
M = snappy.Manifold("preliminary_results/friends/conway_1_friend_A.tri")
```

The `png` directory contains a consistently styled 1600-by-1600 PNG for each
projection and a labeled contact sheet, `all_six_diagrams.png`. Knot diagrams
are charcoal; the RBG diagrams use red, blue, and green for their components
in that order. Regenerate the images with:

```bash
PYTHONPATH=preliminary_results sage -python preliminary_results/render_pngs.py
```

From the repository root, verify every saved file with:

```bash
PYTHONPATH=. sage -c "from preliminary_results.verify_results import main; main()"
```

`generate_friend_diagrams.py` records the search isosigs and regenerates
short, correctly oriented friend diagrams. The RBG files are the actual links
returned by the successful searches, rather than copies of the target pairs.
