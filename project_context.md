This code base is used in Dunfield and Gong's paper on using computers to tabulate
knots sliceness and ribbonness properties for the vast majority of knots up to 19 crossings.
The part of their paper that is relevant for what I want to do is in dunfield_and_gong_abridged.pdf.

Here is some more context for this project. The focus of this project is on the techniques for finding
0-friends of knots and how these can be used to show sliceness for knots. The ultimate goal will be
to generalize this to n-friends, where n is a positive integer. As an overview, the authors search
for 0-friends of a knot and then use pairs of 0-friends to compute RBG links (section 5.1). If the RBG link is
super-special, then we know the traces of K_G and K_B will be diffeomorphic (Theorem 5.8). This gives us a way
to categorize K_B and K_G as not smoothly slice when one of K_G or K_B have non-zero s-invariant (Theorem 5.9).

The ultimate goal of this project is to replicate Theorem 5.10 for n = 0 and for larger n. The first step
and first task is to be able to find n-friends of a given knot. Check the code/find_n_friends.py function
for the current implementation. To test this function, I have focused on the small example of the knot
6_2. The 1-surgery of 6_2 is isometric to the 1-surgery of K13n3596 (as is verified in is_isometric.py),
but the current code does not find the 1-surgery of 6_2 when searching. As is, this code is able to find
small examples of n-friends so it should be safe to assume the code is right.

The next goal for this project is to get the RBG search code working for n > 0. Some progress has
already been made and you can see the code for that in n_rbg.py, which generalizes much of the code
from rbg.py. To check the code, I am using the link diagram in 6_2_and_K13n3596_example which saves
a known RBG link from Figure 6 of Qin et al. This link diagram can be verified using the existing
code but the existing RBG generalization fails to actually find this link.
At the moment, the search runs for a long time and fails to find the link. This is
either because the link is being found but is failing one of the verification steps along the way
of the search (sprinkled throughout the n_blue_green_exteriors() are lots of verification steps to
make sure that the exteriors we find along the way are on the right track). Alternatively, it could
be that these verification steps are all correct and that we are just not searching long enough to
find this link. Your task is to run experiments and try to get this search working with the goal of
finding the afformentioned example.

More details about this can be found in dunfield_and_gong_abridged.pdf and Qin_et_al.pdf.

I'm using a conda environment to run python code in sage. The environment name is called sage.
Keep comments and docstrings concise and do not ramble excessively. Keep implementations simple
where possible. Prioritize code readability and clarity.