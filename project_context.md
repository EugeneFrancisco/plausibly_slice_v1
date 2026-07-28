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
and first task was to be able to find n-friends of a given knot. Check the find_n_friends.py function
for the current implementation. At this point, this code seems like it generally works.

The next part of this project is to get the RBG search code working for n > 0. Some progress has
already been made and you can see the code for that in n-rbg.py, which generalizes much of the code
from rbg.py. To check the code, I am using the link diagram in 6_2_and_K13n3596_example which saves
a known RBG link from Figure 6 of Qin et al. This link diagram can be verified using the existing
code but the existing RBG generalization fails to actually find this link.

More details about this can be found in dunfield_and_gong_abridged.pdf and Qin_et_al.pdf.

I'm using a conda environment to run python code in sage. The environment name is called sage.
Keep comments and docstrings concise and do not ramble excessively. Keep implementations simple
where possible. Prioritize code readability and clarity.

The goal now of this project is to start using the code to try and conclude something about the
knots for which Dunfield and Gong could not conclude slice status. There were some 11,000 knots
for which this was the case. 