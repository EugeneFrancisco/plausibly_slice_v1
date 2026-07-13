'''
This test file is used by me to test the various scripts. To run this yourself,
use sage directly from the terminal and load the files based on their pathways
on your computer.
'''

#load('/Users/henrigreamo/Desktop/plausibly_slice_v1/rbg.py')
load('/Users/henrigreamo/Desktop/plausibly_slice_v1/n-rbg.py')
#load('/Users/henrigreamo/Desktop/plausibly_slice_v1/n_rbg_Eugene_version.py')
load('/Users/henrigreamo/Desktop/plausibly_slice_v1/find_0_friends.py')

def check_common_surgery(E1,E2,n):
    S1=E1.copy()
    S1.dehn_fill(n_surgery_slope(S1,n),0)
    S2=E2.copy()
    S2.dehn_fill(n_surgery_slope(S2,n),0)
    print("Volume difference: " + str(abs(S1.volume()-S2.volume())))
    if abs(S1.volume()-S2.volume())>0.0001:
        return false
    print("Is isometric: " + str(S1.is_isometric_to(S2)))
    return S1.is_isometric_to(S2)
    


#Test n-RBG link is_n_super_special function for n=0:
def zero_test():
    blue_ex=snappy.Manifold('K11n34')
    ans=find_common_zero_surgery_via_words(blue_ex,3)
    green_ex=snappy.Manifold(ans[0][3])
    print("Knots share common 0-surgery:" + str(check_common_surgery(blue_ex,green_ex,n)))
    return forms_super_special_NRBG_link(0,blue_ex,green_ex)



#Test modified n-RBG link is_n_super_special function
def n_test(k):
    examples=[(1,snappy.Manifold('6_2'),snappy.Manifold('K13n3596')),
              (3,snappy.Manifold('6_2'),snappy.Manifold('K14n10164')),
              (3,snappy.Manifold('6_3'),snappy.Manifold('K14n15962')),
              (3,snappy.Manifold('10_125'),snappy.Manifold('10_132'))  
        ]
    n = examples[k][0]
    blue_ex=examples[k][1]
    green_ex=examples[k][2]

    
    print(f"Knots share common {n}-surgery:" + str(check_common_surgery(blue_ex,green_ex,n)))
    return forms_super_special_NRBG_link(n,blue_ex,green_ex)







