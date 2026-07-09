'''
This test file is used by me to test the various scripts. To run this yourself,
use sage directly from the terminal and load the files based on their pathways
on your computer.
'''

#load('/Users/henrigreamo/Desktop/plausibly_slice_v1/rbg.py')
load('/Users/henrigreamo/Desktop/plausibly_slice_v1/n-rbg.py')
load('/Users/henrigreamo/Desktop/plausibly_slice_v1/find_0_friends.py')

def check_isometric(E1,E2):
    if abs(E1.volume()-E2.volume())>0.0001:
        return false
    return E1.is_isometric_to(E2)
    

'''
#Test n-RBG link is_n_super_special function for n=0:
blue_ex=snappy.Manifold('K11n34')
ans=find_common_zero_surgery_via_words(blue_ex,3)
green_ex=snappy.Manifold(ans[0][3])
E1=blue_ex.copy()
E1.dehn_fill((0,1))
E2=green_ex.copy()
E2.dehn_fill((0,1))
print("Knots share common 0-surgery:" + str(check_isometric(E1,E2)))
forms_super_special_NRBG_link(blue_ex,green_ex)
'''


#Test modified n-RBG link is_n_super_special function
n = 3
blue_ex=snappy.Manifold('6_2')
green_ex=snappy.Manifold('K14n10164')
E1=blue_ex.copy()
E1.dehn_fill((n,1))
E2=green_ex.copy()
E2.dehn_fill((n,1))
print("Knots share common n-surgery:" + str(check_isometric(E1,E2)))
forms_super_special_NRBG_link(blue_ex,green_ex)








