'''
This test file is used by me to test the various scripts. To run this yourself,
use sage directly from the terminal and load the files based on their pathways
on your computer.
'''

#load('/Users/henrigreamo/Desktop/plausibly_slice_v1/rbg.py')
#load('/Users/henrigreamo/Desktop/plausibly_slice_v1/n-rbg.py')
load('/Users/henrigreamo/Desktop/plausibly_slice_v1/n-rbg.py')
load('/Users/henrigreamo/Desktop/plausibly_slice_v1/find_n_friends.py')

#Having the same hyperbolic volume doesn't prevent the manifolds from being the same
def check_common_surgery(E1,E2,n):
    S1=E1.copy()
    S1.dehn_fill(n_surgery_slope(S1,n),0)
    S2=E2.copy()
    S2.dehn_fill(n_surgery_slope(S2,n),0)
    A=snappy.ManifoldHP(S1)
    B=snappy.ManifoldHP(S2)
    print("Volume difference check: " + str(abs(A.volume()-B.volume())<0.0001))
    try:
        is_isometric = A.is_isometric_to(B)
        return is_isometric
    except:
        #print(f"Struggling to determine if {n}-surgeries are isometric")
        for i in range(1,10):
            A.randomize()
            B.randomize()
            try :
                return A.is_isometric_to(B)
            except:
                continue
        print(f"Could not determine if {n}-surgeries are isometric")

def check_mirrors(L1,L2,n):
    for i in range (0,4):
        E1=L1.exterior()
        E2=L2.exterior()
        if (i > 0 and i < 3):
            E1=L1.mirror().exterior()
        if (i > 1):
            E2=L2.mirror().exterior()
        print("Check #" + str(i) + ": " + str(check_common_surgery(E1,E2,n)))


#Test n-RBG link is_n_super_special function for n=0:
def zero_test():
    blue_ex=snappy.Manifold('K11n34')
    ans=find_common_n_surgery_via_words(blue_ex,0)
    green_ex=snappy.Manifold(ans[0][3])
    print("K_G identified as: " + str(green_ex.identify()))
    print("Knots share common 0-surgery:" + str(check_common_surgery(blue_ex,green_ex,0)))
    return forms_special_NRBG_link(0,blue_ex,green_ex)



#Test modified n-RBG link is_n_super_special function
def n_test(k):
    #Only example k=0 seems to work properly
    examples=[(1,snappy.Link('6_2'),snappy.Link('K13n3596')),
              (3,snappy.Link('6_2'),snappy.Link('K14n10164').mirror()),
              (3,snappy.Link('6_3'),snappy.Link('K14n15962').mirror()),
              (3,snappy.Link('10_125'),snappy.Link('10_132').mirror())]
    n = examples[k][0]
    blue_ex=examples[k][1].exterior()
    green_ex=examples[k][2].exterior()

    #The n-friends search needs to be refined to account for mirrors
    #ans=find_common_n_surgery_via_words(blue_ex,n)
    #if ans is None:
    #   print(f"Could not find any {n}-friends.")
    #   return
    #green_ex=snappy.Manifold(ans[0][3])
    
    print("K_G identified as: " + str(green_ex.identify()))

    
    print(f"Knots share common {n}-surgery:" + str(check_common_surgery(blue_ex,green_ex,n)))
    return forms_special_NRBG_link(n,blue_ex,green_ex)






