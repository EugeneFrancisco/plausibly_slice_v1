'''
This test file is used to test the various scripts. To run this yourself,
use sage directly from the terminal and load the files based on their pathways
on your computer.
'''

load('/Users/henrigreamo/Desktop/plausibly_slice_v1/Code/n_rbg.py')
load('/Users/henrigreamo/Desktop/plausibly_slice_v1/Code/find_n_friends.py')

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

def search_test():
    pd_code = [[4,2,5,1],[2,9,3,10],[8,3,9,4],[10,6,11,5],
               [6,15,7,16],[14,7,15,8],[11,21,12,20],
               [25,13,26,12],[13,18,14,19],[23,17,24,16],
               [17,23,18,22],[19,1,20,26],[21,24,22,25]]
    K=snappy.Link(pd_code)
    n = 1
    
    E = snappy.ManifoldHP(K.exterior())
    M = E.copy()
    # the n-surgery Z_K^{(n)}
    M.dehn_fill(n_surgery_slope(E, n))

    # drill_word needs solution_type 1
    M = _geometric_triangulation(M)           
    if M is None:
        # non-hyperbolic or no geometric
        print("Failed to find a geometric triangulation.") 
        return                                

    # triangulation reachable
    G = M.fundamental_group(False, False, False)
    phi = nsagetools.MapToAbelianization(G)

    geodesics = safe_length_spectrum(M, 3.0)
    if geodesics is None:
        print("No geodesics found.")
        return

    print(f"Testing {len(geodesics)} geodesics for {n} friends.")
    for g in tqdm(geodesics, desc=f"Testing geodesics for {n} friends"):
        if g.length.real() < 0.0:
            continue
        if not _word_is_homology_generator(phi, g.word):
            continue
        
        F = None
        try:
            F = M.drill_word(g.word).filled_triangulation()
        except:
            print("Error: failed to drill curve")
            continue
        if F.solution_type(enum=True) not in [1, 2]:
            continue
        if safe_is_isometric_to(E, F):        # the geodesic recovers K itself
            continue
        slope = is_knot_exterior(F)
        if not slope:
            continue
        
        # Reframe F as the exterior of K' (meridian = its S^3 slope).
        F.dehn_fill(slope)
        F.set_peripheral_curves('fillings')
        F.dehn_fill((0, 0))
        if n != 0 and not _n_surgery_recovers(F, n, M):
            continue                          # rational-surgery coincidence
        N = snappy.ManifoldHP(F)
        slope = n_surgery_slope(N,n)
        print("Drilled manifold slope: " + str(slope))
        print(check_common_surgery(N,M,n))






