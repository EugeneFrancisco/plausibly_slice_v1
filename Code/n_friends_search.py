'''
This code uses 'find_n_friends.py' to search for n-friends in 'plausibly_unknown.csv'. To run this yourself,
use sage directly from the terminal and load the files based on their pathways
on your computer.
'''
import pandas as pd
import ast
load('/Users/henrigreamo/Desktop/plausibly_slice_v1/Code/find_n_friends.py')

#Checks that the manifolds have the same n-surgery
def check_common_surgery(E1,E2,n):
    S1=E1.copy()
    S1.dehn_fill(n_surgery_slope(S1,n),0)
    S2=E2.copy()
    S2.dehn_fill(n_surgery_slope(S2,n),0)
    return isometric(S1,S2)

#Safer check that the manifolds are isometric
def isometric(E1,E2):
    A=snappy.Manifold(E1)
    B=snappy.Manifold(E2)

    #Checks that the hyperbolic volumes match
    if abs(A.volume()-B.volume())>0.000001:
        return False
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
        return None

#Checks all combinations of mirroring the knots
def check_mirrors(L1,L2,n):
    for i in range (0,4):
        E1=L1.exterior()
        E2=L2.exterior()
        if (i > 0 and i < 3):
            E1=L1.mirror().exterior()
        if (i > 1):
            E2=L2.mirror().exterior()
        print("Check #" + str(i) + ": " + str(check_common_surgery(E1,E2,n)))

def n_friends_search(start: int = 1,end: int = 1,max_n: int = 1, double_check=False):
    out_file_path = "/Users/henrigreamo/Desktop/plausibly_slice_v1/Data/n_friends.csv"
    in_file_path = "/Users/henrigreamo/Desktop/plausibly_slice_v1/Data/plausibly_unknown.csv"
    data_out = pd.read_csv(out_file_path)
    data_in = pd.read_csv(in_file_path)

    #Warning: if there's only one entry the code will override it (it assumes the first row is empty)
    id_num = 1
    if len(data_out) > 1:
        id_num = data_out.at[len(data_out) - 1,"id_num"] + 1
    else:
        data_out=pd.DataFrame([])

    new_data=[]
    
    #df.loc[len(df)]=[1,5,2.6,3,2,"K14n001",[[1,0]],[[2,1]]]
    #df.to_csv(file_path, index=False)
    for i in range (start,end + 1):
        
        
        knot_name = data_in.at[i-1,"name"]
        knot_PD_code = ast.literal_eval(data_in.at[i-1,"PD_codes"])
        knot_volume = data_in.at[i-1,"volume"]
        knot = snappy.Link(knot_PD_code)
        knot_ex = knot.exterior()
        
        if knot_volume < 0:
            print("Knot is not hyperbolic")
            continue
        for n in range(1, max_n+1):
            print(f"Checking knot {i} with n={n}: " + str(knot_name))
            ans = find_common_n_surgery_via_words(knot_ex,n)
            if len(ans)==0 and double_check:
                print(f"No {n}-friends found; trying again")
                ans = find_common_n_surgery_via_words(knot_ex,n)

            if len(ans)==0:
                print(f"No {n}-friends found")
            else:
                print(f"# of {n}-friends found: {len(ans)}")
            
            for j in range(len(ans)):
                #I know this step seems weird but otherwise it fails the check common surgery
                friend_ex = snappy.Manifold(ans[j][3])
                friend_knot = friend_ex.exterior_to_link()
                friend_ex=friend_knot.exterior()

                verify = True
                
                if not check_common_surgery(knot_ex,friend_ex,n):
                    print(f"Friend {j} failed surgery verification check")
                    friend_knot = knot.mirror()
                    friend_ex = friend_knot.exterior()
                    if not check_common_surgery(knot_ex,friend_ex,n):
                        print(f"Mirror of friend {j} failed surgery verification check")
                        verify = False
                
                friend_PD_code = friend_knot.PD_code()
                friend_num_crossings = len(friend_knot.crossings)
                friend_volume = float(friend_ex.volume())
                
                new_data.append({
                    "id_num":id_num,
                    "num_crossings":friend_num_crossings,
                    "volume":friend_volume,
                    "n":n,
                    "verification":verify,
                    "n_friend_name":knot_name,
                    "n_friend_index": i,
                    "knot_PD_code":friend_PD_code,
                    "n_friend_PD_code":knot_PD_code})
                id_num += 1
            print()

        #This might be slow, but it regularly updates the output data
        new_rows = pd.DataFrame(new_data)
        out = pd.concat([data_out,new_rows], ignore_index=True)
        out.to_csv(out_file_path, index=False)
                
            
        
    
          






