'''
This code uses 'find_n_friends.py' to search for n-friends in 'plausibly_unknown.csv'. To run this yourself,
use sage directly from the terminal and load the files based on their pathways
on your computer.
'''
import pandas as pd
import ast
import caffeine

#Change these file paths to match your own computer
load('/Users/henrigreamo/Desktop/plausibly_slice_v1/Code/find_n_friends.py')

out_file_path = "/Users/henrigreamo/Desktop/plausibly_slice_v1/Data/n_friends.csv"
in_file_path = "/Users/henrigreamo/Desktop/plausibly_slice_v1/Data/plausibly_unknown.csv"

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

def search(knot_name, knot_PD_code, knot_volume, knot, knot_ex, n, id_num, i, double_check=False):
    ans = find_common_n_surgery_via_words(knot_ex,n)
    if not ans:
        ans = []
    if len(ans)==0 and double_check:
        print(f"No {n}-friends found; trying again")
        ans = find_common_n_surgery_via_words(knot_ex,n)

    if len(ans)==0:
        print(f"No {n}-friends found")
    else:
        print(f"# of {n}-friends found: {len(ans)}")

    result = []
    for j in range(len(ans)):
        friend_ex = snappy.Manifold(ans[j][3])
        friend_knot = friend_ex.exterior_to_link()

        #I know this step seems weird but otherwise it fails the check common surgery
        friend_ex=friend_knot.exterior()

        verify = True
                
        if not check_common_surgery(knot_ex,friend_ex,n):
            print(f"Friend {j} failed surgery verification check")
            friend_knot = friend_knot.mirror()
            friend_ex = friend_knot.exterior()
            if not check_common_surgery(knot_ex,friend_ex,n):
                print(f"Mirror of friend {j} failed surgery verification check")
                verify = False
                
        friend_PD_code = friend_knot.PD_code()
        friend_num_crossings = len(friend_knot.crossings)
        friend_volume = float(friend_ex.volume())
                
        result.append({
            "id_num":id_num,
            "num_crossings":friend_num_crossings,
            "volume":friend_volume,
            "n":n,
            "verification":verify,
            "n_friend_name":knot_name,
            "n_friend_index": i,
            "knot_PD_code":str(friend_PD_code),
            "n_friend_PD_code":str(knot_PD_code)})
    return result

#Warning: this code reruns and overrides whatever data is at entry id_num
def rerun(id_num: int, n_friend_index: int, n: int, result_index: int = 0, double_check: bool = False):
    data_out = pd.read_csv(out_file_path)
    data_in = pd.read_csv(in_file_path)

    #Sometimes SageMath typecasts integers as its own special class of Ring Integers, so this code accounts for that
    id_num = int(id_num)
    i = int(n_friend_index)
    n = int(n)
    result_index = int(result_index)
    
    
    knot_name = data_in.at[i-1,"name"]
    knot_PD_code = ast.literal_eval(data_in.at[i-1,"PD_codes"])
    knot_volume = data_in.at[i-1,"volume"]
    knot = snappy.Link(knot_PD_code)
    knot_ex = knot.exterior()
    
    verified = data_out.at[id_num-1,"verification"]
    #print(data_out.loc[id_num-1])

    if verified:
        print(f"Line {id_num} has been verified as correct.")
        confirm = input("Confirm overriding data(Y/N): ")
        if confirm != "Y":
            return
        
    if knot_volume < 0:
        print("Knot is not hyperbolic")
        return
    print(f"Rerunning knot {i} with n={n} at ID_num={id_num}: " + str(knot_name))
    result = search(knot_name, knot_PD_code, knot_volume, knot, knot_ex, n, id_num, i, double_check)
    print()

    #Overrides the data at id_num
    if len(result) > 0:
        data_out.iloc[int(id_num) - 1] = pd.Series(result[result_index])
        data_out.to_csv(out_file_path, index=False)

def rerun_all_unverified(double_check=False):
    data_out = pd.read_csv(out_file_path)
    unverified = data_out.loc[data_out["verification"]==False]
    #print(unverified)
    for i in range(len(unverified)):
        row = unverified.iloc[i]
        id_num = row["id_num"]
        n_friend_index = row["n_friend_index"]
        n = row["n"]

        #Note that if this process would normally return multiple n-friends, then this will only show the first one
        rerun(id_num, n_friend_index,n,double_check=double_check)

def n_friends_search(start: int = 1,end: int = 1,max_n: int = 1, double_check=False):
    data_out = pd.read_csv(out_file_path)
    data_in = pd.read_csv(in_file_path)

    caffeine.on()

    #Warning: if there's only one entry the code will override it (it assumes the first row is empty)
    id_num = 1
    if len(data_out) > 1:
        id_num = data_out.at[len(data_out) - 1,"id_num"] + 1
    else:
        data_out=pd.DataFrame([])

    new_data=[]

    for i in range (start,end + 1):
        '''
        For now, I'm going to keep these non-slice knots as a test for the program
        sliceness = int(data_in.at[i-1,"slice"])
        if sliceness == -1:
            print("Knot is not slice")
            continue
        '''
        
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
            
            result = []
            result = search(knot_name, knot_PD_code, knot_volume, knot, knot_ex, n, id_num, i, double_check=False)
            new_data = new_data + result
            id_num += len(result)
            print()

        #This might be slow, but it regularly updates the output data
        new_rows = pd.DataFrame(new_data)
        out = pd.concat([data_out,new_rows], ignore_index=True)
        out.to_csv(out_file_path, index=False)
    caffeine.off()
                
            
        
    
          






