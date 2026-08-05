'''
This code uses 'find_n_friends.py' to search for n-friends in 'plausibly_unknown.csv'. To run this yourself,
use sage directly from the terminal and load the files based on their pathways
on your computer.
'''
import pandas as pd
import ast
import caffeine
from tqdm import tqdm

#Change these file paths to match your own computer
load('/Users/henrigreamo/Desktop/plausibly_slice_v1/Code/find_n_friends.py')

out_file_path = "/Users/henrigreamo/Desktop/plausibly_slice_v1/Data/n_friends.csv"
in_file_path = "/Users/henrigreamo/Desktop/plausibly_slice_v1/Data/plausibly_unknown.csv"

data_types = {
    "id_num": int,
    "num_crossings": int,
    "volume": float,
    "n": int,
    "verification": bool,
    "s": str,
    "s_2": str,
    "s_3": str,
    "nu": str,
    "tau": str,
    "obstructs": bool,
    "n_friend_name": str,
    "n_friend_index": int,
    "knot_PD_code": str,
    "n_friend_PD_code": str
    }

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

#Tries to simplify a knot with so many iterations
def simplify_knot(knot, iterations: int = 5, type_3_limit: int = 1000):
    K = knot.copy()
    for i in tqdm(range(iterations), desc="Attempting to simplify knot"):
        reduction = K.simplify('global',type_3_limit)
    crossings_reduction = len(knot.crossings)-len(K.crossings)
    print("Reduced crossings by: " + str(crossings_reduction))
    return K

#Keeps iterating until the knot can no longer be simplified
def super_simplify(knot, iterations: int = 5, type_3_limit: int = 5000):
    K= knot.copy()
    for i in range(iterations):
        reduced = True
        while reduced:
            reduced = False
            reduced = K.simplify('global',type_3_limit)

    crossings_reduction = len(knot.crossings)-len(K.crossings)
    print("Reduced crossings by: " + str(crossings_reduction))
    return K

#Simplifies the table entries as much as possible (adjust parameters for less/more simplification)
#Note that this doesn't reverify anything (the knots should still be the same knot)
def simplify_table(start: int, end: int,iterations: int = 5, type_3_limit: int = 1000):
    data = pd.read_csv(out_file_path,dtype=data_types)

    caffeine.on()
    for i in range(start,end+1):
        print(f"Simplifying knot {i}")
        #print(data.iloc[i-1])
        PD_code = ast.literal_eval(data.at[i-1,"knot_PD_code"])
        K=snappy.Link(PD_code)
        K=simplify_knot(K, iterations, type_3_limit)
        data.at[i-1,"knot_PD_code"]=str(K.PD_code())
        data.at[i-1,"num_crossings"]=int(len(K.crossings))
        data.to_csv(out_file_path, index=False)
    caffeine.off()

def simplify_specified_knots(indices: list, iterations: int = 5, type_3_limit: int = 1000):
    data = pd.read_csv(out_file_path,dtype=data_types)

    caffeine.on()
    for i in indices:
        i = int(i)
        print(f"Simplifying knot {i}")
        #print(data.iloc[i-1])
        PD_code = ast.literal_eval(data.at[i-1,"knot_PD_code"])
        K=snappy.Link(PD_code)
        K=simplify_knot(K, iterations, type_3_limit)
        data.at[i-1,"knot_PD_code"]=str(K.PD_code())
        data.at[i-1,"num_crossings"]=int(len(K.crossings))
        data.to_csv(out_file_path, index=False)
    caffeine.off()

def super_simplify_table(start: int, end: int,iterations: int = 5, type_3_limit: int = 1000, max_crossings: int = 1000):
    data = pd.read_csv(out_file_path,dtype=data_types)

    caffeine.on()
    for i in range(start,end+1):
        #print(data.iloc[i-1])
        PD_code = ast.literal_eval(data.at[i-1,"knot_PD_code"])
        K=snappy.Link(PD_code)
        if (len(K.crossings) > max_crossings):
            continue
        print(f"Simplifying knot {i}")
        K=super_simplify(K, iterations, type_3_limit)
        data.at[i-1,"knot_PD_code"]=str(K.PD_code())
        data.at[i-1,"num_crossings"]=int(len(K.crossings))
        data.to_csv(out_file_path, index=False)
    caffeine.off()

def super_simplify_specified(indices: list, iterations: int= 5, type_3_limit: int = 5000):
    data = pd.read_csv(out_file_path,dtype=data_types)

    caffeine.on()
    for i in indices:
        i = int(i)
        print(f"Simplifying knot {i}")
        #print(data.iloc[i-1])
        PD_code = ast.literal_eval(data.at[i-1,"knot_PD_code"])
        K=snappy.Link(PD_code)
        K=super_simplify(K, iterations, type_3_limit)
        data.at[i-1,"knot_PD_code"]=str(K.PD_code())
        data.at[i-1,"num_crossings"]=int(len(K.crossings))
        data.to_csv(out_file_path, index=False)
    caffeine.off()
    
#Searchs for n-friends of knot and returns all the relevant information about each found n-friend
def search(knot_name, knot_PD_code, knot_volume, knot, knot_ex, n, id_num, i, knot_sliceness,double_check=False):
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
            "id_num":(id_num + j),
            "num_crossings":friend_num_crossings,
            "volume":friend_volume,
            "n":n,
            "verification":verify,
            "s": "",
            "s_2": "",
            "s_3": "",
            "nu":"",
            "tau": "",
            "obstructs": False,
            "n_friend_name":knot_name,
            "n_friend_index": i,
            "slice": knot_sliceness,
            "knot_PD_code":str(friend_PD_code),
            "n_friend_PD_code":str(knot_PD_code)})
    return result

#This code reruns and overrides whatever data is at entry id_num
def rerun(id_num: int, n_friend_index: int, n: int, result_index: int = 0, double_check: bool = False):
    data_out = pd.read_csv(out_file_path,dtype=data_types)
    data_in = pd.read_csv(in_file_path)

    #Sometimes SageMath typecasts integers as its own special class of Ring Integers, so this code accounts for that
    id_num = int(id_num)
    i = int(n_friend_index)
    n = int(n)
    result_index = int(result_index)
    
    
    knot_name = data_in.at[i-1,"name"]
    knot_PD_code = ast.literal_eval(data_in.at[i-1,"PD_codes"])
    knot_volume = data_in.at[i-1,"volume"]
    knot_sliceness = data_in.at[i-1,"slice"]
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
    result = search(knot_name, knot_PD_code, knot_volume, knot, knot_ex, n, id_num, i, knot_sliceness,double_check)
    print()

    #Makes sure the result has the correct id_num (it can be off sometimes)
    result[result_index]["id_num"]=id_num

    #Overrides the data at id_num
    if len(result) > 0:
        data_out.iloc[int(id_num) - 1] = pd.Series(result[result_index])
        data_out.to_csv(out_file_path, index=False)

#Reruns and overrides the data at each row which has "False" under verification
def rerun_all_unverified(double_check=False):
    data_out = pd.read_csv(out_file_path,dtype=data_types)
    unverified = data_out.loc[data_out["verification"]==False]
    #print(unverified)
    for i in range(len(unverified)):
        row = unverified.iloc[i]
        id_num = row["id_num"]
        n_friend_index = row["n_friend_index"]
        n = row["n"]

        #Note that if this process would normally return multiple n-friends, then this will only show the first one
        rerun(id_num, n_friend_index,n,double_check=double_check)

#Searches for n-friends of the knots with index in [start,end]
def n_friends_search(start: int = 1,end: int = 1,min_n: int = 1, max_n: int = 5, double_check=False):
    data_out = pd.read_csv(out_file_path,dtype=data_types)
    data_in = pd.read_csv(in_file_path)

    #Turns on caffeine to suppress macOS sleep functions
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

        knot_sliceness = data_in.at[i-1,"slice"]
        
        if knot_volume < 0:
            print("Knot is not hyperbolic")
            continue
        for n in range(min_n, max_n+1):
            print(f"Checking knot {i} with n={n}: " + str(knot_name))
            
            result = []
            result = search(knot_name, knot_PD_code, knot_volume, knot, knot_ex, n, id_num, i, knot_sliceness, double_check=False)
            new_data = new_data + result
            id_num += len(result)
            print()

        #Regularly updates the output data
        new_rows = pd.DataFrame(new_data)
        out = pd.concat([data_out,new_rows], ignore_index=True)
        out.to_csv(out_file_path, index=False)
    caffeine.off()

#Can be reworked to add any sort of data that you might want to transfer
def add_sliceness():
    data_out = pd.read_csv(out_file_path,dtype=data_types)
    data_in = pd.read_csv(in_file_path)

    for i in range(len(data_out)):
        index = int(data_out.at[i,"n_friend_index"])
        sliceness = data_in.at[index-1,"slice"]
        data_out.at[i,"slice"]=sliceness
    data_out.to_csv(out_file_path, index=False)
        
    
          






