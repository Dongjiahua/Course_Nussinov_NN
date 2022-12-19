import tempfile
import subprocess 
import random
from algo import Nussinov
import numpy as np
import matplotlib.pyplot as plt
import math
nu = Nussinov(max_att=5)
times =5
def generate(n):
    rna_list = [random.choice(["U","A","G","C"]) for i in range(n)]
    ran_data = "".join(rna_list)
    return ran_data

def get_score(seq):
    return seq.count("(")

def test():
    mxs = []
    nus = []

    diffs=[]
    for length in range(5,100,5):
        mx_score = []
        nu_score = []
        diff=[]
        print(f"process {length}")
        for ti in range(times):
            rna = generate(length)
            sc, nu.matrix = nu.solve(rna)
            nu.stru, nu.trace = nu.backward(rna, nu.matrix)
            nu_st = nu.stru[0]
            with tempfile.NamedTemporaryFile(suffix='.fa') as tempf:
                tempf.write(f">S\n{rna}".encode())
                tempf.seek(0)
                output = subprocess.check_output(f'mxfold2 predict {tempf.name}', shell=True).decode()
                mx_st = output.split("\n")[2].split(" ")[0]
            
            mx_score.append(get_score(mx_st))
            nu_score.append(get_score(nu_st))
            d = [1 for s in range(length) if mx_st[s]!=nu_st[s]]
            diff.append(math.fabs(sum(d)/length))
        mxs.append(sum(mx_score)/times)
        nus.append(sum(nu_score)/times)
        diffs.append(sum(diff)/times)
    xs = np.array(list(range(5,100,5)))
    # plt.xlabel('length')
    # plt.ylabel('pairs')
    # plt.title('xlabels() function')
    # plt.plot(xs, np.array(nus),label='Nussinov')
    # plt.plot(xs, np.array(mxs),label='MXFold2')
    # plt.legend()
    # plt.show()

    plt.xlabel('length')
    plt.ylabel('L1 difference')
    plt.plot(xs, np.array(diffs))
    plt.legend()
    plt.show()

if __name__ == "__main__":
    test()

        