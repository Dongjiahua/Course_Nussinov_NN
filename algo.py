import numpy as np
import pandas as pd



from collections import OrderedDict

def slice(ret,end):
    left = ret[end:]
    for l in left:
        ret.remove(l)

class Nussinov:
    def __init__(self, mode="normal", max_att=100):
        if mode == "normal":
            self.max_att = max_att
            self.tau = {"A":"U", "U":"A", "C":"G", "G":"C"}
            self.rna = "GGUCCAC"
            sc, self.matrix = self.solve(self.rna)
            self.stru, self.trace = self.backward(self.rna, self.matrix)

            
            

    def dfs(self,rna, matrix, array, optimals, i, j, steps):
        # print(f"i:{i},j:{j}, steps:{steps},array:{''.join(array)}")
        if j<=i:
            if i==j:
                array[i]="."
                steps+=1
            if matrix[i,j]==0 and steps == len(array):
                optimals.add("".join(array))
            return
        if len(optimals)>self.max_att:
            return
        if self.tau[rna[i]]==rna[j] and i<len(rna)-1 and j>0 and matrix[i,j] == matrix[i+1,j-1]+1:
            array[i] = "("
            array[j] = ")"
            self.dfs(rna, matrix, array, optimals, i+1,j-1, steps+2)
        if j>0 and matrix[i,j]==matrix[i,j-1]:
            array[j] = "."
            self.dfs(rna, matrix, array, optimals, i,j-1, steps+1)
        for k in range(i+1,j):
            if matrix[i,j]==matrix[i,k]+matrix[k+1,j]:
                self.dfs(rna, matrix, array, optimals,i,k, steps)
                self.dfs(rna, matrix, array, optimals,k+1,j, steps+k-i+1)
                break
            
    def solve(self,rna):
        rna = rna.upper()
        # print(rna)
        length = len(rna)
        matrix = np.zeros((length,length))
        for l in range(1,length):
            for i in range(length-l):
                j=i+l
                cur_max = 0
                if self.tau[rna[i]]==rna[j]:
                    cur_max = max(cur_max, matrix[i+1,j-1]+1)
                cur_max = max(cur_max, matrix[i,j-1])
                for k in range(i+1,j):
                    cur_max = max(cur_max, matrix[i,k]+matrix[k+1,j])

                matrix[i,j] = cur_max 
        
        return matrix[0,length-1], matrix 

    def dfs_backtrace(self, i, j, rna, matrix,ret,structure):
        length = len(ret)
        success = False
        if i>=j:
            if matrix[i,j]==0:
                if i==j:
                    ret.append([i,j])
                success = True
            else:
                return False

        if not success and self.tau[rna[i]]==rna[j] and matrix[i,j] == matrix[i+1,j-1]+1 and structure[i]=="(" and structure[j]==")":
            ret.append([i+1,j-1])
            success = self.dfs_backtrace(i+1,j-1, rna,matrix,ret,structure)
        if not success and matrix[i,j]==matrix[i,j-1] and structure[j]==".":
            slice(ret,length)
            ret.append([i,j-1])
            success = self.dfs_backtrace(i,j-1, rna,matrix,ret,structure)  
        elif not success:
            for k in range(i+1,j-1):
                if matrix[i,j]==matrix[i,k]+matrix[k+1,j]: 
                    slice(ret,length)
                    ret.append([i,k])
                    ret.append([k+1,j])  
                    s1 = self.dfs_backtrace(i,k,rna,matrix,ret,structure) 
                    s2 = self.dfs_backtrace(k+1,j,rna,matrix,ret,structure)
                    success = s1 and s2 
                    if success:
                        break    
        if not success:
            slice(ret,length)
        return success    



    def backward(self,rna, matrix):
        rna = rna.upper()
        array = [""]*len(rna)
        optimals = set()
        self.dfs(rna, matrix, array, optimals, 0, len(rna)-1, 0, )
        trace = [self.backtrace(rna,opt,matrix) for opt in optimals]
        optimals = list(optimals)
        trace = {optimals[i]:trace[i] for i in range(len(optimals))}        
        return optimals, trace

    def backtrace(self, rna, structure, matrix):
        n = len(structure)
        j = n-1
        i = 0
        ret = [[0,j]]
        self.dfs_backtrace(i,j,rna,matrix,ret, structure)
        return ret
        


data = OrderedDict(
    [
        ("Date", ["2015-01-01", "2015-10-24", "2016-05-10", "2017-01-10", "2018-05-10", "2018-08-15"]),
        ("Region", ["Montreal", "Toronto", "New York City", "Miami", "San Francisco", "London"]),
        ("Temperature", [1, -20, 3.512, 4, 10423, -441.2]),
        ("Humidity", [10, 20, 30, 40, 50, 60]),
        ("Pressure", [2, 10924, 3912, -10, 3591.2, 15]),
    ]
)
df = pd.DataFrame(np.array(np.eye(4)).tolist())

if __name__ == "__main__":
    nu = Nussinov()
    RNA = "GGUCCACGGUGGC"
    stru="(..)((().))()"
    s, matrix = nu.solve(RNA)
    # structure, ret = nu.backward(RNA,matrix)
    ret=  nu.backtrace(RNA,stru,matrix)
    print(ret)

