import numpy as np
import scipy as sp
import pandas as df

from scipy.sparse.linalg import splu

A_sp = df.read_csv("A_Sparse.txt")
A_row = A_sp["Row(int)"].to_numpy(dtype=int) - 1
A_col = A_sp["Col(int)"].to_numpy(dtype=int) - 1
A_data = A_sp["Data(float)"].to_numpy(dtype=float)

A_sp = sp.sparse.csc_matrix((A_data,(A_row,A_col)))

b_sp = df.read_csv("b_Sparse.txt")
b_row = b_sp["Row(int)"].to_numpy(dtype=int) - 1
b_col = b_sp["Col(int)"].to_numpy(dtype=int) - 1
b_data = b_sp["Data(float)"].to_numpy(dtype=float)

b_sp = sp.sparse.csc_matrix((b_data,(b_row,b_col)),shape=(A_sp.shape[1],max(b_col)+1)).toarray()

A_LU = splu(A_sp)
x = A_LU.solve(b_sp)

v_points = df.read_csv("v_Points.txt").to_numpy()
N = df.read_csv("info.txt",header=None).to_numpy()[0]
M = df.read_csv("info.txt",header=None).to_numpy()[1]

Dv = np.zeros(v_points.shape[0])
for k in range(0,v_points.shape[0]) :
    v1_i: int = v_points[k][0]
    v2_i: int = v_points[k][1]
    v1_j: int = v_points[k][2]
    v2_j: int = v_points[k][3]
    row_v1: int = (v1_i + ((N) * (v1_j)))[0]
    row_v2: int = (v2_i + ((N) * (v2_j)))[0]
    Dv[k] = (x[row_v2,k] - x[row_v1,k])
np.savetxt("Dv.txt",Dv,delimiter='\n')

np.savetxt("x.txt",x,delimiter='\n',fmt=(x.shape[1]*'%lf,').rstrip(','))