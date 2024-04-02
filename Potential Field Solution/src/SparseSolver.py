import numpy as np
import scipy as sp
import pandas as df

from pypardiso import spsolve

A_sp = df.read_csv("A_Sparse.txt")
A_row = A_sp["Row(int)"].to_numpy(dtype=int) - 1
A_col = A_sp["Col(int)"].to_numpy(dtype=int) - 1
A_data = A_sp["Data(float)"].to_numpy(dtype=float)

A_sp = sp.sparse.csr_matrix((A_data,(A_row,A_col)))

b = np.fromfile("b.txt",sep='\n')

x = spsolve(A_sp, b)
np.savetxt("x.txt",x,delimiter='\n')