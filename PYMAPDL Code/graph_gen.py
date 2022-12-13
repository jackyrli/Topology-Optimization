#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 15:51:37 2022

@author: ltlworkstation
"""

from APDL_class import *
import numpy as np
import matplotlib.pyplot as plt

data_result = np.loadtxt("results.txt")

best_iter = data_result[48,:].T
best_cur = 0
for i in range(7176):
    data_cur = best_iter[i]
    if data_cur > best_cur:
        best_cur = data_cur
    best_iter[i] = best_cur

csfont = {'fontname':'Times New Roman'}
plt.figure()
plt.plot(list(range(1,7177)), best_iter)
plt.xlabel("Iteration Number",**csfont)
plt.ylabel("Cost of Optimization",**csfont)
plt.title("Korali Optimization Cost Over Iterations",**csfont)
plt.rcParams["font.family"] = "Times New Roman"

#%%
np.savetxt("first.txt",data_result[:,0])
resulting_cost = data_result[48,:]
ind = np.where(resulting_cost == best_iter[740])
np.savetxt("second.txt", data_result[:,483])
np.savetxt("third.txt", data_result[:,7165])
#%% first
from APDL_class import APDL
from ansys.mapdl.core import launch_mapdl
import numpy as np
input = np.loadtxt("first.txt")
exec_loc = "/home/ltlworkstation/Ansys/ansys_inc/v222/ansys/bin/ansys222"
mapdl = launch_mapdl(exec_loc)
apdl_class = APDL(mapdl)
a, result = apdl_class.FEA_compute(input)
result.plot_principal_nodal_stress(0,'seqv',background='w',show_edges=True,text_color='k',add_text=True)
apdl_class.finish_exit_pymapdl()
#%% second
from APDL_class import APDL
from ansys.mapdl.core import launch_mapdl
import numpy as np
input = np.loadtxt("second.txt")
exec_loc = "/home/ltlworkstation/Ansys/ansys_inc/v222/ansys/bin/ansys222"
mapdl = launch_mapdl(exec_loc)
apdl_class = APDL(mapdl)
a, result = apdl_class.FEA_compute(input)
result.plot_principal_nodal_stress(0,'seqv',background='w',show_edges=True,text_color='k',add_text=True)
apdl_class.finish_exit_pymapdl()

#%% third
from APDL_class import APDL
from ansys.mapdl.core import launch_mapdl
import numpy as np
input = np.loadtxt("third.txt")
exec_loc = "/home/ltlworkstation/Ansys/ansys_inc/v222/ansys/bin/ansys222"
mapdl = launch_mapdl(exec_loc)
apdl_class = APDL(mapdl)
a, result = apdl_class.FEA_compute(input)
result.plot_principal_nodal_stress(0,'seqv',background='w',show_edges=True,text_color='k',add_text=True)
apdl_class.finish_exit_pymapdl()