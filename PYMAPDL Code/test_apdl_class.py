# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 14:28:30 2022

@author: Jacky Li
"""
import numpy as np
from APDL_class import APDL
from ansys.mapdl.core import launch_mapdl
exec_loc = 'C:/Program Files/ANSYS Inc/ANSYS Student/v222/ansys/bin/winx64/ANSYS222.exe'
mapdl = launch_mapdl(exec_loc)
apdl_class = APDL(mapdl)

input = np.loadtxt("test.txt")

#output = apdl_class.parse_input_12(input)
apdl_class.FEA_compute(input)

#%%'
apdl_class.finish_exit_pymapdl()
#%%
from APDL_class import APDL
from ansys.mapdl.core import launch_mapdl
import numpy as np
exec_loc = 'C:/Program Files/ANSYS Inc/ANSYS Student/v222/ANSYS/bin/winx64/ANSYS222.exe'
mapdl = launch_mapdl(exec_loc)
apdl_class = APDL(mapdl)
input = np.loadtxt("test.txt")
max_stress, result = apdl_class.FEA_compute(input)
apdl_class.plot_principal_nodal(result)
