#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 14:51:09 2022

@author: ltlworkstation
"""

from APDL_class import APDL
from ansys.mapdl.core import launch_mapdl
import numpy as np
input = np.loadtxt("sim_results/error1.txt")
exec_loc = "/home/ltlworkstation/Ansys/ansys_inc/v222/ansys/bin/ansys222"
mapdl = launch_mapdl(exec_loc)
apdl_class = APDL(mapdl)

apdl_class.Lz = 0.06
apdl_class.depth = 0.06

parsed_input = apdl_class.parse_input_12(input) #parse the input: subject to change
apdl_class.setup_geometry(parsed_input) #introducing defects as well
apdl_class.mapdl.vsbv(1,'all')
#apdl_class.mapdl.vplot()



parsed_input = apdl_class.parse_input_12(input) #parse the input: subject to change

apdl_class.setup_geometry(parsed_input) #introducing defects as well
apdl_class.mapdl.vsbv(1,'all')

apdl_class.mesh_apdl()
mapdl = apdl_class.mapdl
mapdl.eplot()



apdl_class.apply_loads_boundary_conditions(1e-3)
apdl_class.mapdl.slashsolu()
output = apdl_class.mapdl.solve()

print("successfully generated output")
apdl_class.mapdl.finish()
result = apdl_class.mapdl.result
max_stress = apdl_class.get_max_stress(result)
print("max stress is: ", max_stress)
apdl_class.finish_one_iteration(input, max_stress)

apdl_class.finish_exit_pymapdl()