#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 17:18:59 2022

@author: ltlworkstation
"""
from APDL_class import APDL
from ansys.mapdl.core import launch_mapdl
import numpy as np
input = np.loadtxt("sim_results/sim8.txt")
exec_loc = "/home/ltlworkstation/Ansys/ansys_inc/v222/ansys/bin/ansys222"
mapdl = launch_mapdl(exec_loc)
apdl_class = APDL(mapdl)

parsed_input = apdl_class.parse_input_12(input) #parse the input: subject to change
apdl_class.setup_geometry(parsed_input) #introducing defects as well
apdl_class.mapdl.vsbv(1,'all')

apdl_class.mesh_apdl()

apdl_class.apply_loads_boundary_conditions(1e-3)
apdl_class.mapdl.slashsolu()
output = apdl_class.mapdl.solve()
mapdl = apdl_class.mapdl
mapdl.post1()
mapdl.set("last", "last")
mapdl.nsel("S", "LOC", "X", 0)
mapdl.fsum()
reaction_1_bottom_X = mapdl.get("REAC_1", "FSUM", "", "ITEM", "FX")
reaction_2 = mapdl.get("REAC_2", "FSUM", "", "ITEM", "FY")
print(reaction_1_bottom_X)

mapdl.nsel("S", "LOC", "X", 2)

mapdl.fsum()
reaction_2_top_X = mapdl.get("REAC_2", "FSUM", "", "ITEM", "FX")
print(reaction_2_top_X)
print(reaction_2)
stiffness = (abs(reaction_1_bottom_X)+abs(reaction_2_top_X))/2/abs(1e-3)
print(stiffness)
rforces, nnum, dof = mapdl.result.nodal_reaction_forces(0)
dof_ref = mapdl.result.result_dof(0)
sum=0
count=0
rforces[:3], nnum[:7], dof[:7], dof_ref
for i in range(len(dof)):
    if dof[i]==2:
        sum=sum+rforces[i-1]
        count=count+1
print(sum)
print(count)
