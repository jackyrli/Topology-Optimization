#!/usr/bin/env python

# Testfunction form CCI'2006: g09
# import matlab function
import time
import numpy as np
from APDL_class import APDL
from ansys.mapdl.core import launch_mapdl
import numpy as np

defects_num = 12

exec_loc = "/home/ltlworkstation/Ansys/ansys_inc/v222/ansys/bin/ansys222"
mapdl = launch_mapdl(exec_loc)
apdl_class = APDL(mapdl)
print('launched pymapdl')

apdl_class.gen_korali_iteration_plot()
#prev line under construction

def model(k):
  #global iterations
  vdata = k["Parameters"]
  print("\n+++++++++++++++++++++++++++++++++++++++++++++++++++")
  print('vdata',vdata[0:12])
  res,result = apdl_class.FEA_compute(vdata)
  apdl_class.plot_korali_iterations(res-1)
  
  time.sleep(0.5)
  print('res from model.py : ',res)
  k["F(x)"] = res
  print("+++++++++++++++++++++++++++++++++++++++++++++++++++\n")

# plot iterations

# Constraints

def g1(k):
  v = k["Parameters"]
  temp = apdl_class.solid_volume_percentage(v)
  k["F(x)"] = temp - 0.85
  print('solid area percentage is from model.py: ', temp)


def g3(k):
  v = k["Parameters"]
  temp = find_num_defects_from_input(v)
  k["F(x)"] = temp - 12
  print('the number of defects is : ', temp)      
	
# helper

 
def find_num_defects_from_input(v1):
	v1 = np.array(v1)
	temp = np.round(v1[0:25])
	num_defects_true = np.ones((1,25)) @ temp.T
	return num_defects_true[0]

