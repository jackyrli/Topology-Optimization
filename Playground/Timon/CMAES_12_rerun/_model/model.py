#!/usr/bin/env python

# Testfunction form CCI'2006: g09
# import matlab function
import time
import numpy as np
import matlab.engine
eng = matlab.engine.start_matlab()
print('imported matlab engine')
nex = 300.0
ney = 300.0
LX = 2.0
LY = 2.0
defects_num = 12
optimize_struct = eng.VM_structure(nex,ney,LX,LY,defects_num)
iterations = eng.generate_plot()


def model(k):
  global iterations
  vdata = k["Parameters"]
  print('vdata',vdata[0:25])
  res = eng.compute(optimize_struct, np.array(vdata))
  iterations = eng.cont_plot(-res,iterations)
  
  time.sleep(0.5)
  print('res from model.py : ',-res)
  k["F(x)"] = -res


# plot iterations

# Constraints

def g1(k):
  v = k["Parameters"]
  temp = find_solid_area(v)
  temp = -1.0*temp-0.851
  k["F(x)"] = temp
  print('solid area percentage is : ', temp)

def g3(k): #Needs to be greater than 0
  v = k["Parameters"]
  temp = find_num_defects_from_input(v)
  k["F(x)"] = temp - 10.9
  print('the number of defects is : ', temp)      
  #print(k["F(x)"])		
# helper
def find_solid_area(v1):
  v1 = np.array(v1)
  #print("v1: ", v1[])
  void_area = np.pi*v1[25:25+defects_num] @ v1[25+defects_num:25+2*defects_num].T;
  total_area = 2*2;
  solid_area_percentage = (total_area - void_area)/total_area;
  print("Python: ", solid_area_percentage)
  return solid_area_percentage
 
def find_num_defects_from_input(v1):
  v1 = np.array(v1)
  temp = np.round(v1[0:25])
  num_defects_true = np.ones((1,25)) @ temp.T
  return num_defects_true[0]

