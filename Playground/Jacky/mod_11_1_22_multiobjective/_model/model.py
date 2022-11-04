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
## TODO ##
num_defects = 12
## END TODO ##
#optimize_struct = eng.VM_structure(nex,ney,LX,LY,num_defects)
#iterations = eng.generate_plot()
#print('before started model')

def model(p):
  print('started model')
  global iterations
  vdata = p["Parameters"]
  print('vdata',vdata[0:num_defects])
  #res1, res2 = eng.compute(optimize_struct, np.array(vdata))
  res1 = eng.compute(optimize_struct, np.array(vdata))
  iterations = eng.cont_plot(res,iterations)
  
  time.sleep(0.5)
  print('res1 from model.py : ',res1, 'res2 :', res2)
  p["F(x)"] = [res1, res2]
  #p["F(x)"] = [res1]
	
# Constraints

def area85(p):
  v = p["Parameters"]
  temp = find_solid_area(v)
  p["F(x)"] = temp - 0.85
  print('solid area percentage is : ', temp)

def g2(p):
  v = p["Parameters"]
  p["F(x)"] = -1 * find_solid_area(v) + 0.69


# helper
def find_solid_area(v1):
  v1 = np.array(v1)
  void_area = np.pi*v1[num_defects:num_defects+num_defects] @ v1[num_defects+num_defects:num_defects+2*num_defects].T;
  total_area = 2*2;
  solid_area_percentage = (total_area - void_area)/total_area;
  return solid_area_percentage


