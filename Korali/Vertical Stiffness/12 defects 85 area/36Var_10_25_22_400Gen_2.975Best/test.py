#!/usr/bin/env python3

# In this example, we demonstrate how Korali finds values for the
# variables that maximize the objective function, given by a
# user-provided computational model, subject to a set of
# constraints.

# Importing the computational model and constraints
import sys
import math
sys.path.append('_model') # the model is defined in model.py in _model folder
from model import *
# Creating new experiment
import korali
e = korali.Experiment()
print('Successfully set up the experiment')

## TODO ##
num_defects = 12
## END TODO ##

# Selecting problem type

e["Problem"]["Type"] = "Optimization"
e["Problem"]["Objective Function"] = model
e["Problem"]["Constraints"] = [area85]
print('setting up problem')

# Creating 61 variables and setting their CCMA-ES bounds
for i in range(num_defects * 3): # first location, second l1, third l2, fourth l3
  e["Variables"][i]["Name"] = "X" + str(i)
  if i in list(range(num_defects)): # first 25 bound
    e["Variables"][i]["Lower Bound"] = 1
    e["Variables"][i]["Upper Bound"] = 25-i
    e["Variables"][i]["Granularity"] = 1
    #e["Variables"][i]["Initial Mean"]
    #e["Variables"][i]["Initial Standard Deviation"] = 0.5
    #added 10/17/22
    #e["Variables"][i]["Values"] = [0,1]
  elif i in list(range(num_defects,num_defects*3)): # l1 and l2 length bound
    e["Variables"][i]["Lower Bound"] = 0.02
    e["Variables"][i]["Upper Bound"] = 0.18
    e["Variables"][i]["Initial Value"] = 0.13
#   else: # theta bound
#     e["Variables"][i]["Lower Bound"] = 0
#     e["Variables"][i]["Upper Bound"] = 0
#     #e["Variables"][i]["Initial Mean"] = 0
print('finished setting up variables')

# Configuring the constrained optimizer CCMA-ES
e["Solver"]["Type"] = "Optimizer/CMAES"
e["Solver"]["Is Sigma Bounded"] = True
e["Solver"]["Population Size"] = 8 # each generation will have around 8 attempts after passing the viability regime
e["Solver"]["Viability Population Size"] = 3 # increase if constraints violated by too much times
e["Solver"]["Termination Criteria"]["Max Value"] = 3.1e7
e["Solver"]["Termination Criteria"]["Max Generations"] = 400


# added by myself
e["Solver"]["Termination Criteria"]["Max Infeasible Resamplings"] = 1e8

#possible for change
#e[“Solver”][“Mu Value”]
#e[“Solver”][“Mu Type”] 
#e["Solver"]["Mirrored Sampling"] = True # not applicable to problems with constraints
#e[“Solver”][“Max Covariance Matrix Corrections”]
#e[“Solver”][“Target Success Rate”]
#e[“Solver”][“Covariance Matrix Adaption Strength”]
#e[“Solver”][“Global Success Learning Rate”] 

# Starting Korali's Engine and running experiment
k = korali.Engine()
k.run(e)
