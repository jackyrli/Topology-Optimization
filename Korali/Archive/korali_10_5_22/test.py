#!/usr/bin/env python3
import sys
sys.path.append('_model')
from model import *
import math

# In this example, we demonstrate how Korali finds values for the
# variables that maximize the objective function, given by a
# user-provided computational model, subject to a set of
# constraints.

# Importing the computational model and constraints
import sys
sys.path.append('_model') # the model is defined in model.py in _model folder
from model import *
# Creating new experiment
import korali
e = korali.Experiment()
print('Successfully set up the experiment')


# Selecting problem type
e["Problem"]["Type"] = "Optimization"
e["Problem"]["Objective Function"] = model
e["Problem"]["Constraints"] = [g3,g1]
print('setting up problem')

# Creating 61 variables and setting their CCMA-ES bounds
for i in range(61):
  e["Variables"][i]["Name"] = "X" + str(i)
  if i in list(range(25)): # first 25 bound
    e["Variables"][i]["Lower Bound"] = 0
    e["Variables"][i]["Upper Bound"] = 1
    e["Variables"][i]["Values"] = [0,1]
    #e["Variables"][i]["Initial Standard Deviation"] = 0.5
    #added 10/17/22
    #e["Variables"][i]["Values"] = [0,1]
  elif i in list(range(25,49)): # l1 and l2 length bound
    e["Variables"][i]["Lower Bound"] = 0.02
    e["Variables"][i]["Upper Bound"] = 0.18
    e["Variables"][i]["Initial Value"] = 0.14
  else: # theta bound
    e["Variables"][i]["Lower Bound"] = -1*math.pi/2
    e["Variables"][i]["Upper Bound"] = math.pi/2
    e["Variables"][i]["Initial Value"] = 0
print('finished setting up variables')

# Configuring the constrained optimizer CCMA-ES
#e["Solver"]["Type"] = "Optimizer/CMAES"
#e["Solver"]["Is Sigma Bounded"] = True
#e["Solver"]["Population Size"] = 8
#e["Solver"]["Viability Population Size"] = 4
#e["Solver"]["Termination Criteria"]["Max Value"] = 3.1e7
#e["Solver"]["Termination Criteria"]["Max Generations"] = 3000


# Rprop
e["Solver"]["Type"] = "Optimizer/DEA"
e["Solver"]["Population Size"] = 100
e["Solver"]["Termination Criteria"]["Min Value Difference Threshold"] = 1e-12
e["Solver"]["Termination Criteria"]["Max Generations"] = 20000

# added by myself
e["Solver"]["Termination Criteria"]["Max Infeasible Resamplings"] = 1e10

# Starting Korali's Engine and running experiment
k = korali.Engine()
k.run(e)
