#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 09:31:48 2022

@author: ltlworkstation
COMMENT: to be run once after each completion of optimization
"""
import matplotlib.pyplot as plt
import numpy as np

def save_data_into_one_file(iteration_number):
    entries = np.loadtxt("sim_results/sim1.txt")
    entries = entries.shape
    data_tot = np.zeros((entries[0],iteration_number))
    for i in range(iteration_number):
        data_cur = np.loadtxt("sim_results/sim"+str(i+1)+".txt")
        data_tot[:,i] = data_cur
    np.savetxt("results.txt", data_tot)
    return data_tot

def plot_iterations(iteration_number, data_tot, row_of_cost):
    plt.figure()
    x = list(range(1,iteration_number + 1))
    y = data_tot[row_of_cost-1,:]
    plt.scatter(x, y)
    plt.xlabel("Iterations")
    plt.ylabel("Cost")
    plt.title("KORALI OPTIMIZATION ITERATIONS PLOT")
    plt.show()
    plt.savefig('Korali_iterations_plot.png')
    max_value = np.max(y)
    max_index = np.argmax(y) + 1
    print("the largest cost appears at sim"+str(max_index)+".txt with a value of "+str(max_value)+".")
    return max_value, max_index

def plot_coupling(data_tot):
    plt.figure()
    plt.scatter(data_tot[50-1,:],data_tot[49-1,:])
    plt.xlabel("area percentage")
    plt.ylabel("cost value")
    plt.title("coupling data")
    
if __name__ == "__main__":
    iteration_number = 7176##CHANGE TO TOT NUMBER OF SIM FILES
    row_of_cost = 49 ## ROW OF COST AS IT APPEARS ON THE TXT FILE
    data_tot = save_data_into_one_file(iteration_number)
    max_value, max_index = plot_iterations(iteration_number, data_tot, row_of_cost)
    plot_coupling(data_tot)
    #%% geometry
    from APDL_class import plot_fun
    plot_fun("sim"+str(max_index)+".txt")
    
    #%%
    from APDL_class import plot_fun
    plot_fun("sim"+str(7167)+".txt")