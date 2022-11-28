# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 15:54:26 2022

NOTE TO SELF!!!! korali plot: python3 -m korali.plot

@author: Jacky Li
"""
import numpy as np
import matplotlib.pyplot as plt

class APDL:

    def __init__(self, mapdl):
        ##TODO##
        # finish setting up class attributes
        self.Lx = 2
        self.Ly = 2
        self.Lz = 0.05
        self.E = 1.086e8
        self.NU = 0.36
        self.depth = self.Lz
        # mapdl define
        self.mapdl = mapdl
        # iteration number and recording
        self.iter = 0
        #about optimization
        self.defects = 12
        self.errors_encountered = 0
        return

    def FEA_compute(self, input):
        ##TODO: TESTING##
        # input: angles in degrees!!!!
        # output: the maximum stress and the result variable
        
        try:
            parsed_input = self.parse_input_12(input) #parse the input: subject to change
            
            self.setup_geometry(parsed_input) #introducing defects as well
            self.mapdl.vsbv(1,'all')

            self.mesh_apdl()
            self.apply_loads_boundary_conditions(1e-3)
            self.mapdl.slashsolu()
            self.mapdl.allsel()
            output = self.mapdl.solve()
            
            print("successfully generated output")
            self.mapdl.finish()
            result = self.mapdl.result
            max_stress = self.get_max_stress(result)
            print("max stress is: ", max_stress)
            self.finish_one_iteration(input, max_stress)
            
            return max_stress, result
        
        except Exception as e:
            
            print("Exception encountered on iteration ", str(self.iter))
            print(e)
            self.errors_encountered += 1
            ##TODO##
            cost = 0
            result = None
            ##END TODO##
            print("\nresults saved in error"+str(self.errors_encountered)+".txt")
            np.savetxt('./sim_results/error'+str(self.errors_encountered)+'.txt', np.hstack((input,cost)), delimiter='/n')
            
            return cost, result
        
        
    def parse_input_12(self, input):
        ##parsed input for 48 inputs##
        # input : from Korali Opt of size 3*self.defects
        # output : parsed input
        input = np.array(input)
        input[0:12] = np.round(input[0:12])
        
        #start parsing first 12
        index_array = np.array(range(25))
        #pre line: array we are deleting from, the remaining is the 0, deleted is 1
        void_indices = np.zeros(self.defects)
        v25 = np.zeros(25)
        v_def = input[0:self.defects]
        v_def = v_def.astype(int)
        for i in range(self.defects):
            void_indices[i] = index_array[v_def[i]]
            index_array = np.delete(index_array, v_def[i])
        void_indices = void_indices.astype(int)
        v25[void_indices] = 1
        #end parsing first 12
        
        vdata = np.zeros(100)
        vdata[0:25] = v25
        count_defects = 0;
        for i in range(25):
            if vdata[i] >= 0.5: #defect in place
                vdata[2*i+25:2*i+27] = [input[self.defects+count_defects], input[2*self.defects+count_defects]]
                vdata[75+i] = input[3 * self.defects + count_defects] #angle
                count_defects += 1
            else:
                vdata[2*i+25:2*i+27] = [0,0]
                vdata[75+i] = 0;

        return vdata

    def mesh_apdl(self):
        ##TODO: comment on code##
        # input : object
        # output : finished meshing object
        self.mapdl.mshkey(0)
        self.mapdl.mshape(1, "3D")
        
        self.mapdl.et(1,'SOLID186')
        
        self.mapdl.lesize("ALL", self.depth, layer1=1)
        self.mapdl.mp('ex',1,self.E)
        self.mapdl.mp('nuxy',1,self.NU)

        self.mapdl.vmesh('all')
        #self.mapdl.eplot()
        return

    def apply_loads_boundary_conditions(self, displacement):
        ##TODO##
        # input:
        # output:
        self.mapdl.nsel('s','loc','x',0)
        self.mapdl.d('all','all',0)

        self.mapdl.nsel('s','loc','x', self.Lx)
        self.mapdl.d('all','ux',displacement)
        #mapdl.sf('all','pres',PRESSURE)

        self.mapdl.allsel()
        self.mapdl.finish()
        return



    # getting results
    def get_stiffness(self, result):
        ##TODO##
        # input:
        # output:
        return

    def get_max_stress(self, result):
        ##TODO##
        # input:
        # output:
        nnum, stress = result.principal_nodal_stress(0)
        von_mises = stress[:, -1]
        max_stress = np.nanmax(von_mises)
        print(max_stress)
        return max_stress

    def get_stiffnesspervolume(self, result):
        ##TODO##
        # input:
        # output: 
        return

    def finish_one_iteration(self, input, result_value):
        ##TODO: implement recording and plotting##
        # input:
        # output: clear up any unsaved changes in mapdl object and record iterations
        self.iter += 1
        print("sim", str(self.iter), " is finished.")
        percentage = self.solid_volume_percentage(input)
        print("Solid Area Percentage is : ", str(percentage))
        np.savetxt('./sim_results/sim'+str(self.iter)+'.txt', np.hstack((input,result_value, percentage)), delimiter='/n')
        self.mapdl.clear()
        return

    def setup_geometry(self, parsed_input):
        ##TODO: TEST##
        # input:
        # output: prepare the geometry of the base geometry
        #self.mapdl.btol(0.1)
        self.mapdl.clear()
        self.mapdl.prep7()
        self.mapdl.csys(0)
        self.mapdl.block(0, self.Lx, 0, self.Ly, 0, self.depth)
        
        #start setting up defects
        xcenter = [0.2,0.2,0.2,0.2,0.2, 0.6,0.6,0.6,0.6,0.6, 1.0,1.0,1.0,1.0,1.0, 1.4,1.4,1.4,1.4,1.4, 1.8,1.8,1.8,1.8,1.8]
        ycenter = [0.2,0.6,1.0,1.4,1.8, 0.2,0.6,1.0,1.4,1.8, 0.2,0.6,1.0,1.4,1.8, 0.2,0.6,1.0,1.4,1.8, 0.2,0.6,1.0,1.4,1.8]
        
        for i in range(25):
            if parsed_input[i] == 1:
                l1 = parsed_input[2*i+25]
                l2 = parsed_input[2*i+26]
                theta = parsed_input[i+75]
                self.draw_circle(xcenter[i], ycenter[i], l1, l2, theta)
        return self.mapdl
    
    def finish_exit_pymapdl(self):
        ##TEST DONE##
        self.mapdl.exit()
        print("exit mapdl successful")
    
    def plot_principal_nodal(self, result):
        ##TEST DONE##
        result.plot_principal_nodal_stress(0,'seqv',background='w',show_edges=True,text_color='k',add_text=True)
    
    def plot_korali_iterations(self, cost):
        plt.scatter(self.iter, cost)
        plt.pause(0.5)
        return
    
    def gen_korali_iteration_plot(self):
        fig = plt.figure()
        plt.title('Korali Optimization')
        plt.ylabel('cost of optimization')
        plt.xlabel('Iterations Number')
        
        return 1

    # utils for creating geometry
    def draw_circle(self, center_x, center_y, l1, l2, theta):
        ##TODO: I think what is better is to change this cylindrical clocal## 
        # input : necessary ellipse parameters
        # output : object with created ellipse
        self.mapdl.csys(0)
        k1 = self.mapdl.cyl4(center_x, center_y,l1,'','','',self.depth) #create circle in local coords
        self.mapdl.clocal(20, 0, center_x, center_y, 0, theta) #TODO TIMON
        k2 = self.mapdl.vlscale(k1, '', '', 1, l2/l1, '','','',1) #scale circle to ellipse
        self.mapdl.csys(0) #return to original coords sys
        #self.mapdl.vsbv(1,'all') #subtract volume
        return
    
    def solid_volume_percentage(self, v1):
        v1 = np.array(v1) 
        void_area = np.pi*v1[self.defects:2*self.defects] @ v1[2*self.defects:3*self.defects].T;
        total_area = self.Lx * self.Ly;
        solid_area_percentage = (total_area - void_area)/total_area;
        return solid_area_percentage
    
    
def plot_fun(input):
    from APDL_class import APDL
    from ansys.mapdl.core import launch_mapdl
    import numpy as np
    
    try:
        input = np.loadtxt("sim_results/"+input)
        exec_loc = "/home/ltlworkstation/Ansys/ansys_inc/v222/ansys/bin/ansys222"
        mapdl = launch_mapdl(exec_loc)
        apdl_class = APDL(mapdl)
        parsed_input = apdl_class.parse_input_12(input) #parse the input: subject to change
        apdl_class.setup_geometry(parsed_input) #introducing defects as well
        apdl_class.mapdl.vsbv(1,'all')
        apdl_class.mesh_apdl()
        apdl_class.mapdl.eplot()
        apdl_class.finish_exit_pymapdl()
    except Exception:
        print("exception")
        apdl_class.finish_exit_pymapdl()
