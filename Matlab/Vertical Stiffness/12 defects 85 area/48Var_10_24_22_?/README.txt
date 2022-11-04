VM_structure.m
    input: 
        [# locational defects + trailing l1 l2 theta]
        # locational defects: integers
        l1 [1xnum_defects]
        l2 [1xnum_defects]
        l3 [1xnum_defects]
    output (cost):
        if the optimizer runs into any FEA problems --> 0
        otherwise:
            0 < cost < 3.12
        NOTE: NO PENALIZATION OF COST

model.py
    

test.py
    automated defects update
    NEED TO CHANGE num_defects for updated input vector
    solver configurations:
        Type = "Optimizer/CMAES"
        Is Sigma Bounded = True
        Population Size = 8 
        Viability Population Size = 3
        Max Value = 3.1e7
        Max Generations = 400
    CONSTRAINTS:
        area85
        0 < theta < 45

results:
    