start_time = 0.0
end_time = 1800.0
timestep = 0.1
resultsfile = 'riboflavin.csv'
 
reactions = '''
ribA, GTP > C01304
ribD1, C01304 > C01268
ribD2, C01268 > C04454
E1, C04454 > C04732
ribH, C04732 > C04332
ribE, C04332 > riboflavin
E2, riboflavin > FMN
'''
 
# <enzyme>, <degradation rate>, <k1>, <k2>, <k-1> where Km = <k-1>/<k1>
enzyme_kinetics = '''
ribA, 0.0, 1e-6, 1e-3, 1.0
ribD1, 0.0, 1e-6, 1e-3, 1.0
ribD2, 0.0, 1e-6, 1e-3, 1.0
E1, 0.0, 1e-6, 1e-3, 1.0
ribH, 0.0, 1e-6, 1e-3, 1.0
ribE, 0.0, 1e-6, 1e-3, 1.0
E2, 0.0, 1e-6, 1e-3, 1.0
'''
 
initial_conditions = '''
ribA = 1e-3
ribD1 = 1e-3
GTP = 1.0
ribD2 = 1e-3
E1 = 1e-3
ribH = 1e-3
ribE = 1e-3
E2 = 1e-3
'''
 
compound_gammas = '''
'''