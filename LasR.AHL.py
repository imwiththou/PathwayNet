import ode
 
 
def inducer_AHL(t, y):
    degradation = 0 * y[0]
    return 0.0 - degradation
         
 
def mRNA_pTetR(t, y):
    production = 0.062233 #Vmax
    degradation = 0.0058 * y[1]
    return production - degradation
         
 
def peptide_RBS1(t, y):
    production = 0.595 * y[1]
    degradation = 0.00116 * y[2]
    return production - degradation
 
 
def mRNA_pLasR(t, y):
    production = 0.0 + ((6.817e-08*(y[5]**2))/(7.702e-05**2 + y[5]**2))
    degradation = 0.0058 * y[3]
    return production - degradation
         
 
def peptide_RBS2(t, y):
    production = 0.1785 * y[3]
    degradation = 0.0061 * y[4]
    return production - degradation
         
 
def interaction_LasR_AHL(t, y):
    association = 16000 * y[0] * y[2]
    dissociation = 0 * y[5]
    return association - dissociation
  
circuitODE = range(6)
y = range(6)        
circuitODE[0] = inducer_AHL
circuitODE[1] = mRNA_pTetR
circuitODE[2] = peptide_RBS1
circuitODE[3] = mRNA_pLasR
circuitODE[4] = peptide_RBS2
circuitODE[5] = interaction_LasR_AHL
y[0] = 1e-09    # [inducer_AHL]
y[1] = 1e-06    # [mRNA_pTetR]
y[2] = 1e-06    # [peptide_RBS1]
y[3] = 0.0    # [mRNA_pLasR]
y[4] = 0.0    # [peptide_RBS2]
y[5] = 0.0    # [interaction_LasR_AHL]
 
 
t0 = 0.0
tmax = 1800.0
dt = 0.1
outfile = 'lasR_system.csv'
f = open(outfile, 'w')
header = ['time', 'AHL', 'LasR_mRNA', 'LasR', 'GFP_mRNA', 'GFP', 'LasR.AHL']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()
 
 
t0 = 0.0
tmax = 600.0
dt = 0.1
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    print ','.join([str(item) for item in x])
 
''' Corresponding SciPy version
from scipy.integrate import odeint
from numpy import linspace, array
 
def geneticCircuit(y, t):
    return(0.0 - (0*y[0]),
           0.062233 - (0.0058*y[1]),
           (0.595*y[1]) - (0.00116*y[2]),
           (0.0+((6.817e-08*(y[5]**2))/(7.702e-05**2+y[5]**2))) - (0.0058*y[3]),
           (0.1785*y[3]) - (0.0061*y[4]),
           (16000*y[0]*y[2]) - (0*y[5]))
 
yinit = array([1e-09, 1e-06, 1e-06, 0.0, 0.0, 0.0])
time = linspace(0.0, 2000.0, 20000)
y = odeint(geneticCircuit, yinit, time)
'''
