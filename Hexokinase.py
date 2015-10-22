import ode


def Hexokinase (t, y):
    production = 1e-05 * y[2] + 1e-02 * y[6]
    degradation = 1e-02 * y[1] * y[0]
    return production - degradation
    

def Glucose (t, y):
    production = 1e-05 * y[0]
    degradation = 1e-02 * y[1]
    return production - degradation
        

def Glucose_Hexokinase (t, y):
    production =  1e-02 * y[0] * y[1] + 1e-05 * y[4]
    degradation = 1e-05 * y[2] + 1e-02 * y[2] * y[3]
    return production - degradation


def ATP (t, y):
    production = 1e-05 * y[4]
    degradation = 1e-02 * y[2] * y[3]
    return production - degradation
        

def Glucose_ATP_Hexokinase (t, y):
    production = 1e-02 * y[2] * y[3]
    degradation = 1e-05 * y[4] + 1e-02 * y[4]
    return production - degradation
        

def ADP (t, y):
    association = 1e-02 * y[4]
    dissociation = 0 * y[5]
    return association - dissociation


def Hexokinase_G6P (t,y):
    production = 1e-02 * y[4]
    degradation = 1e-02 * y[6]
    return production - degradation       

def G6P (t,y):
    production = 1e-02 * y[6]
    degradation = 0 * y[7]
    return production - degradation

circuitODE = range(8)
y = range(8)        
circuitODE[0] = Hexokinase
circuitODE[1] = Glucose
circuitODE[2] = Glucose_Hexokinase
circuitODE[3] = ATP
circuitODE[4] = Glucose_ATP_Hexokinase
circuitODE[5] = ADP
circuitODE[6] = Hexokinase_G6P
circuitODE[7] = G6P

y[0] = 1e-09    # [Hexokinase]
y[1] = 1e-06    # [Glucose]
y[2] = 0.0   # [Glucose_Hexokinase]
y[3] = 1e-06    # [ATP]
y[4] = 0.0    # [Glucose_ATP_Hexokinase]
y[5] = 0.0    # [ADP]
y[6] = 0.0  #[Hexokinase_G6P]
y[7] = 0.0  #[G6P]

t0 = 0.0
tmax = 1800
dt = 0.1
outfile = 'Hexokinase_system.csv'
f = open(outfile, 'w')
header = ['time', 'Hk', 'Glu', 'Glu_Hk', 'ATP', 'Glu_ATP_Hk', 'ADP', 'Hk_G6P', 'G6P']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    print ', '.join([str(item) for item in x])
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()

''' Corresponding SciPy version
from scipy.integrate import odeint
from numpy import linspace, array
def geneticCircuit(y, t):
    return(y'[0] = 0.0 - (0*y[0]),
           y`[1] = 0.062233 - (0.0058*y[1]),
           y'[2] = (0.595*y[1]) - (0.00116*y[2]),
           (0.0+((6.817e-08*(y[5]**2))/(7.702e-05**2+y[5]**2))) - (0.0058*y[3]),
           (0.1785*y[3]) - (0.0061*y[4]),
           (16000*y[0]*y[2]) - (0*y[5]))

yinit = array([1e-09, 1e-06, 1e-06, 0.0, 0.0, 0.0])
time = linspace(0.0, 2000.0, 20000)
y = odeint(geneticCircuit, yinit, time)
'''
