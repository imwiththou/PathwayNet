import ode


def Complex_Glucose_Hexokinase(t,y):
	production = k1*y[0]*y[1]
	degradation = k2*y[2]
	return production - degradation

def Glucose(t,y):
	production=0
	degradation=k1*y[2]
	return production - degradation

def G6P(t,y):
	production=k2*y[2]
	degradation=0
	return production - degradation
	
def Hexokinase(t,y):
	production=k2*y[1]
	degradation=k1*y[1]
	return production - degradation


#ODE representation

circuitODE=range(3)
circuitODE[0]=Glucose
circuitODE[1]=Hexokinase
circuitODE[2]=Complex_Glucose_Hexokinase

y=range(3)
y[0]=1e-06 #glucose concentration
y[1]=1e-06 #hexokinase concentration
y[2]=0 #glucose-hexokinase complex concentration
k1=3e-07
k2=1e-06

t0 = 0.0
tmax = 1800.0
dt = 1
outfile = 'lasR_system.csv'
f = open(outfile, 'w')
header = ['time', 'Glucose', 'Hexokinase', 'Complex_G_H']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()


t0 = 0.0
tmax = 600.0
dt = 1
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

yinit = array([1e-06, 1e-06, 0.0])
time = linspace(0.0, 2000.0, 20000)
y = odeint(geneticCircuit, yinit, time)
'''
