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
y[1]=3e-06 #hexokinase concentration
y[2]=0 #glucose-hexokinase complex concentration
k1=3e-04
k2=1e-04

t0 = 0.0
tmax = 18000.0
dt = 1
outfile = 'Complex_G_H_new.csv'
f = open(outfile, 'w')
header = ['time', 'Glucose', 'Hexokinase', 'Complex_G_H']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    print ','.join([str(item) for item in x])
