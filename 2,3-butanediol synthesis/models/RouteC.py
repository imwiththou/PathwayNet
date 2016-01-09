import ode
#-----------------------------------------------------
#RouteC:Glu-Pv-AA-AT-BDO

#enzyme setup

#parameter not confirmed yet

GluPvc_vmax = 0.0
GluPvc_rate = 7.32
GluPvc_km = 1.1
GluPvc_gamma = 0.0

PvAAc_vmax = 9.2
PvAAc_rate = 0.22
PvAAc_km = 0.717
PvAAc_gamma = 0.0

AAATc_vmax = 0
AAATc_rate = 4.2
AAATc_km = 2.11
AAATc_gamma = 0.0

ATBDOc_vmax = 0.9
ATBDOc_rate = 3.2
ATBDOc_km = 1.5
ATBDOc_gamma = 0.0

#substrate setup

Gluc_vmax = 0
Gluc_gamma = 0

Pvc_gamma = 0
AAc_gamma = 0
ATc_gamma = 0
BDOc_gamma = 0

#ODE system setup, reactant and product concentrations input respectively

y = range(9)
#[GluPvc]
y[0] = 1.1e-05
#[PvAAc]
y[1] = 1e-06
#[AAATc]
y[2] = 1e-05
#[ATBDOc]
y[3] = 7e-06

#[Gluc]
y[4] = 1e-05
#[Pvc]
y[5] = 0
#[AAc]
y[6] = 0
#[ATc]
y[7] = 0
#[BDOc]
y[8] = 0

#individual definition

def GluPvc(t, y):
	production = GluPvc_vmax
	degradation = GluPvc_gamma * y[0]
	usage = 0.0 
	return production - degradation - usage

def PvAAc(t, y):
	production = PvAAc_vmax
	degradation = PvAAc_gamma * y[1]
	usage = 0.0
	return production - degradation - usage

def AAATc(t, y):
	production = AAATc_vmax
	degradation = AAATc_gamma * y[2]
	usage = 0.0
	return production - degradation - usage
	
def ATBDOc(t, y):
	production = ATBDOc_vmax
	degradation = ATBDOc_gamma * y[3]
	usage = 0.0
	return production - degradation - usage
	
	
def Gluc(t, y):
	production = Gluc_vmax
	degradation = Gluc_gamma * y[4]
	usage = (y[4] * y[0] * GluPvc_rate)/(y[4] + GluPvc_km)
	return production - degradation - usage
	
def Pvc(t, y):
	production = (y[4] * y[0] * GluPvc_rate)/(y[4] + GluPvc_km)
	degradation = Pvc_gamma * y[5]
	usage = (y[5] * y[1] * PvAAc_rate)/(y[5] + PvAAc_km)
	return production - degradation - usage
	
def AAc(t, y):
	production = (y[5] * y[1] * PvAAc_rate)/(y[5] + PvAAc_km)
	degradation = AAc_gamma * y[6]
	usage = (y[6] * y[2] * AAATc_rate)/(y[6] + AAATc_km)
	return production - degradation - usage
	
def ATc(t, y):
	production = (y[6] * y[2] * AAATc_rate)/(y[6] + AAATc_km)
	degradation = ATc_gamma * y[7]
	usage = (y[7] * y[3] * ATBDOc_rate)/(y[7] + ATBDOc_km)
	return production - degradation - usage
	
def BDOc(t, y):
	production = (y[7] * y[3] * ATBDOc_rate)/(y[7] + ATBDOc_km)
	degradation = BDOc_gamma * y[8]
	usage = 0
	return production - degradation - usage
	

#circuit ODE

circuitODE = range(9)
circuitODE[0] = GluPvc
circuitODE[1] = PvAAc
circuitODE[2] = AAATc
circuitODE[3] = ATBDOc

circuitODE[4] = Gluc
circuitODE[5] = Pvc
circuitODE[6] = AAc
circuitODE[7] = ATc
circuitODE[8] = BDOc


#iteration setup

t0 = 0.0
tmax = 2000.0
dt = 0.1
outfile = 'RouteC.csv'
f = open(outfile, 'w')
header = ['time', 'GluPvc', 'PvAAc', 'AAATc', 'ATBDOc', 'Gluc', 'Pvc', 'AAc', 'ATc', 'BDOc']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()