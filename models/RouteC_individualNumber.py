import ode
#-----------------------------------------------------
#RouteC:Glu-GA3P-Pv-AA-AT-BDO

#enzyme setup

#parameter not confirmed yet

GluGA3Pc_vmax = 0.0
GluGA3Pc_rate = 7.32
GluGA3Pc_km = 1.1
GluGA3Pc_gamma = 0.0

GA3PPvc_vmax = 0.0
GA3PPvc_rate = 1.32
GA3PPvc_km = 0.9
GA3PPvc_gamma = 0.1

PvAAc_vmax = 0
PvAAc_rate = 0.22
PvAAc_km = 0.717
PvAAc_gamma = 0.0

AAATc_vmax = 0
AAATc_rate = 4.2
AAATc_km = 2.11
AAATc_gamma = 0.0

ATBDOc_vmax = 0
ATBDOc_rate = 3.2
ATBDOc_km = 1.5
ATBDOc_gamma = 0.0

#substrate setup

Gluc_vmax = 0
Gluc_gamma = 0

GA3Pc_gamma = 0
Pvc_gamma = 0
AAc_gamma = 0
ATc_gamma = 0
BDOc_gamma = 0

#ODE system setup, reactant and product concentrations input respectively

y = range(11)
#[GluGA3Pc]
y[0] = 1.1e-05
#[GA3PPvc]
y[1] = 1e-05
#[PvAAc]
y[2] = 1e-06
#[AAATc]
y[3] = 1e-05
#[ATBDOc]
y[4] = 7e-06

#[Gluc]
y[5] = 1e-05
#[GA3Pc]
y[6] = 0
#[Pvc]
y[7] = 0
#[AAc]
y[8] = 0
#[ATc]
y[9] = 0
#[BDOc]
y[10] = 0

#individual definition

def GluGA3Pc(t, y):
	production = GluGA3Pc_vmax
	degradation = GluGA3Pc_gamma * y[0]
	usage = 0.0 
	return production - degradation - usage

def GA3PPvc(t, y):
	production = GA3PPvc_vmax
	degradation = GA3PPvc_gamma * y[1]
	usage = 0
	return production - degradation - usage

def PvAAc(t, y):
	production = PvAAc_vmax
	degradation = PvAAc_gamma * y[2]
	usage = 0.0
	return production - degradation - usage

def AAATc(t, y):
	production = AAATc_vmax
	degradation = AAATc_gamma * y[3]
	usage = 0.0
	return production - degradation - usage
	
def ATBDOc(t, y):
	production = ATBDOc_vmax
	degradation = ATBDOc_gamma * y[4]
	usage = 0.0
	return production - degradation - usage
	
	
def Gluc(t, y):
	production = Gluc_vmax
	degradation = Gluc_gamma * y[5]
	usage = (y[5] * y[0] * GluGA3Pc_rate)/(y[5] + GluGA3Pc_km)
	return production - degradation - usage
	
def GA3Pc(t, y):
	production = (y[5] * y[0] * GluGA3Pc_rate)/(y[5] + GluGA3Pc_km)
	degradation = GA3Pc_gamma * y[6]
	usage = (y[6] * y[1] * GA3PPvc_rate)/(y[6] + GA3PPvc_km)
	return production - degradation - usage

def Pvc(t, y):
	production = (y[6] * y[1] * GA3PPvc_rate)/(y[6] + GA3PPvc_km)
	degradation = Pvc_gamma * y[7]
	usage = (y[7] * y[8] * PvAAc_rate)/(y[7] + PvAAc_km)
	return production - degradation - usage
	
def AAc(t, y):
	production = (y[7] * y[2] * PvAAc_rate)/(y[7] + PvAAc_km)
	degradation = AAc_gamma * y[8]
	usage = (y[8] * y[9] * AAATc_rate)/(y[8] + AAATc_km)
	return production - degradation - usage
	
def ATc(t, y):
	production = (y[8] * y[3] * AAATc_rate)/(y[8] + AAATc_km)
	degradation = ATc_gamma * y[9]
	usage = (y[9] * y[10] * ATBDOc_rate)/(y[9] + ATBDOc_km)
	return production - degradation - usage
	
def BDOc(t, y):
	production = (y[9] * y[4] * ATBDOc_rate)/(y[9] + ATBDOc_km)
	degradation = BDOc_gamma * y[10]
	usage = 0
	return production - degradation - usage
	

#circuit ODE

circuitODE = range(11)
circuitODE[0] = GluGA3Pc
circuitODE[1] = GA3PPvc
circuitODE[2] = PvAAc
circuitODE[3] = AAATc
circuitODE[4] = ATBDOc

circuitODE[5] = Gluc
circuitODE[6] = GA3Pc
circuitODE[7] = Pvc
circuitODE[8] = AAc
circuitODE[9] = ATc
circuitODE[10] = BDOc


#iteration setup

t0 = 0.0
tmax = 20000.0
dt = 0.1
outfile = 'RouteC.csv'
f = open(outfile, 'w')
header = ['time', 'GluGA3Pc', 'GA3PPvc', 'PvAAc', 'AAATc', 'ATBDOc', 'Gluc', 'GA3Pc', 'Pvc', 'AAc', 'ATc', 'BDOc']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()