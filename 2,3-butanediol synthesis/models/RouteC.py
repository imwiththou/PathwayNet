import ode
#-----------------------------------------------------
#RouteC:Glu-Pv-AA-AT-BDO

#enzyme setup

#parameter not confirmed yet

GluPvc_vmax = 0.0
GluPvc_rate = 
GluPvc_km = 
GluPvc_gamma = 0.0

PvAAc_vmax = 
PvAAc_rate = 
PvAAc_km = 
PvAAc_gamma = 0.0

AAATc_vmax = 0
AAATc_rate = 
AAATc_km = 
AAATc_gamma = 0.0

ATBDOc_vmax = 0.9
ATBDOc_rate = 
ATBDOc_km = 
ATBDOc_gamma = 0.0

#substrate setup

Gluc_vmax = 0
Gluc_gamma = 0

Pvc_gamma = 0
AAc_gamma = 0
ATc_gamma = 0
BDO_gamma = 0

#ODE system setup, reactant and product concentrations input respectively

y = range(9)
#[GluPvc]
y[0] = 
#[PvAAc]
y[1] = 
#[AAATc]
y[2] = 
#[ATBDOc]
y[3] = 

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
	degradation = ATBDOc * y[3]
	usage = 0.0
	return production - degradation - usage
	
	
def Gluc(t, y):
	production = Gluc_vmax
	degradation = Gluc_gamma
	usage = ([Gluc] * [GluPvc] * GluPvc_rate)/([Gluc] + GluPvc_km)
	return production - degradation - usage
	
def Pvc(t, y):
	production = ([Gluc] * [GluPvc] * GluPvc_rate)/([Gluc] + GluPvc_km)
	degradation = Pvc_gamma * [Pvc]
	usage = ([Pvc] * [PvAAc] * PvAAc_rate)/([Pvc] + PvAAc_km)
	return production - degradation - usage
	
def AAc(t, y):
	production = ([Pvc] * [PvAAc] * PvAAc_rate)/([Pvc] + PvAAc_km)
	degradation = AAc_gamma * [AAc]
	usage = ([AAc] * [AAATc] * AAATc_rate)/([AAc] + AAATc_km)
	return production - degradation - usage
	
def ATc(t, y):
	production = ([AAc] * [AAATc] * AAATc_rate)/([AAc] + AAATc_km)
	degradation = ATc_gamma * [ATc]
	usage = ([ATc] * [ATBDOc] * ATBDOc_rate)/([ATc] + ATBDOc_km)
	return production - degradation - usage
	
def BDOc(t, y):
	production = ([ATc] * [ATBDOc] * ATBDOc_rate)/([ATc] + ATBDOc_km)
	degradation = BDOc_gamma * [BDOc]
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
tmax = 20000.0
dt = 0.1
outfile = 'RouteB.csv'
f = open(outfile, 'w')
header = ['time', 'GluPvc', 'PvAAc', 'AAATc', 'ATBDOc', 'Gluc', 'Pvc', 'AAc', 'ATc', 'BDOc']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()