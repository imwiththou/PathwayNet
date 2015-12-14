import ode
#The whole pathways networks related with 2,3-butanediol synthesis can be divided into four dependent paths
#-----------------------------------------------------
#RouteA: Glu-GA3P-DHAP-G3P

#enzyme setup

GluGA3P_vmax =
GluGA3P_rate =
GluGA3P_km =
GluGA3P_gamma = 

GA3PDHAP_vmax =
GA3PDHAP_rate = 
GA3PDHAP_km = 
GA3PDHAP_gamma =

DHAPG3P_vmax =
DHAPG3P_rate =
DHAPG3P_km = 
DHAPG3P_gamma =

#substrate setup

Glu_vmax = 
Glu_gamma = 

GA3P_gamma = 
DHAP_gamma = 
G3P_gamma = 

#ode system setup, reactant and product concentrations input respectively

y = range(7)
#[GluGA3P]
y[1] = 
#[GA3PDHAP]
y[2] = 
#[DHAPG3P]
y[3] = 
#[Glu]
y[4] = 1e-06
#[GA3P]
y[5] = 0.0
#[DHAP]
y[6] = 0.0
#[G3P]
y[7] = 0.0

#ode individual definition

def GluGA3P(t, y):
	production = GluGA3P_vmax
	degradation = GluGA3P_gamma * y[4]
	usage = 0
	return production - degradation - usage

def GA3PDHAP(t, y):
	production = GA3PDHAP_vmax
	degradation = GA3PDHAP_gamma * y[5]
	usage = 0
	return production - degradation - usage
	
def DHAPG3P(t, y):
	production = DHAPG3P_vmax
	degradation = DHAPG3P_gamma * y[6]
	usage = 0
	return production - degradation - usage
	
def Glu(t, y):
	production = Glu_vmax
	degradation = Glu_gamma * y[4]
	usage = (y[4] * y[1] * GluGA3P_rate) / (y[4] + GluGA3P_km)
	return production - degradation - usage
	
def GA3P(t, y):
	production = (y[4] * y[1] * GluGA3P_rate) / (y[4] + GluGA3P_km)
	degradation = GA3P_gamma * y[5]
	usage = (y[5] * y[2] * GA3PDHAP_rate) / (y[5] +GA3PDHAP_km)
	return production - degradation - usage
	
def DHAP(t, y):
	production = (y[5] * y[2] * GA3PDHAP_rate) / (y[5] +GA3PDHAP_km)
	degradation = DHAP_gamma * y[6]
	usage = (y[6] * y[3] * DHAPG3P_rate) / (y[6]+DHAPG3P_km)
	return production - degradation - usage	

def G3P(t, y):
	production = (y[6] * y[3] * DHAPG3P_rate) / (y[6]+DHAPG3P_km)
	degradation = G3P_gamma * y[7]
	usage = 0.0
	return production - degradation - usage
	
#circuit ODE
circuitODE = range(7)
circuitODE[1] = GluGA3P
circuitODE[2] = GA3PDHAP
circuitODE[3] = DHAPG3P
circuitODE[4] = Glu
circuitODE[5] = GA3P
circuitODE[6] = DHAP
circuitODE[7] = G3P

#iteration setup

t0 = 0.0
tmax = 2000.0
dt = 0.1
outfile = 'RouteA.csv'
f = open(outfile, 'w')
header = ['time', 'GluGA3P', 'GA3PDHAP', 'DHAPG3P', 'Glu', 'GA3P', 'DHAP', 'G3P']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()

#-----------------------------------------------------
#RouteB:Glu-GA3P-Pv-alphaAL-DA-AT-BDO


#-----------------------------------------------------
#RouteC:Glu-GA3P-Pv-AA-AT-BDO


#-----------------------------------------------------
#RouteD:Glu-GA3P-Pv-AA-ACT-ETH

