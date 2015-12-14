import ode
#The whole pathways networks related with 2,3-butanediol synthesis can be divided into four dependent paths
#-----------------------------------------------------
#RouteA: Glu-GA3P-DHAP-G3P

#enzyme setup

GluGA3P_vmax = 0.0
GluGA3P_rate = 4.14
GluGA3P_km = 0.451
GluGA3P_gamma = 0.0

GA3PDHAP_vmax = 0.0
GA3PDHAP_rate = 564.3791
GA3PDHAP_km = 6.4539905
GA3PDHAP_gamma = 0.0

DHAPG3P_vmax = 0.0
DHAPG3P_rate = 0.48
DHAPG3P_km = 0.46
DHAPG3P_gamma = 0.0

#substrate setup

Glu_vmax = 0.0
Glu_gamma = 0.0

GA3P_gamma = 0.0
DHAP_gamma = 0.0
G3P_gamma = 0.0

#ode system setup, reactant and product concentrations input respectively

y = range(7)
#[GluGA3P]
y[0] = 4.06e-05
#[GA3PDHAP]
y[1] = 1.5e-06
#[DHAPG3P]
y[2] = 0.05
#[Glu]
y[3] = 1e-06
#[GA3P]
y[4] = 0.0
#[DHAP]
y[5] = 0.0
#[G3P]
y[6] = 0.0

#ode individual definition

def GluGA3P(t, y):
	production = GluGA3P_vmax
	degradation = GluGA3P_gamma * y[3]
	usage = 0
	return production - degradation - usage

def GA3PDHAP(t, y):
	production = GA3PDHAP_vmax
	degradation = GA3PDHAP_gamma * y[4]
	usage = 0
	return production - degradation - usage
	
def DHAPG3P(t, y):
	production = DHAPG3P_vmax
	degradation = DHAPG3P_gamma * y[5]
	usage = 0
	return production - degradation - usage
	
def Glu(t, y):
	production = Glu_vmax
	degradation = Glu_gamma * y[3]
	usage = (y[3] * y[0] * GluGA3P_rate) / (y[3] + GluGA3P_km)
	return production - degradation - usage
	
def GA3P(t, y):
	production = (y[3] * y[0] * GluGA3P_rate) / (y[3] + GluGA3P_km)
	degradation = GA3P_gamma * y[4]
	usage = (y[4] * y[1] * GA3PDHAP_rate) / (y[4] +GA3PDHAP_km)
	return production - degradation - usage
	
def DHAP(t, y):
	production = (y[4] * y[1] * GA3PDHAP_rate) / (y[4] +GA3PDHAP_km)
	degradation = DHAP_gamma * y[5]
	usage = (y[5] * y[2] * DHAPG3P_rate) / (y[5]+DHAPG3P_km)
	return production - degradation - usage	

def G3P(t, y):
	production = (y[5] * y[2] * DHAPG3P_rate) / (y[5]+DHAPG3P_km)
	degradation = G3P_gamma * y[6]
	usage = 0.0
	return production - degradation - usage
	
#circuit ODE
circuitODE = range(7)
circuitODE[0] = GluGA3P
circuitODE[1] = GA3PDHAP
circuitODE[2] = DHAPG3P
circuitODE[3] = Glu
circuitODE[4] = GA3P
circuitODE[5] = DHAP
circuitODE[6] = G3P

#iteration setup

t0 = 0.0
tmax = 20000.0
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

