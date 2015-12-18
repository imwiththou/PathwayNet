import ode
#The whole pathways networks related with 2,3-butanediol synthesis can be divided into four dependent paths
#-----------------------------------------------------
#RouteA: Glu-GA3P-DHAP-G3P-Gly

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

G3PGly_vmax = 0
G3PGly_rate = 47
G3PGly_km = 0.047
G3PGly_gamma = 0 #data for Homo sapiens

#substrate setup

Glu_vmax = 0.0
Glu_gamma = 0.0

GA3P_gamma = 0.0
DHAP_gamma = 0.0
G3P_gamma = 0.0
Gly_gamma = 0.0

#ode system setup, reactant and product concentrations input respectively

y = range(9)
#[GluGA3P]
y[0] = 4.06e-05
#[GA3PDHAP]
y[1] = 1.5e-06
#[DHAPG3P]
y[2] = 0.05
#[G3PGly]
y[3] = 1e-06
#[Glu]
y[4] = 1e-06
#[GA3P]
y[5] = 0.0
#[DHAP]
y[6] = 0.0
#[G3P]
y[7] = 0.0
#[Gly]
y[8] = 0.0


#individual definition

def GluGA3P(t, y):
	production = GluGA3P_vmax
	degradation = GluGA3P_gamma * y[0]
	usage = 0
	return production - degradation - usage

def GA3PDHAP(t, y):
	production = GA3PDHAP_vmax
	degradation = GA3PDHAP_gamma * y[1]
	usage = 0
	return production - degradation - usage
	
def DHAPG3P(t, y):
	production = DHAPG3P_vmax
	degradation = DHAPG3P_gamma * y[2]
	usage = 0
	return production - degradation - usage
	
def G3PGly(t, y):
	production = G3PGly_vmax
	degradation = G3PGly_gamma * y[3]
	usage = 0.0
	return production - degradation -usage
	
def Glu(t, y):
	production = Glu_vmax
	degradation = Glu_gamma * y[4]
	usage = (y[4] * y[5] * GluGA3P_rate) / (y[4] + GluGA3P_km)
	return production - degradation - usage
	
def GA3P(t, y):
	production = (y[4] * y[5] * GluGA3P_rate) / (y[4] + GluGA3P_km)
	degradation = GA3P_gamma * y[5]
	usage = (y[5] * y[6] * GA3PDHAP_rate) / (y[5] + GA3PDHAP_km)
	return production - degradation - usage
	
def DHAP(t, y):
	production = (y[5] * y[6] * GA3PDHAP_rate) / (y[5] +GA3PDHAP_km)
	degradation = DHAP_gamma * y[6]
	usage = (y[6] * y[7] * DHAPG3P_rate) / (y[6] + DHAPG3P_km)
	return production - degradation - usage	

def G3P(t, y):
	production = (y[6] * y[7] * DHAPG3P_rate) / (y[6] + DHAPG3P_km)
	degradation = G3P_gamma * y[7]
	usage = (y[7] * y[8] * G3PGly_rate) / (y[7] + G3PGly_km)
	return production - degradation - usage
	
def Gly(t, y):
	production = (y[7] * y[8] * G3PGly_rate) / (y[7] + G3PGly_km)
	degradation = Gly_gamma * y[8]
	usage = 0.0
	return production - degradation - usage
	

	
#circuit ODE
circuitODE = range(9)
circuitODE[0] = GluGA3P
circuitODE[1] = GA3PDHAP
circuitODE[2] = DHAPG3P
circuitODE[3] = G3PGly
circuitODE[4] = Glu
circuitODE[5] = GA3P
circuitODE[6] = DHAP
circuitODE[7] = G3P
circuitODE[8] = Gly

#iteration setup

t0 = 0.0
tmax = 20000.0
dt = 0.1
outfile = 'RouteA_updated.csv'
f = open(outfile, 'w')
header = ['time', 'GluGA3P', 'GA3PDHAP', 'DHAPG3P', 'G3PGly', 'Glu', 'GA3P', 'DHAP', 'G3P', 'Gly']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()

#-----------------------------------------------------
#RouteB:Glu-GA3P-Pv-AlphaAL-DA-AT-BDO


#-----------------------------------------------------
#RouteC:Glu-GA3P-Pv-AA-AT-BDO


#-----------------------------------------------------
#RouteD:Glu-GA3P-Pv-AA-ACT-ETH

