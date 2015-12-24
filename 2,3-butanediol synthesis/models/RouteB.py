import ode
#-----------------------------------------------------
#RouteB:Glu-GA3P-Pv-AlphaAL-DA-AT-BDO

#enzyme setup

GluGA3Pb_vmax = 0.0
GluGA3Pb_rate = 4.14
GluGA3Pb_km = 0.451
GluGA3Pb_gamma = 0.0

#not confirmed yet

GA3PPvb_vmax = 0.51
GA3PPvb_rate = 4.14
GA3PPvb_km = 0.7
GA3PPvb_gamma = 0.0

PvAlphaALb_vmax = 0
PvAlphaALb_rate = 21.3
PvAlphaALb_km = 0.92
PvAlphaALb_gamma = 0.0

AlphaALDAb_vmax = 0.9
AlphaALDAb_rate = 21.5
AlphaALDAb_km = 0.672
AlphaALDAb_gamma = 0.0

DAATb_vmax = 0
DAATb_rate = 59.1
DAATb_km = 0.77
DAATb_gamma = 0.0

ATBDOb_vmax = 8.2
ATBDOb_rate = 0.11
ATBDOb_km = 0.732
ATBDOb_gamma = 0.0


#substrate setup

Glub_vmax = 0
Glub_gamma = 0

GA3Pb_gamma = 0
Pvb_gamma = 0
AlphaALb_gamma = 0
DAb_gamma = 0
ATb_gamma = 0
BDOb_gamma = 0


#ODE system setup, reactant and product concentrations input respectively

y = range(13)
#[GluGA3Pb]
y[0] = 4.07e-05
#[GA3PPvb]
y[1] = 1.5e-05
#[PvAlphaALb], i.e. acetolactate synthase
y[2] = 0.03
#[AlphaALDAb]
y[3] = 1e-06
#[DAATb]
y[4] = 1e-06
#[ATBDOb]
y[5] = 1.2e-06

#[Glub]
y[6] = 1e-05
#[GA3Pb]
y[7] = 0
#[Pvb]
y[8] = 0
#[AlphaALb]
y[9] = 0
#[DAb]
y[10] = 0
#[ATb]
y[11] = 0
#[BDOb]
y[12] = 0

#individual definition

def GluGA3Pb(t, y):
	production = GluGA3Pb_vmax
	degradation = GluGA3Pb_gamma * y[0]
	usage = 0.0 
	return production - degradation - usage

def GA3PPvb(t, y):
	production = GA3PPvb_vmax
	degradation = GA3PPvb_gamma * y[1]
	usage = 0.0
	return production - degradation - usage

def PvAlphaALb(t, y):
	production = PvAlphaALb_vmax
	degradation = PvAlphaALb_gamma * y[2]
	usage = 0.0
	return production - degradation - usage
	
def AlphaALDAb(t, y):
	production = AlphaALDAb_vmax
	degradation = AlphaALDAb_gamma * y[3]
	usage = 0.0
	return production - degradation - usage
	
def DAATb(t, y):
	production = DAATb_vmax
	degradation = DAATb_gamma * y[4]
	usage = 0.0
	return production - degradation - usage
	
def ATBDOb(t, y):
	production = ATBDOb_vmax
	degradation = ATBDOb_gamma * y[5]
	usage = 0.0
	return production - degradation - usage
	

def Glub(t, y):
	production = Glub_vmax
	degradation = Glub_gamma * y[6]
	usage = (y[6] * y[0] * GluGA3Pb_rate)/(y[6] + GluGA3Pb_km)
	return production - degradation - usage
	
def GA3Pb(t, y):
	production = (y[6] * y[0] * GluGA3Pb_rate)/(y[6] + GluGA3Pb_km)
	degradation = GA3Pb_gamma * y[7]
	usage = (y[7] * y[1] * GA3PPvb_rate)/(y[7] + GA3PPvb_km)
	return production - degradation - usage
	
def Pvb(t, y):
	production = (y[7] * y[1] * GA3PPvb_rate)/(y[7] + GA3PPvb_km)
	degradation = Pvb_gamma * y[8]
	usage = (y[8] * y[2] * PvAlphaALb_rate)/(y[8] + PvAlphaALb_km)
	return production - degradation - usage
	
def AlphaALb(t, y):
	production = (y[8] * y[2] * PvAlphaALb_rate)/(y[8] + PvAlphaALb_km)
	degradation = AlphaALb_gamma * y[9]
	usage = (y[9] * y[3] * AlphaALDAb_rate)/(y[9] + AlphaALDAb_km)
	return production - degradation - usage
	
def DAb(t, y):
	production = (y[9] * y[3] * AlphaALDAb_rate)/(y[9] + AlphaALDAb_km)
	degradation = DAb_gamma * y[10]
	usage = (y[10] * y[4] * DAATb_rate)/(y[10] + DAATb_km)
	return production - degradation - usage
	
def ATb(t, y):
	production = (y[10] * y[4] * DAATb_rate)/(y[10] + DAATb_km)
	degradation = ATb_gamma * y[11]
	usage = (y[11] * y[5] * ATBDOb_rate)/(y[11] + ATBDOb_km)
	return production - degradation - usage
	
def BDOb(t, y):
	production = (y[11] * y[5] * ATBDOb_rate)/(y[11] + ATBDOb_km)
	degradation = BDOb_gamma * y[12]
	usage = 0.0
	return production - degradation - usage

#circuit ODE

circuitODE = range(13)
circuitODE[0] = GluGA3Pb
circuitODE[1] = GA3PPvb
circuitODE[2] = PvAlphaALb
circuitODE[3] = AlphaALDAb
circuitODE[4] = DAATb
circuitODE[5] = ATBDOb

circuitODE[6] = Glub
circuitODE[7] = GA3Pb
circuitODE[8] = Pvb
circuitODE[9] = AlphaALb
circuitODE[10] = DAb
circuitODE[11] = ATb
circuitODE[12] = BDOb


#iteration setup

t0 = 0.0
tmax = 20000.0
dt = 0.1
outfile = 'RouteB.csv'
f = open(outfile, 'w')
header = ['time', 'GluGA3Pb', 'GA3PPvb', 'PvAlphaALb', 'AlphaALDAb', 'DAATb', 'ATBDOb', 'Glub', 'GA3Pb', 'Pvb', 'AlphaALb', 'DAb', 'ATb', 'BDOb']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()