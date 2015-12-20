import ode
#-----------------------------------------------------
#RouteB:Glu-GA3P-Pv-AlphaAL-DA-AT-BDO

#enzyme setup

GluGA3Pb_vmax =
GluGA3Pb_rate =
GluGA3Pb_km = 
GluGA3Pb_gamma =

GA3PPvb_vmax =
GA3PPvb_rate = 
GA3PPvb_km = 
GA3PPvb_gamma =

PvAlphaALb_vmax =
PvAlphaALb_rate =
PvAlphaALb_km =
PvAlphaALb_gamma =

AlphaAlDAb_vmax =
AlphaAlDAb_rate =
AlphaAlDAb_km =
AlphaAlDAb_gamma =

DAATb_vmax = 
DAATb_rate =
DAATb_km =
DAATb_gamma =

ATBDOb_vmax =
ATBDOb_rate =
ATBDOb_km =
ATBDOb_gamma =


#substrate setup

Glub_vmax =
Glub_gamma =

GA3Pb_gamma =
Pvb_gamma =
AlphaALb_gamma =
DAb_gamma =
ATb_gamma =
BDOb_gamma =


#ODE system setup, reactant and product concentrations input respectively

y = range(13)
#[GluGA3Pb]
y[0] = 
#[GA3PPvb]
y[1] =
#[PvAlphaALb]
y[2] =
#[AlphaALDAb]
y[3] =
#[DAATb]
y[4] =
#[ATBDOb]
y[5] =

#[Glub]
y[6] =
#[GA3Pb]
y[7] =
#[Pvb]
y[8] =
#[AlphaALb]
y[9] =
#[DAb]
y[10] =
#[ATb]
y[11] =
#[BDOb]
y[12] =

#individual definition

def GluGA3Pb(t, y):
	production = GluGA3Pb_vmax
	degradation = GluGA3Pb_gamma * [GluGA3Pb]
	usage = 0.0 
	return production - degradation - usage

def GA3PPvb(t, y):
	production = GA3PPvb_vmax
	degradation = GA3PPvb_gamma * [GA3PPvb]
	usage = 0.0
	return production - degradation - usage

def PvAlphaALb(t, y):
	production = PvAlphaALb_vmax
	degradation = PvAlphaALb_gamma * [PvAlphaALb]
	usage = 0.0
	return production - degradation - usage
	
def AlphaALDAb(t, y):
	production = AlphaALDAb_vmax
	degradation = AlphaALDAb_gamma * [AlphaALDAb]
	usage = 0.0
	return production - degradation - usage
	
def DAATb(t, y):
	production = DAATb_vmax
	degradation = DAATb_gamma * [DAATb]
	usage = 0.0
	return production - degradation - usage
	
def ATBDOb(t, y):
	production = ATBDOb_vmax
	degradation = ATBDOb_gamma * [ATBDOb]
	usage = 0.0
	return production - degradation - usage
	

def Glub(t, y):
	production = Glub_vmax
	degradation = Glub_gamma * [Glub]
	usage = ([Glub] * [GA3Pb] * GluGA3Pb_rate)/([Glub] + GluGA3Pb_km)
	return production - degradation - usage
	
def GA3Pb(t, y):
	production = ([Glub] * [GA3Pb] * GluGA3Pb_rate)/([Glub] + GluGA3Pb_km)
	degradation = GA3Pb_gamma * [GA3Pb]
	usage = ([GA3Pb] * [Pvb] * GA3PPvb_rate)/([GA3P] + GA3PPvb_km)
	return production - degradation - usage
	
def Pvb(t, y):
	production = ([GA3Pb] * [Pvb] * GA3PPvb_rate)/([GA3P] + GA3PPvb_km)
	degradation = Pvb_gamma * [Pvb]
	usage = ([Pvb] * [AlphaALb] * PvbAlphaALb_rate)/([Pvb] + PvbAlphaALb_km)
	return production - degradation - usage
	
def AlphaALb(t, y):
	production = ([Pvb] * [AlphaALb] * PvbAlphaALb_rate)/([Pvb] + PvbAlphaALb_km)
	degradation = AlphaALb_gamma * [AlphaALb]
	usage = ([AlphaALb] * [DAb] * AlphaALDAb_rate)/([AlphaALb] + AlphaALDAb_km)
	return production - degradation - usage
	
def DAb(t, y):
	production = ([AlphaALb] * [DAb] * AlphaALDAb_rate)/([AlphaALb] + AlphaALDAb_km)
	degradation = DAb_gamma * [DAb]
	usage = ([DAb] * [ATb] * DAATb_rate)/([DAb] + DAATb_km)
	return production - degradation - usage
	
def ATb(t, y):
	production = ([DAb] * [ATb] * DAATb_rate)/([DAb] + DAATb_km)
	degradation = ATb_gamma * [ATb]
	usage = ([ATb] * [BDOb] * ATBDOb_rate)/([ATb] + ATBDOb_km)
	return production - degradation - usage
	
def BDOb(t, y):
	production = ([ATb] * [BDOb] * ATBDOb_rate)/([ATb] + ATBDOb_km)
	degradation = BDOb_gamma * [BDOb]
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