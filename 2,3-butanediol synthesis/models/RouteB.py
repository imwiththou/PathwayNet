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

def GluFA3Pb(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage

def GA3PPvb(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage

def PvAlphaALb(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage
	
def AlphaALDAb(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage
	
def DAATb(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage
	
def ATBDOb(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage
	

def Glub(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage
	
def GA3Pb(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage
	
def Pvb(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage
	
def AlphaALb(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage
	
def DAb(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage
	
def ATb(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage
	
def BDOb(t, y):
	production =
	degradation =
	usage =
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

