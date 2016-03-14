import ode
#--------------------------------------------
#RouteD: Glu-Pv-AA-ACT-ETH

#Jan 5, 2016

#enzyme setup

#parameter not confirmed yet

GluGA3Pd_vmax = 0.0
GluGA3Pd_rate = 7.32
GluGA3Pd_km = 1.1
GluGA3Pd_gamma = 0.0

GA3PPvd_vmax = 0.0
GA3PPvd_rate = 1.32
GA3PPvd_km = 0.9
GA3PPvd_gamma = 0.1

PvAAd_vmax = 1e-05
PvAAd_rate = 0.22
PvAAd_km = 0.717
PvAAd_gamma = 0.0

AAACTd_vmax = 0
AAACTd_rate = 12.1
AAACTd_km = 0.91
AAACTd_gamma = 1.2

ACTETHd_vmax = 0
ACTETHd_rate = 9.1
ACTETHd_km = 0.618
ACTETHd_gamma = 0.0

#substrate setup

Glud_vmax = 0
Glud_gamma = 0

GA3Pd_gamma = 0
Pvd_gamma = 0
AAd_gamma = 0
ACTd_gamma = 0
ETHd_gamma = 0

#ODE system setup, reactant and product concentrations input respectively

y = range(11)
#[GluGA3Pd]
y[0] = 3e-05
#[GA3PPvd]
y[1] = 1e-06
#[PvAAd]
y[2] = 1e-05
#[AAACTd]
y[3] = 9e-06
#[ACTETHd]
y[4] = 1.2e-05

#[Glud]
y[5] = 1e-05
#[GA3Pd]
y[6] = 1e-06
#[Pvd]
y[7] = 0.0
#[AAd]
y[8] = 0.0
#[ACTd]
y[9] = 0.0
#[ETHd]
y[10] = 0.0

#individual definition

def GluGA3Pd(t, y):
	production = GluGA3Pd_vmax
	degradation = GluGA3Pd_gamma * [GluGA3P]
	usage = 0.0
	return	production - degradation - usage
	
def GA3PPvd(t, y):
	production = GA3PPvd_vmax
	degradation = GA3PPvd_gamma * [GA3PPv]
	usage = 0.0
	return production - degradation - usage
	
def PvAAd(t, y):
	production = PvAAd_vmax
	degradation = PvAAd_gamma * [PvAAd]
	usage = 0.0
	return	production - degradation - usage
	
def AAACTd(t, y):
	production = AAACTd_vmax
	degradation = AAACTd_gamma * [AAACTd]
	usage = 0.0
	return	production - degradation - usage
	
def ACTETHd(t, y):
	production = ACTETHd_vmax
	degradation = ACTETHd_gamma * [ACTETHd]
	usage = 0.0
	return	production - degradation - usage
	

def Glud(t, y):
	production = Glud_vmax
	degradation = Glud_gamma * [Glud]
	usage = ([Glud] * [GluGA3Pd] * GluGA3Pd_rate)/([Glud] + GluGA3P_km)
	return	production - degradation - usage
	
def GA3Pd(t, y):
	production = ([Glud] * [GluGA3Pd] * GluGA3Pd_rate)/([Glud] + GluGA3P_km)
	degradation = GA3Pd_gamma * [GA3Pd]
	usage = ([GA3Pd] * [GA3PPvd] * GA3PPvd_rate)/([GA3Pd] + GA3PPvd_km)
	
def Pvd(t, y):
	production = ([GA3Pd] * [GA3PPvd] * GA3PPvd_rate)/([GA3Pd] + GA3PPvd_km)
	degradation = Pvd_gamma * [Pvd] 
	usage = ([Pvd] * [PvAAd] * PvAAd_rate)/([Pvd] + PvAAd_km)
	return	production - degradation - usage
	
def AAd(t, y):
	production = ([Pvd] * [PvAAd] * PvAAd_rate)/([Pvd] + PvAAd_km)
	degradation = AAd_gamma * [AAd]
	usage = ([AAd] * [AAACTd] * AAACTd_rate)/([AAd] + AAACTd_km)
	return	production - degradation - usage
	
def ACTd(t, y):
	production = ([AAd] * [AAACTd] * AAACTd_rate)/([AAd] + AAACTd_km)
	degradation = ACTd_gamma * [ACTd]
	usage = ([ACTd] * [ACTETHd] * ACTETHd_rate)/([ACTd] + ACTETHd_km)
	return	production - degradation - usage
	
def ETHd(t, y):
	production = ([ACTd] * [ACTETHd] * ACTETHd_rate)/([ACTd] + ACTETHd_km)
	degradation = ETHd_gamma * [ETHd]
	usage = 0.0
	return	production - degradation - usage
	
#circuit ODE

circuitODE = range(11)
circuitODE[0] = GluGA3Pd
circuitODE[1] = GA3PPvd
circuitODE[2] = PvAAd
circuitODE[3] = AAACTd
circuitODE[4] = ACTETHd

circuitODE[5] = Glud
circuitODE[6] = GA3Pd
circuitODE[7] = Pvd
circuitODE[8] = AAd
circuitODE[9] = ACTd
circuitODE[10] = ETHd

#iteration setup

t0 = 0.0
tmax = 20000.0
dt = 0.1
outfile = 'RouteD.csv'
f = open(outfile, 'w')
header = ['time', 'Glud', 'GA3Pd', 'Pvd', 'AAd', 'ACTd', 'ETHd']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()
