import ode
#--------------------------------------------
#RouteD: Glu-Pv-AA-ACT-ETH

#Jan 5, 2016

#enzyme setup

#parameter not confirmed yet

GluPvd_vmax = 0.0
GluPvd_rate = 3.14
GluPvd_km = 0.112
GluPvd_gamma = 0

PvAAd_vmax = 0
PvAAd_rate = 5.2
PvAAd_km = 1.414
PvAAd_gamma = 0.2

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

Pvd_gamma = 0
AAd_gamma = 0
ACTd_gamma = 0
ETHd_gamma = 0

#ODE system setup, reactant and product concentrations input respectively

y = range(9)
#[GluPvd]
y[0] = 3e-05
#[PvAAd]
y[1] = 1e-05
#[AAACTd]
y[2] = 9e-06
#[ACTETHd]
y[3] = 1.2e-05

#[Glud]
y[4] = 1e-05
#[Pvd]
y[5] = 0.0
#[AAd]
y[6] = 0.0
#[ACTd]
y[7] = 0.0
#[ETHd]
y[8] = 0.0

#individual definition

def GluPvd(t, y):
	production = GluPvd_vmax
	degradation = GluPvd_gamma * y[0]
	usage = 0.0
	return	production - degradation - usage
	
def PvAAd(t, y):
	production = PvAAd_vmax
	degradation = PvAAd_gamma * y[1]
	usage = 0.0
	return	production - degradation - usage
	
def AAACTd(t, y):
	production = AAACTd_vmax
	degradation = AAACTd_gamma * y[2]
	usage = 0.0
	return	production - degradation - usage
	
def ACTETHd(t, y):
	production = ACTETHd_vmax
	degradation = ACTETHd_gamma * y[3]
	usage = 0.0
	return	production - degradation - usage
	

def Glud(t, y):
	production = Glud_vmax
	degradation = Glud_gamma * y[4]
	usage = (y[4] * y[0] * GluPvd_rate)/(y[4] + GluPvd_km)
	return	production - degradation - usage
	
def Pvd(t, y):
	production = (y[4] * y[0] * GluPvd_rate)/(y[4] + GluPvd_km)
	degradation = Pvd_gamma * y[5] 
	usage = (y[5] * y[1] * PvAAd_rate)/(y[5] + PvAAd_km)
	return	production - degradation - usage
	
def AAd(t, y):
	production = (y[5] * y[1] * PvAAd_rate)/(y[5] + PvAAd_km)
	degradation = AAd_gamma * y[6]
	usage = (y[6] * y[2] * AAACTd_rate)/(y[6] + AAACTd_km)
	return	production - degradation - usage
	
def ACTd(t, y):
	production = (y[6] * y[2] * AAACTd_rate)/(y[6] + AAACTd_km)
	degradation = ACTd_gamma * y[7]
	usage = (y[7] * y[3] * ACTETHd_rate)/(y[7] + ACTETHd_km)
	return	production - degradation - usage
	
def ETHd(t, y):
	production = (y[7] * y[3] * ACTETHd_rate)/(y[7] + ACTETHd_km)
	degradation = ETHd_gamma * y[8]
	usage = 0.0
	return	production - degradation - usage
	
#circuit ODE

circuitODE = range(9)
circuitODE[0] = GluPvd
circuitODE[1] = PvAAd
circuitODE[2] = AAACTd
circuitODE[3] = ACTETHd

circuitODE[4] = Glud
circuitODE[5] = Pvd
circuitODE[6] = AAd
circuitODE[7] = ACTd
circuitODE[8] = ETHd

#iteration setup

t0 = 0.0
tmax = 20000.0
dt = 0.1
outfile = 'RouteD.csv'
f = open(outfile, 'w')
header = ['time', 'GluPvd', 'PvAAd', 'AAACTd', 'ACTETHd', 'Glud', 'Pvd', 'AAd', 'ACTd', 'ETHd']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()
