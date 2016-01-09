import ode
#--------------------------------------------
#RouteD: Glu-Pv-AA-ACT-ETH

#Jan 5, 2016

#enzyme setup

#parameter not confirmed yet

GluPvd_vmax = 0.0
GluPvd_rate = 
GluPvd_km =
GluPvd_gamma =

PvAAd_vmax = 
PvAAd_rate =
PvAAd_km =
PvAAd_gamma =

AAACTd_vmax =
AAACTd_rate =
AAACTd_km =
AAACTd_gamma =

ACTETHd_vmax = 
ACTETHd_rate =
ACTETHd_km =
ACTETHd_gamma =

#substrate setup

Glud_vmax =
Glud_gamma =

Pvd_gamma =
AAd_gamma =
ACTd_gamma =
ETHd_gamma =

#ODE system setup, reactant and product concentrations input respectively

y = range(9)
#[GluPvd]
y[0] = 
#[PvAAd]
y[1] = 
#[AAACTd]
y[2] = 
#[ACTETHd]
y[3] = 

#[Glud]
y[4] =
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
	degradation = GluPvd_gamma * [GluPvd]
	usage = 0.0
	return	production - degradation - usage
	
def PvAAd(t, y):
	production = PvAAd_vmax
	degradation = PvAAd_gamma * [PvAAd]
	usage = 0.0
	return	production - degradation - usage
	
def AAACTd(t, y):
	production = AAACT_vmax
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
	usage = ([Glud] * [GluPvd] * GluPvd_rate)/([Glud] + GluPvd_km)
	return	production - degradation - usage
	
def Pvd(t, y):
	production = ([Glud] * [GluPvd] * GluPvd_rate)/([Glud] + GluPvd_km)
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
	degradation = ETHd_gamma * [Ethd]
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
outfile = 'RouteB.csv'
f = open(outfile, 'w')
header = ['time', 'GluPvd', 'PvAAd', 'AAACTd', 'ACTETHd', 'Glud', 'Pvd', 'AAd', 'ACTd', 'ETHd']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()
