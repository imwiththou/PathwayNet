import ode
#The whole pathways networks related with 2,3-butanediol synthesis can be divided into four dependent paths
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

#ode system setup

y = range(7)
#GluGA3P
y[1] =
#GA3PDHAP
y[2] =
#DHAPG3P
y[3] =
#Glu
y[4] =
#GA3P
y[5] =
#DHAP
y[6] =
#G3P
y[7] =

#ode individual definition

def GluGA3P(t, y):
	production =
	degradation =
	usage =
	return production - degradation - usage

def GA3PDHAP(t, y):
	production =
	degradation =
	useage =
	return production - degradation - usage
	
def DHAPG3P(t, y):
	production =
	degradation =
	useage =
	return production - degradation - usage
	
def Glu(t, y):
	production =
	degradation =
	useage =
	return production - degradation - usage
	
def GA3P(t, y):
	production =
	degradation =
	useage =
	return production - degradation - usage
	
def DHAP(t, y):
	production =
	degradation =
	useage =
	return production - degradation - usage	

def G3P(t, y):
	production =
	degradation =
	useage =
	return production - degradation - usage
	
#circuit ODE
circuitODE = range(7)




#RouteB:Glu-GA3P-Pv-alphaAL-DA-AT-BDO


#RouteC:Glu-GA3P-Pv-AA-AT-BDO


#RouteD:Glu-GA3P-Pv-AA-ACT-ETH

