import ode

#enzyme initial condition setup
GlucoseGA3P_vmax = 0
GlucoseGA3P_gamma = 0
GlucoseGA3P_rate = 1.1
GlucoseGA3P_km = 4.4e-04 #confirmed

GA3PDHAP_vmax = 0
GA3PDHAP_gamma = 0
GA3PDHAP_rate = 4.4
GA3PDHAP_km = 3e-04 #confirmed

DHAPG3P_vmax = 0
DHAPG3P_gamma = 0
DHAPG3P_rate = 3.14
DHAPG3P_km = 1.5e-04 #confirmed

G3PGlycerol_vmax = 0
G3PGlycerol_gamma = 0
G3PGlycerol_rate = 3.14
G3PGlycerol_km = 5e-03 #confirmed

GA3PPyruvate_vmax = 0
GA3PPyruvate_gamma = 0
GA3PPyruvate_rate = 5.6
GA3PPyruvate_km = 8.86e-02 #confirmed

PyruvateAlphaAcetolactate_vmax = 7.5e-11 #confirmed, conflict with model assumption
PyruvateAlphaAcetolactate_gamma = 0
PyruvateAlphaAcetolactate_rate = 3.14
PyruvateAlphaAcetolactate_km = 6.33e-03 #confirmed

AlphaAcetolactateDiacetyl_vmax = 0
AlphaAcetolactateDiacetyl_gamma = 0
AlphaAcetolactateDiacetyl_rate = 3.14
AlphaAcetolactateDiacetyl_km = 0.269

DiacetylAcetoin_vmax = 0
DiacetylAcetoin_gamma = 0
DiacetylAcetoin_rate = 41.0 #confirmed
DiacetylAcetoin_km = 6.7e-04 #confirmed

AcetoinBDO_vmax =  0 #3.9e-04 #confirmed, conflict with model assuption
AcetoinBDO_gamma = 0
AcetoinBDO_rate = 16 #confirmed
AcetoinBDO_km = 1.175e-01 #confirmed

PyruvateAcetylaldehyde_vmax = 0
PyruvateAcetylaldehyde_gamma = 0
PyruvateAcetylaldehyde_rate = 12.9 #confirmed
PyruvateAcetylaldehyde_km = 5.61e-04 #confirmed

AcetylaldehydeAcetoin_vmax = 0
AcetylaldehydeAcetoin_gamma = 0
AcetylaldehydeAcetoin_rate = 3.14
AcetylaldehydeAcetoin_km = 0.269

AcetylaldehydeAcetate_vmax = 0
AcetylaldehydeAcetate_gamma = 0
AcetylaldehydeAcetate_rate = 20
AcetylaldehydeAcetate_km = 2e-04

AcetateEthanol_vmax = 0
AcetateEthanol_gamma = 0
AcetateEthanol_rate = 17.0 #confirmed as 170 but the value is ridiculiously high
AcetateEthanol_km = 4.5e-05 #confirmed

#reactant initial condition setup

Glucose_vmax = 0
Glucose_gamma = 0

GA3P_gamma = 0
DHAP_gamma = 0
G3P_gamma = 0
Glycerol_gamma = 0
Pyruvate_gamma = 0
AlphaAcetolactate_gamma = 0
Diacetyl_gamma = 0
Acetoin_gamma = 0
BDO_gamma = 0
Acetylaldehyde_gamma = 0
Acetate_gamma = 0
Ethanol_gamma = 0

#ode system setup, reactant and product concentrations input respectively
y = range(26)
y[0] = 1e-06 #[GlucoseGA3P]
y[1] = 1e-06#[GA3PDHAP]
y[2] = 1e-06#[DHAPG3P]
y[3] = 1e-06#[G3PGlycerol]
y[4] = 1e-06#[GA3PPyruvate]
y[5] = 1e-06#[PyruvateAlphaAcetolactate]
y[6] = 1e-06#[AlphaAcetolactateDiacetyl]
y[7] = 1e-06#[DiacetylAcetoin]
y[8] = 1e-06#[AcetoinBDO]
y[9] = 1e-06#[PyruvateAcetylaldehyde]
y[10] = 1e-06#[AcetylaldehydeAcetoin]
y[11] = 1e-06#[AcetylaldehydeAcetate]
y[12] = 1e-06#[AcetateEthanol]
#-------
y[13] = 1e-05#[Glucose]
y[14] = 0#[GA3P]
y[15] = 0#[DHAP]
y[16] = 0#[G3P]
y[17] = 0#[Glycerol]
y[18] = 0#[Pyruvate]
y[19] = 0#[AlphaAcetolactate]
y[20] = 0#[Diacetyl]
y[21] = 0#[Acetoin]
y[22] = 0#[BDO]
y[23] = 0#[Acetylaldehyde]
y[24] = 0#[Acetate]
y[25] = 0#[Ethanol]


#individual ordinary differential equations' definition

def Glucose(t, y):
	production = Glucose_vmax
	degradation = Glucose_gamma * y[13]
	usage = (y[13] * y[0] * GlucoseGA3P_rate) / (y[13] + GlucoseGA3P_km)
	return production - degradation - usage
	
def GA3P(t, y):
	production = (y[13] * y[0] * GlucoseGA3P_rate) / (y[13] + GlucoseGA3P_km)
	degradation = GA3P_gamma * y[14]
	usage1 = (y[14] * y[1] * GA3PDHAP_rate) / (y[14] + GA3PDHAP_km)
	usage2 = (y[14] * y[4] * GA3PPyruvate_rate) / (y[14] + GA3PPyruvate_km)
	usage = usage1 + usage2
	return production - degradation - usage
	
	
def DHAP(t, y):
	production = (y[14] * y[1] * GA3PDHAP_rate) / (y[14] + GA3PDHAP_km)
	degradation = y[15] * DHAP_gamma
	usage = (y[15] * y[2] * DHAPG3P_rate) / (y[15] + DHAPG3P_km)
	return production - degradation - usage
	
	
def G3P(t, y):
	production = (y[15] * y[2] * DHAPG3P_rate) / (y[15] + DHAPG3P_km)
	degradation = y[16] * G3P_gamma
	usage = (y[16] * y[3] * G3PGlycerol_rate) / (y[16] + G3PGlycerol_km)
	return production - degradation - usage

	
def Glycerol(t, y):
	production = (y[16] * y[3] * G3PGlycerol_rate) / (y[16] + G3PGlycerol_km)
	degradation = y[17] * Glycerol_gamma
	usage = 0
	return production - degradation - usage
	
		
#-----

def Pyruvate(t, y):
	production = (y[14] * y[4] * GA3PPyruvate_rate) / (y[14] + GA3PPyruvate_km)
	degradation = y[18] * Pyruvate_gamma
	usage1 = (y[18] * y[5] * PyruvateAlphaAcetolactate_rate)/(y[18] + PyruvateAlphaAcetolactate_km)
	usage2 = (y[18] * y[9] * PyruvateAcetylaldehyde_rate)/(y[18] + PyruvateAcetylaldehyde_km)
	usage = usage1 + usage2
	return production - degradation - usage
	
	
def AlphaAcetolactate(t, y):
	production = (y[18] * y[5] * PyruvateAlphaAcetolactate_rate)/(y[18] + PyruvateAlphaAcetolactate_km)
	degradation = y[19] * AlphaAcetolactate_gamma
	usage = (y[19] * y[6] * AlphaAcetolactateDiacetyl_rate)/(y[19] + AlphaAcetolactateDiacetyl_km)
	return production - degradation - usage
	
	
def Diacetyl(t, y):
	production = (y[19] * y[6] * AlphaAcetolactateDiacetyl_rate)/(y[19] + AlphaAcetolactateDiacetyl_km)
	degradation = y[20] * Diacetyl_gamma
	usage = (y[20] * y[7] * DiacetylAcetoin_rate)/(y[20] + DiacetylAcetoin_km)
	return production - degradation - usage
	
	
def Acetoin(t, y):
	production1 = (y[20] * y[7] * DiacetylAcetoin_rate)/(y[20] + DiacetylAcetoin_km)
	production2 = (y[23] * y[10] * AcetylaldehydeAcetoin_rate)/(y[23] + AcetylaldehydeAcetoin_km)
	production = production1 + production2
	degradation = y[21] * Acetoin_gamma
	usage = (y[21] * y[8] * AcetoinBDO_rate)/(y[21] + AcetoinBDO_km)
	return production - degradation - usage
	
	
def BDO(t, y):
	production = (y[21] * y[8] * AcetoinBDO_rate)/(y[21] + AcetoinBDO_km)
	degradation = y[22] * BDO_gamma
	usage = 0
	return production - degradation - usage
	
	
def Acetylaldehyde(t, y):
	production = (y[18] * y[9] * PyruvateAcetylaldehyde_rate)/(y[18] + PyruvateAcetylaldehyde_km)
	degradation = y[23] * Acetylaldehyde_gamma
	usage1 = (y[23] * y[10] * AcetylaldehydeAcetoin_rate)/(y[23] + AcetylaldehydeAcetoin_km)
	usage2 = (y[23] * y[11] * AcetylaldehydeAcetate_rate)/(y[23] + AcetylaldehydeAcetate_km)
	usage = usage1 + usage2
	return production - degradation - usage
	
	
def Acetate(t, y):
	production = (y[23] * y[11] * AcetylaldehydeAcetate_rate)/(y[23] + AcetylaldehydeAcetate_km)
	degradation = y[24] * Acetate_gamma
	usage = (y[24] * y[12] * AcetateEthanol_rate)/(y[24] + AcetateEthanol_km)
	return production - degradation - usage
	
	
def Ethanol(t, y):
	production = (y[24] * y[12] * AcetateEthanol_rate)/(y[24] + AcetateEthanol_km)
	degradation = y[25] * Ethanol_gamma
	usage = 0
	return production - degradation - usage
	
	
#-----enzymes------

def GlucoseGA3P(t, y):
	production = GlucoseGA3P_vmax
	degradation = GlucoseGA3P_gamma * y[0]
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
	
	
def G3PGlycerol(t, y):
	production = G3PGlycerol_vmax
	degradation = G3PGlycerol_gamma * y[3]
	usage = 0
	return production - degradation - usage
	
	
def GA3PPyruvate(t, y):
	production = GA3PPyruvate_vmax 
	degradation = GA3PPyruvate_gamma * y[4]
	usage = 0
	return production - degradation - usage
	
	
def PyruvateAlphaAcetolactate(t, y):
	production = PyruvateAlphaAcetolactate_vmax
	degradation = PyruvateAlphaAcetolactate_gamma * y[5]
	usage = 0
	return production - degradation - usage
	
	
def AlphaAcetolactateDiacetyl(t, y):
	production = AlphaAcetolactateDiacetyl_vmax
	degradation = AlphaAcetolactateDiacetyl_gamma * y[6]
	usage = 0
	return production - degradation - usage
	
	
def DiacetylAcetoin(t, y):
	production = DiacetylAcetoin_vmax
	degradation = DiacetylAcetoin_gamma * y[7]
	usage = 0
	return production - degradation - usage
	
	
def AcetoinBDO(t, y):
	production = AcetoinBDO_vmax
	degradation = AcetoinBDO_gamma * y[8]
	usage = 0
	return production - degradation - usage
	
	
def PyruvateAcetylaldehyde(t, y):
	production = PyruvateAcetylaldehyde_vmax
	degradation = PyruvateAcetylaldehyde_gamma * y[9]
	usage = 0
	return production - degradation - usage
	
	
def AcetylaldehydeAcetoin(t, y):
	production = AcetylaldehydeAcetoin_vmax
	degradation =AcetylaldehydeAcetoin_gamma * y[10]
	usage =0
	return production - degradation - usage
	
	
def AcetylaldehydeAcetate(t, y):
	production = AcetylaldehydeAcetate_vmax
	degradation = AcetylaldehydeAcetate_gamma * y[11]
	usage = 0
	return production - degradation - usage
	
	
def AcetateEthanol(t, y):
	production = AcetateEthanol_vmax
	degradation = AcetateEthanol_gamma * y[12]
	usage = 0
	return production - degradation - usage

	#circuit ODE definitions

circuitODE = range(26)
circuitODE[0] =  GlucoseGA3P
circuitODE[1] =  GA3PDHAP
circuitODE[2] =  DHAPG3P
circuitODE[3] =  G3PGlycerol
circuitODE[4] =  GA3PPyruvate
circuitODE[5] =  PyruvateAlphaAcetolactate
circuitODE[6] =  AlphaAcetolactateDiacetyl
circuitODE[7] =  DiacetylAcetoin
circuitODE[8] =  AcetoinBDO
circuitODE[9] =  PyruvateAcetylaldehyde
circuitODE[10] =  AcetylaldehydeAcetoin
circuitODE[11] =  AcetylaldehydeAcetate
circuitODE[12] =  AcetateEthanol
circuitODE[13] =  Glucose
circuitODE[14] =  GA3P
circuitODE[15] =  DHAP
circuitODE[16] =  G3P
circuitODE[17] =  Glycerol
circuitODE[18] =  Pyruvate
circuitODE[19] =  AlphaAcetolactate
circuitODE[20] =  Diacetyl
circuitODE[21] =  Acetoin
circuitODE[22] =  BDO
circuitODE[23] =  Acetylaldehyde
circuitODE[24] =  Acetate
circuitODE[25] =  Ethanol


#iteration setup

t0 = 0.0
tmax = 5000.0
dt = 0.1
outfile = 'troubleshoot.csv'
f = open(outfile, 'w')
header = ['time', 'GlucoseGA3P', 'GA3PDHAP', 'DHAPG3P', 'G3PGlycerol', 'GA3PPyruvate', 'PyruvateAlphaAcetolactate', 'AlphaAcetolactateDiacetyl', 'DiacetylAcetoin', 'AcetoinBDO', 'PyruvateAcetylaldehyde', 'AcetylaldehydeAcetoin', 'AcetylaldehydeAcetate', 'AcetateEthanol', 'Glucose', 'GA3P', 'DHAP', 'G3P', 'Glycerol', 'Pyruvate', 'AlphaAcetolactate', 'Diacetyl', 'Acetoin', 'BDO', 'Acetylaldehyde', 'Acetate', 'Ethanol']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()