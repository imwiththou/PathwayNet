import ode

#enzyme initial condition setup
GlucoseGA3P_vmax = 
GlucoseGA3P_gamma =
GlucoseGA3P_rate =
GlucoseGA3P_km =

GA3PDHAP_vmax =
GA3PDHAP_gamma =
GA3PDHAP_rate =
GA3PDHAP_km =

DHAPG3P_vmax =
DHAPG3P_gamma =
DHAPG3P_rate =
DHAPG3P_km =

G3PGlycerol_vmax =
G3PGlycerol_gamma =
G3PGlycerol_rate =
G3PGlycerol_km =

GA3PPyruvate_vmax =
GA3PPyruvate_gamma =
GA3PPyruvate_rate = 
GA3PPyruvate_km =

PyruvateAlphaAcetolactate_vmax =
PyruvateAlphaAcetolactate_gamma =
PyruvateAlphaAcetolactate_rate =
PyruvateAlphaAcetolactate_km =

AlphaAcetolactateDiacetyl_vmax =
AlphaAcetolactateDiacetyl_gamma =
AlphaAcetolactateDiacetyl_rate =
AlphaAcetolactateDiacetyl_km =

DiacetylAcetoin_vmax =
DiacetylAcetoin_gamma =
DiacetylAcetoin_rate =
DiacetylAcetoin_km =

AcetoinBDO_vmax =
AcetoinBDO_gamma =
AcetoinBDO_rate =
AcetoinBDO_km =

PyruvateAcetylaldehyde_vmax =
PyruvateAcetylaldehyde_gamma =
PyruvateAcetylaldehyde_rate =
PyruvateAcetylaldehyde_km =

AcetylaldehydeAcetoin_vmax =
AcetylaldehydeAcetoin_gamma =
AcetylaldehydeAcetoin_rate =
AcetylaldehydeAcetoin_km =

AcetylaldehydeAcetate_vmax =
AcetylaldehydeAcetate_gamma =
AcetylaldehydeAcetate_rate =
AcetylaldehydeAcetate_km =

AcetateEthanol_vmax =
AcetateEthanol_gamma =
AcetateEthanol_rate =
AcetateEthanol_km =

#reactant initial condition setup

Glucose_vmax =
Glucose_gamma =

GA3P_gamma =
DHAP_gamma =
G3P_gamma =
Glycerol_gamma =
Pyruvate_gamma =
AlphaAcetolactate_gamma =
Diacetyl_gamma =
Acetoin_gamma =
BDO_gamma =
Acetylaldehyde_gamma =
Acetate_gamma =
Ethanol_gamma =

#ode system setup, reactant and product concentrations input respectively
y = range[26]
y[0] = #[GlucoseGA3P]
y[1] = #[GA3PDHAP]
y[2] = #[DHAPG3P]
y[3] = #[G3PGlycerol]
y[4] = #[GA3PPyruvate]
y[5] = #[PyruvateAlphaAcetolactate]
y[6] = #[AlphaAcetolactateDiacetyl]
y[7] = #[DiacetylAcetoin]
y[8] = #[AcetoinBDO]
y[9] = #[PyruvateAcetylaldehyde]
y[10] = #[AcetylaldehydeAcetoin]
y[11] = #[AcetylaldehydeAcetate]
y[12] = #[AcetateEthanol]
#-------
y[13] = #[Glucose]
y[14] = #[GA3P]
y[15] = #[DHAP]
y[16] = #[G3P]
y[17] = #[Glycerol]
y[18] = #[Pyruvate]
y[19] = #[AlphaAcetolactate]
y[20] = #[Diacetyl]
y[21] = #[Acetoin]
y[22] = #[BDO]
y[23] = #[Acetylaldehyde]
y[24] = #[Acetate]
y[25] = #[Ethanol]


#individual ordinary differential equations' definition

def Glucose(t, y):
	production = Glucose_vmax
	degradation = Glucose_gamma * [Glucose]
	usage = ([Glucose] * [GlucoseGA3P] * GlucoseGA3P_rate) / ([Glucose] + GlucoseGA3P_km)
	
def GA3P(t, y):
	production = ([Glucose] * [GlucoseGA3P] * GlucoseGA3P_rate) / ([Glucose] + GlucoseGA3P_km)
	degradation = GA3P_gamma * [GA3P]
	usage = (([GA3P] * [GA3PDHAP] * GA3PDHAP_rate) / ([GA3P] + GA3PDHAP_km)) + (([GA3P] * [GA3PPyruvate] * GA3PPyruvate_rate) / ([GA3P] * GA3PPyruvate_km))
	
def DHAP(t, y):
	production = ([GA3P] * [GA3PDHAP] * GA3PDHAP_rate) / ([GA3P] + GA3PDHAP_km)
	degradation = [DHAP] * DHAP_gamma
	usage = ([DHAP] * [DHAPG3P] * DHAPG3P_rate) / ([DHAP] + DHAPG3P_km)
	
def G3P(t, y):
	production = ([DHAP] * [DHAPG3P] * DHAPG3P_rate) / ([DHAP] + DHAPG3P_km)
	degradation = [G3P] * G3P_gamma
	usage = ([G3P] * [G3PGlycerol] * G3PGlycerol_rate) / ([G3P] + G3PGlycerol_km)
	
def Glycerol(t, y):
	production = ([G3P] * [G3PGlycerol] * G3PGlycerol_rate) / ([G3P] + G3PGlycerol_km)
	degradation = [Glycerol] * Glycerol_gamma
	usage = 0
		
#-----

def Pyruvate(t, y):
	production = ([GA3P] * [GA3PPyruvate] * GA3PPyruvate_rate) / ([GA3P] * GA3PPyruvate_km)
	degradation = [Pyruvate] * Pyruvate_gamma
	usage = (([Pyruvate] * [PyruvateAlphaAcetolactate] * PyruvateAlphaAcetolactate_rate)/([Pyruvate] + PyruvateAlphaAcetolactate_km)) + (([Pyruvate] * [PyruvateAcetylaldehyde] * PyruvateAcetylaldehyde_rate)/([Pyruvate] + PyruvateAcetylaldehyde_km))
	
def AlphaAcetolactate(t, y):
	production = ([Pyruvate] * [PyruvateAlphaAcetolactate] * PyruvateAlphaAcetolactate_rate)/([Pyruvate] + PyruvateAlphaAcetolactate_km)
	degradation = [AlphaAcetolactate] * AlphaAcetolactate_gamma
	usage = ([AlphaAcetolactate] * [AlphaAcetolactateDiacetyl] * AlphaAcetolactateDiacetyl_rate)/([AlphaAcetolactate] + AlphaAcetolactateDiacetyl_km)
	
def Diacetyl(t, y):
	production = ([AlphaAcetolactate] * [AlphaAcetolactateDiacetyl] * AlphaAcetolactateDiacetyl_rate)/([AlphaAcetolactate] + AlphaAcetolactateDiacetyl_km)
	degradation = [Diacetyl] * Diacetyl_gamma
	usage = ([Diacetyl] * [DiacetylAcetoin] * DiacetylAcetoin_rate)/([Diacetyl] + DiacetylAcetoin_km)
	
def Acetoin(t, y):
	production = (([Diacetyl] * [DiacetylAcetoin] * DiacetylAcetoin_rate)/([Diacetyl] + DiacetylAcetoin_km)) + (([Acetylaldehyde] * [AcetylaldehydeAcetoin] * AcetylaldehydeAcetoin_rate)/([Acetylaldehyde] + AcetylaldehydeAcetoin_km))
	degradation = [Acetoin] * Acetoin_gamma
	usage = ([Acetoin] * [AcetoinBDO] * AcetoinBDO_rate)/([Acetoin] + AcetoinBDO_km)
	
def BDO(t, y):
	production = ([Acetoin] * [AcetoinBDO] * AcetoinBDO_rate)/([Acetoin] + AcetoinBDO_km)
	degradation = [BDO] * BDO_gamma
	usage = 0
	
def Acetylaldehyde(t, y):
	production = ([Pyruvate] * [PyruvateAcetylaldehyde] * PyruvateAcetylaldehyde_rate)/([Pyruvate] + PyruvateAcetylaldehyde_km)
	degradation = [Acetylaldehyde] * Acetylaldehyde_gamma
	usage = (([Acetylaldehyde] * [AcetylaldehydeAcetoin] * AcetylaldehydeAcetoin_rate)/([Acetylaldehyde] + AcetylaldehydeAcetoin_km)) + (([Acetylaldehyde] * [AcetylaldehydeAcetate] * AcetylaldehydeAcetate_rate)/([Acetylaldehyde] + AcetylaldehydeAcetate_km))
	
def Acetate(t, y):
	production = ([Acetylaldehyde] * [AcetylaldehydeAcetate] * AcetylaldehydeAcetate_rate)/([Acetylaldehyde] + AcetylaldehydeAcetate_km)
	degradation = [Acetate] * Acetate_gamma
	usage = ([Acetate] * [AcetateEthanol] * AcetateEthanol_rate)/([Acetate] + AcetateEthanol_km)
	
def Ethanol(t, y):
	production = ([Acetate] * [AcetateEthanol] * AcetateEthanol_rate)/([Acetate] + AcetateEthanol_km)
	degradation = [Ethanol] * Ethanol_gamma
	usage = 0
	
#-----enzymes------

def GlucoseGA3P(t, y):
	production = GlucoseGA3P_vmax
	degradation = GlucoseGA3P_gamma * [GlucoseGA3P]
	usage = 0
	
def GA3PDHAP(t, y):
	production = GA3PDHAP_vmax
	degradation = GA3PDHAP_gamma * [GA3PDHAP]
	usage = 0
	
def DHAPG3P(t, y):
	production = DHAPG3P_vmax
	degradation = DHAPG3P_gamma * [DHAPG3P]
	usage = 0
	
def G3PGlycerol(t, y):
	production = G3PGlycerol_vmax
	degradation = G3PGlycerol_gamma * [G3PGlycerol]
	usage = 0
	
def GA3PPyruvate(t, y):
	production = GA3PPyruvate_vmax 
	degradation = GA3PPyruvate_gamma * [GA3PPyruvate]
	usage = 0
	
def PyruvateAlphaAcetolactate(t, y):
	production = PyruvateAlphaAcetolactate_vmax
	degradation = PyruvateAlphaAcetolactate_gamma * [PyruvateAlphaAcetolactate]
	usage = 0
	
def AlphaAcetolactateDiacetyl(t, y):
	production = AlphaAcetolactateDiacetyl_vmax
	degradation = AlphaAcetolactateDiacetyl_gamma * [AlphaAcetolactateDiacetyl]
	usage = 0
	
def DiacetylAcetoin(t, y):
	production = DiacetylAcetoin_vmax
	degradation = DiacetylAcetoin_gamma * [DiacetylAcetoin]
	usage = 0
	
def AcetoinBDO(t, y):
	production = AcetoinBDO_vmax
	degradation = AcetoinBDO_gamma * [AcetoinBDO]
	usage = 0
	
def PyruvateAcetylaldehyde(t, y):
	production = PyruvateAcetylaldehyde_vmax
	degradation = PyruvateAcetylaldehyde_gamma * [PyruvateAcetylaldehyde]
	usage = 0
	
def AcetylaldehydeAcetoin(t, y):
	production = AcetylaldehydeAcetoin_vmax
	degradation =AcetylaldehydeAcetoin_gamma * [AcetylaldehydeAcetoin]
	usage =0
	
def AcetylaldehydeAcetate(t, y):
	production = AcetylaldehydeAcetate_vmax
	degradation = AcetylaldehydeAcetate_gamma * [AcetylaldehydeAcetate]
	usage = 0
	
def AcetateEthanol(t, y):
	production = AcetateEthanol_vmax
	degradation = AcetateEthanol_gamma * [AcetateEthanol]
	usage = 0
	
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
tmax = 10000.0
dt = 0.1
outfile = 'IntegratedSimulationResult.csv'
f = open(outfile, 'w')
header = ['time', 'GlucoseGA3P', 'GA3PDHAP', 'DHAPG3P', 'G3PGlycerol', 'GA3PPyruvate', 'PyruvateAlphaAcetolactate', 'AlphaAcetolactateDiacetyl', 'DiacetylAcetoin', 'AcetoinBDO', 'PyruvateAcetylaldehyde', 'AcetylaldehydeAcetoin', 'AcetylaldehydeAcetate', 'AcetateEthanol', 'Glucose', 'GA3P', 'DHAP', 'G3P', 'Glycerol', 'Pyruvate', 'AlphaAcetolactate', 'Diacetyl', 'Acetoin', 'BDO', 'Acetylaldehyde', 'Acetate', 'Ethanol']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()
