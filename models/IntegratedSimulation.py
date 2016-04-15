import ode

#enzyme initial condition setup


#reactant initial condition setup


#ode system setup, reactant and product concentrations input respectively



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
circuitODE[0] = 
circuitODE[1] = 

