import ode
#The whole pathways networks related with 2,3-butanediol synthesis can be divided into four dependent paths
#-----------------------------------------------------

#Dec 24, 2015

#RouteA: Glu-GA3P-DHAP-G3P-Gly
#enzyme setup

GluGA3Pa_vmax = 0.0
GluGA3Pa_rate = 4.14
GluGA3Pa_km = 0.451
GluGA3Pa_gamma = 0.0

GA3PDHAPa_vmax = 0.0
GA3PDHAPa_rate = 564.3791
GA3PDHAPa_km = 6.4539905
GA3PDHAPa_gamma = 0.0

DHAPG3Pa_vmax = 0.0
DHAPG3Pa_rate = 0.48
DHAPG3Pa_km = 0.46
DHAPG3Pa_gamma = 0.0

G3PGlya_vmax = 0
G3PGlya_rate = 47
G3PGlya_km = 0.047
G3PGlya_gamma = 0 #data for Homo sapiens

#substrate setup

Glua_vmax = 0.0
Glua_gamma = 0.0

GA3Pa_gamma = 0.0
DHAPa_gamma = 0.0
G3Pa_gamma = 0.0
Glya_gamma = 0.0

#RouteB:Glu-GA3P-Pv-AlphaAL-DA-AT-BDO
#enzyme setup

GluGA3Pb_vmax = 0.0
GluGA3Pb_rate = 4.14
GluGA3Pb_km = 0.451
GluGA3Pb_gamma = 0.0

#not confirmed yet

GA3PPvb_vmax = 0
GA3PPvb_rate = 4.14
GA3PPvb_km = 0.7
GA3PPvb_gamma = 0.0

PvAlphaALb_vmax = 0
PvAlphaALb_rate = 21.3
PvAlphaALb_km = 0.92
PvAlphaALb_gamma = 0.0

AlphaALDAb_vmax = 0
AlphaALDAb_rate = 21.5
AlphaALDAb_km = 0.672
AlphaALDAb_gamma = 0.0

DAATb_vmax = 0
DAATb_rate = 59.1
DAATb_km = 0.77
DAATb_gamma = 0.0

ATBDOb_vmax = 0
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

#RouteC:Glu-GA3P-Pv-AA-AT-BDO
#enzyme setup

#parameter not confirmed yet

GluGA3Pc_vmax = 0.0
GluGA3Pc_rate = 7.32
GluGA3Pc_km = 1.1
GluGA3Pc_gamma = 0.0

GA3PPvc_vmax = 0.0
GA3PPvc_rate = 1.32
GA3PPvc_km = 0.9
GA3PPvc_gamma = 0.1

PvAAc_vmax = 1e-05
PvAAc_rate = 0.22
PvAAc_km = 0.717
PvAAc_gamma = 0.0

AAATc_vmax = 0
AAATc_rate = 4.2
AAATc_km = 2.11
AAATc_gamma = 0.0

ATBDOc_vmax = 1e-05
ATBDOc_rate = 3.2
ATBDOc_km = 1.5
ATBDOc_gamma = 0.0

#substrate setup

Gluc_vmax = 0
Gluc_gamma = 0

GA3Pc_gamma = 0
Pvc_gamma = 0
AAc_gamma = 0
ATc_gamma = 0
BDOc_gamma = 0

#RouteD:Glu-GA3P-Pv-AA-ACT-ETH


#ode system setup, reactant and product concentrations input respectively

y = range(9)
#[GluGA3Pa]
y[0] = 4.06e-05
#[GA3PDHAPa]
y[1] = 1.5e-06
#[DHAPG3Pa]
y[2] = 0.05
#[G3PGlya]
y[3] = 1e-06
#[Glua]
y[4] = 1e-06
#[GA3Pa]
y[5] = 0.0
#[DHAPa]
y[6] = 0.0
#[G3Pa]
y[7] = 0.0
#[Glya]
y[8] = 0.0


#individual definition

def GluGA3Pa(t, y):
	production = GluGA3Pa_vmax
	degradation = GluGA3Pa_gamma * y[0]
	usage = 0
	return production - degradation - usage

def GA3PDHAPa(t, y):
	production = GA3PDHAPa_vmax
	degradation = GA3PDHAPa_gamma * y[1]
	usage = 0
	return production - degradation - usage
	
def DHAPG3Pa(t, y):
	production = DHAPG3Pa_vmax
	degradation = DHAPG3Pa_gamma * y[2]
	usage = 0
	return production - degradation - usage
	
def G3PGlya(t, y):
	production = G3PGlya_vmax
	degradation = G3PGlya_gamma * y[3]
	usage = 0.0
	return production - degradation -usage
	
def Glua(t, y):
	production = Glua_vmax
	degradation = Glua_gamma * y[4]
	usage = (y[4] * y[0] * GluGA3Pa_rate) / (y[4] + GluGA3Pa_km)
	return production - degradation - usage
	
def GA3Pa(t, y):
	production = (y[4] * y[0] * GluGA3Pa_rate) / (y[4] + GluGA3Pa_km)
	degradation = GA3Pa_gamma * y[5]
	usage = (y[5] * y[1] * GA3PDHAPa_rate) / (y[5] + GA3PDHAPa_km)
	return production - degradation - usage
	
def DHAPa(t, y):
	production = (y[5] * y[1] * GA3PDHAPa_rate) / (y[5] +GA3PDHAPa_km)
	degradation = DHAPa_gamma * y[6]
	usage = (y[6] * y[2] * DHAPG3Pa_rate) / (y[6] + DHAPG3Pa_km)
	return production - degradation - usage	

def G3Pa(t, y):
	production = (y[6] * y[2] * DHAPG3Pa_rate) / (y[6] + DHAPG3Pa_km)
	degradation = G3Pa_gamma * y[7]
	usage = (y[7] * y[3] * G3PGlya_rate) / (y[7] + G3PGlya_km)
	return production - degradation - usage
	
def Glya(t, y):
	production = (y[7] * y[3] * G3PGlya_rate) / (y[7] + G3PGlya_km)
	degradation = Glya_gamma * y[8]
	usage = 0.0
	return production - degradation - usage
	

	
#circuit ODE
circuitODE = range(44)
circuitODE[0] = GluGA3Pa
circuitODE[1] = GA3PDHAPa
circuitODE[2] = DHAPG3Pa
circuitODE[3] = G3PGlya
circuitODE[4] = Glua
circuitODE[5] = GA3Pa
circuitODE[6] = DHAPa
circuitODE[7] = G3Pa
circuitODE[8] = Glya

#iteration setup

t0 = 0.0
tmax = 20000.0
dt = 0.1
outfile = 'RouteA_updated.csv'
f = open(outfile, 'w')
header = ['time', 'GluGA3Pa', 'GA3PDHAPa', 'DHAPG3Pa', 'G3PGlya', 'Glua', 'GA3Pa', 'DHAPa', 'G3Pa', 'Glya']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()

#-----------------------------------------------------
#RouteB:Glu-GA3P-Pv-AlphaAL-DA-AT-BDO





#-----------------------------------------------------
#RouteC:Glu-GA3P-Pv-AA-AT-BDO


#-----------------------------------------------------
#RouteD:Glu-GA3P-Pv-AA-ACT-ETH

