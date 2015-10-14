import ode
 
ribA_vmax = 0.0
ribA_rate = 1e-09
ribA_km = 1e-09
ribA_gamma = 0.0
 
ribD1_vmax = 0.0
ribD1_rate = 1e-09
ribD1_km = 1e-09
ribD1_gamma = 0.0
 
ribD2_vmax = 0.0
ribD2_rate = 1e-09
ribD2_km = 1e-09
ribD2_gamma = 0.0
 
E1_vmax = 0.0
E1_rate = 1e-09
E1_km = 1e-09
E1_gamma = 0.0
 
ribH_vmax = 0.0
ribH_rate = 1e-09
ribH_km = 1e-09
ribH_gamma = 0.0
 
ribE_vmax = 0.0
ribE_rate = 1e-09
ribE_km = 1e-09
ribE_gamma = 0.0
 
E2_vmax = 0.0
E2_rate = 1e-09
E2_km = 1e-09
E2_gamma = 0.0
 
GTP_vmax = 0.0
 
GTP_gamma = 0.0
C01304_gamma = 0.0
C01268_gamma = 0.0
C04454_gamma = 0.0
C04732_gamma = 0.0
C04332_gamma = 0.0
riboflavin_gamma = 0.0
FMN_gamma = 0.0
 
y = range(15)  
y[0] = 1e-09    # [ribA]
y[1] = 1e-09    # [ribD1]
y[2] = 1e-09    # [ribD2]
y[3] = 1e-09    # [E1]
y[4] = 1e-09    # [ribH]
y[5] = 1e-09    # [ribH]
y[6] = 1e-09    # [E2]
 
y[7] = 1e-06    # [GTP]
y[8] = 0.0    # [C01304]
y[9] = 0.0    # [C01268]
y[10] = 0.0    # [C04454]
y[11] = 0.0    # [C04732]
y[12] = 0.0    # [C04332]
y[13] = 0.0    # [riboflavin]
y[14] = 0.0    # [FMN]

def ribA(t, y):
    production = ribA_vmax
    degradation = ribA_gamma * y[0]
    return production - degradation
 
def ribD1(t, y):
    production = ribD1_vmax
    degradation = ribD1_gamma * y[1]
    return production - degradation
 
def ribD2(t, y):
    production = ribD2_vmax
    degradation = ribD2_gamma * y[2]
    return production - degradation
 
def E1(t, y):
    production = E1_vmax
    degradation = E1_gamma * y[3]
    return production - degradation
 
def ribH(t, y):
    production = ribH_vmax
    degradation = ribH_gamma * y[4]
    return production - degradation
 
def ribE(t, y):
    production = ribE_vmax
    degradation = ribE_gamma * y[5]
    return production - degradation
 
def E2(t, y):
    production = E2_vmax
    degradation = E2_gamma * y[6]
    return production - degradation
 
def GTP(t, y):
    production = GTP_vmax
    degradation = GTP_gamma * y[7]
    usage = (ribA_rate * y[0] * y[7]) / (y[7] + ribA_km)
    return production - degradation - usage
 
def C01304(t, y):
    production = (ribA_rate * y[0] * y[7]) / (y[7] + ribA_km)
    degradation = C01304_gamma * y[8]
    usage = (ribD1_rate * y[1] * y[8]) / (y[8] + ribD1_km)
    return production - degradation - usage
 
def C01268(t, y):
    production = (ribD1_rate * y[1] * y[8]) / (y[8] + ribD1_km)
    degradation = C01268_gamma * y[9]
    usage = (ribD2_rate * y[2] * y[9]) / (y[9] + ribD2_km)
    return production - degradation - usage
 
def C04454(t, y):
    production = (ribD2_rate * y[2] * y[9]) / (y[9] + ribD2_km)
    degradation = C04454_gamma * y[10]
    usage = (E1_rate * y[3] * y[10]) / (y[10] + E1_km)
    return production - degradation - usage
 
def C04732(t, y):
    production = (E1_rate * y[3] * y[10]) / (y[10] + E1_km)
    degradation = C04732_gamma * y[11]
    usage = (ribH_rate * y[4] * y[11]) / (y[11] + ribH_km)
    return production - degradation - usage
 
def C04332(t, y):
    production = (ribH_rate * y[4] * y[11]) / (y[11] + ribH_km)
    degradation = C04332_gamma * y[12]
    usage = (ribE_rate * y[5] * y[12]) / (y[12] + ribE_km)
    return production - degradation - usage
 
def riboflavin(t, y):
    production = (ribE_rate * y[5] * y[12]) / (y[12] + ribE_km)
    degradation = riboflavin_gamma * y[13]
    usage = (E2_rate * y[6] * y[13]) / (y[13] + E2_km)
    return production - degradation - usage
 
def FMN(t, y):
    production = (E2_rate * y[6] * y[13]) / (y[13] + E2_km)
    degradation = FMN_gamma * y[14]
    return production - degradation
 
circuitODE = range(15)
circuitODE[0] = ribA
circuitODE[1] = ribD1
circuitODE[2] = ribD1
circuitODE[3] = E1
circuitODE[4] = ribH
circuitODE[5] = ribE
circuitODE[6] = E2
 
circuitODE[7] = GTP
circuitODE[8] = C01304
circuitODE[9] = C01268
circuitODE[10] = C04454
circuitODE[11] = C04732
circuitODE[12] = C04332
circuitODE[13] = riboflavin
circuitODE[14] = FMN
 
t0 = 0.0
tmax = 1800.0
dt = 0.1
outfile = 'riboflavin.csv'
f = open(outfile, 'w')
header = ['time', 'ribA', 'ribD1', 'ribD2', 'E1', 'ribH', 'ribE', 'E2',
          'GTP', 'C01304', 'C01268', 'C04454', 'C04732', 'C04332',
          'riboflavin', 'FMN']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()
