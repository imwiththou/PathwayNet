import ode

hexokinase_vmax = 0.0
hexokinase_rate = 1e-09
hexokinase_km = 1e-09
hexokinase_gamma = 0.0

phosphoglucoisomerase_vmax = 0.0
phosphoglucoisomerase_rate = 1e-09
phosphoglucoisomerase_km = 1e-09
phosphoglucoisomerase_gamma = 0.0

phosphoglucoisomerase2_vmax = 0.0
phosphoglucoisomerase2_rate = 1e-09
phosphoglucoisomerase2_km = 1e-09
phosphoglucoisomerase2_gamma = 0.0

glucose_vmax = 0.0
glucose_gamma = 0.0
glucose6phosphate_gamma = 0.0
fructose6phosphate_gamma = 0.0
glucose6phosphatereversible_gamma = 0.0


y = range(7)
y[0] = 1e-09    #[hexokinase]
y[1] = 1e-06    #[glucose]
y[2] = 0.0      #[glucose6phosphate]
y[3] = 1e-06    #[phosphoglucoisomerase]
y[4] = 0.0      #[fructose6phosphate]
y[5] = 1e-06    #[phosphoglucoisomerase2]
y[6] = 0.0      #[glucose6phosphatereversible]

def hexokinase(t, y):
    production = hexokinase_vmax
    degradation = hexokinase_gamma * y[0]
    return production - degradation

def glucose(t, y):
    production = glucose_vmax
    degradation = glucose_gamma * y[1]
    usage = (hexokinase_rate * y[0] * y[1]) / (y[1] + hexokinase_km)
    return production - degradation - usage

def glucose6phosphate(t, y):
    production = (hexokinase_rate * y[0] * y[1]) / (y[1] + hexokinase_km) 
    degradation = glucose6phosphate_gamma * y[2]
    return production - degradation

def phosphoglucoisomerase(t, y):
    production = phosphoglucoisomerase_vmax
    degradation = phosphoglucoisomerase_gamma * y[3]
    return production - degradation

def fructose6phosphate(t, y):
    production = (phosphoglucoisomerase_rate * y[3] * y[2]) / (y[2] + phosphoglucoisomerase_km)
    degradation = fructose6phosphate_gamma * y[4]
    usage = (phosphoglucoisomerase2_rate * y[5] * y[4]) / (y[4] + phosphoglucoisomerase2_km)
    return production - degradation - usage

def phosphoglucoisomerase2(t, y):
    production = phosphoglucoisomerase2_vmax
    degradation = phosphoglucoisomerase2_gamma * y[5]
    return production - degradation

def glucose6phosphatereversible(t, y):
    production = (phosphoglucoisomerase2_rate * y[5] * y[4]) / (y[4] + phosphoglucoisomerase2_km)
    degradation = glucose6phosphatereversible_gamma * y[6]
    return production - degradation


circuitODE = range(7)
circuitODE[0] = hexokinase
circuitODE[1] = glucose
circuitODE[2] = glucose6phosphate
circuitODE[3] = phosphoglucoisomerase
circuitODE[4] = fructose6phosphate
circuitODE[5] = phosphoglucoisomerase2
circuitODE[6] = glucose6phosphatereversible

t0 = 0.0
tmax = 1800.0
dt = 0.1
outfile = 'fructose6phosphate.csv'
f = open(outfile, 'w')
header = ['time', 'hexokinase', 'glucose', 'glucose6phosphate', 'phosphoglucoisomerase', 'fructose6phosphate',
          'phosphoglucoisomerase2', 'glucose6phosphatereversible' ]
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()
