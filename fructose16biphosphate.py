import ode

glucose_vmax = 0.0
glucose_rate = 1e-09
glucose_km = 1e-09
glucose_gamma = 0.0

phosphoglucoisomerase_vmax = 0.0
phosphoglucoisomerase_rate = 1e-09
phosphoglucoisomerase_km = 1e-09
phosphoglucoisomerase_gamma = 0.0

phosphoglucoisomerase2_vmax = 0.0
phosphoglucoisomerase2_rate = 1e-09
phosphoglucoisomerase2_km = 1e-09
phosphoglucoisomerase2_gamma = 0.0

phosphofructokinase_vmax = 0.0
phosphofructokinase_rate = 1e-09
phosphofructokinase_km = 1e-09
phosphofructokinase_gamma = 0.0

hexokinase_vmax = 0.0
hexokinase_gamma = 0.0
glucose6phosphate_gamma = 0.0
fructose6phosphate_gamma = 0.0
glucose6phosphatereversible_gamma = 0.0
fructose16biphosphate_gamma = 0.0


y = range(9)
y[0] = 4e-09    #[glucose]
y[1] = 1e-06    #[hexokinase]
y[2] = 0.0      #[glucose6phosphate]
y[3] = 1e-06    #[phosphoglucoisomerase]
y[4] = 0.0      #[fructose6phosphate]
y[5] = 1e-06    #[phosphoglucoisomerase2]
y[6] = 0.0      #[glucose6phosphatereversible]
y[7] = 1e-06    #[phosphofructokinase]
y[8] = 0.0      #[fructose16biphosphate]

def glucose(t, y):
    production = glucose_vmax
    degradation = glucose_gamma * y[0]
    return production - degradation

def hexokinase(t, y):
    production = hexokinase_vmax
    degradation = hexokinase_gamma * y[1]
    usage = (glucose_rate * y[0] * y[1]) / (y[1] + glucose_km)
    return production - degradation - usage

def glucose6phosphate(t, y):
    production = (glucose_rate * y[0] * y[1]) / (y[1] + glucose_km) 
    degradation = glucose6phosphate_gamma * y[2]
    usage = (phosphoglucoisomerase_rate * y[3] * y[2]) / (y[2] + phosphoglucoisomerase_km)
    return production - degradation - usage

def phosphoglucoisomerase(t, y):
    production = phosphoglucoisomerase_vmax
    degradation = phosphoglucoisomerase_gamma * y[3]
    return production - degradation

def fructose6phosphate(t, y):
    production = (phosphoglucoisomerase_rate * y[3] * y[2]) / (y[2] + phosphoglucoisomerase_km)
    degradation = fructose6phosphate_gamma * y[4]
    usage = (phosphoglucoisomerase2_rate * y[5] * y[4]) / (y[4] + phosphoglucoisomerase2_km)
    usage2 = (phosphofructokinase_rate * y[7] * y[8]) / (y[8] + phosphofructokinase_km)
    return production - degradation - usage - usage2

def phosphoglucoisomerase2(t, y):
    production = phosphoglucoisomerase2_vmax
    degradation = phosphoglucoisomerase2_gamma * y[5]
    return production - degradation

def glucose6phosphatereversible(t, y):
    production = (phosphoglucoisomerase2_rate * y[5] * y[4]) / (y[4] + phosphoglucoisomerase2_km)
    degradation = glucose6phosphatereversible_gamma * y[6]
    return production - degradation

def phosphofructokinase(t, y):
    production = phosphofructokinase_vmax
    degradation = phosphofructokinase_gamma * y[7]
    return production - degradation

def fructose16biphosphate(t, y):
    production = (phosphoglucoisomerase2_rate * y[5] * y[4]) / (y[4] + phosphoglucoisomerase2_km) + (phosphofructokinase_rate * y[7] * y[8]) / (y[8] + phosphofructokinase_km)
    degradation = fructose16biphosphate_gamma * y[8]
    return production - degradation


circuitODE = range(9)
circuitODE[0] = glucose
circuitODE[1] = hexokinase
circuitODE[2] = glucose6phosphate
circuitODE[3] = phosphoglucoisomerase
circuitODE[4] = fructose6phosphate
circuitODE[5] = phosphoglucoisomerase2
circuitODE[6] = glucose6phosphatereversible
circuitODE[7] = phosphofructokinase
circuitODE[8] = fructose16biphosphate

t0 = 0.0
tmax = 1800.0
dt = 0.1
outfile = 'fructose16biphosphate.csv'
f = open(outfile, 'w')
header = ['time', 'glucose', 'hexokinase', 'glucose6phosphate', 'phosphoglucoisomerase', 'fructose6phosphate',
          'phosphoglucoisomerase2', 'glucose6phosphatereversible', 'phosphofructokinase', 'fructose16biphosphate']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()
