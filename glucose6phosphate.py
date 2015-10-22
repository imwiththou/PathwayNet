import ode

hexokinase_vmax = 0.0
hexokinase_rate = 1e-09
hexokinase_km = 1e-09
hexokinase_gamma = 0.0

glucose_vmax = 0.0
glucose_gamma = 0.0
glucose6phosphate_gamma = 0.0

y = range(3)
y[0] = 1e-09    #[hexokinase]
y[1] = 1e-06    #[glucose]
y[2] = 0.0      #[glucose6phosphate]

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

circuitODE = range(3)
circuitODE[0] = hexokinase
circuitODE[1] = glucose
circuitODE[2] = glucose6phosphate

t0 = 0.0
tmax = 1800.0
dt = 0.1
outfile = 'glucose6phosphate.csv'
f = open(outfile, 'w')
header = ['time', 'hexokinase', 'glucose', 'glucose6phosphate']
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()
