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

phosphofructokinase_vmax = 0.0
phosphofructokinase_rate = 1e-09
phosphofructokinase_km = 1e-09
phosphofructokinase_gamma = 0.0

fructosephosphatealdolase_vmax = 0.0
fructosephosphatealdolase_rate = 1e-09
fructosephosphatealdolase_km = 1e-09
fructosephosphatealdolase_gamma = 0.0

fructosephosphatealdolase2_vmax = 0.0
fructosephosphatealdolase2_rate = 1e-09
fructosephosphatealdolase2_km = 1e-09
fructosephosphatealdolase2_gamma = 0.0

g6pdehydrogenase_vmax = 0.0
g6pdehydrogenase_rate = 1e-09
g6pdehydrogenase_km = 1e-09
g6pdehydrogenase_gamma = 0.0

glucose_vmax = 0.0
glucose_gamma = 0.0
g6p_gamma = 0.0
f6p_gamma = 0.0
f16p_gamma = 0.0
g3p_gamma = 0.0
dhp_gamma = 0.0
pg6_gamma = 0.0

y = range(14)
y[0] = 1e-09    #[hexokinase]
y[1] = 1e-06    #[glucose]
y[2] = 0.0      #[g6p]
y[3] = 1e-09    #[phosphoglucoisomerase]
y[4] = 1e-09    #[phosphoglucoisomerase2]
y[5] = 0.0      #[f6p]
y[6] = 1e-09    #[phosphofructokinase]
y[7] = 0.0      #[f16p]
y[8] = 1e-09    #[fructosephosphatealdolase]
y[9] = 1e-09    #[fructosephosphatealdolase2]
y[10] = 0.0     #[g3p]
y[11] = 0.0     #[dhp]
y[12] = 1e-09   #[g6pdehydrogenase]
y[13] = 0.0     #[pg6]

def hexokinase(t, y):
    production = hexokinase_vmax
    degradation = hexokinase_gamma * y[0]
    return production - degradation

def glucose(t, y):
    production = glucose_vmax
    degradation = glucose_gamma * y[1]
    usage = (hexokinase_rate * y[0] * y[1]) / (y[1] + hexokinase_km)
    return production - degradation - usage

def g6p(t, y):
    production = (hexokinase_rate * y[0] * y[1]) / (y[1] + hexokinase_km) + (phosphoglucoisomerase2_rate * y[4] * y[5]) / (y[5] + phosphoglucoisomerase2_km)
    degradation = g6p_gamma * y[2]
    usage = (phosphoglucoisomerase_rate * y[3] * y[2]) / (y[2] + phosphoglucoisomerase_km) + (g6pdehydrogenase_rate * y[12] * y[2]) / (y[2] + g6pdehydrogenase_km)
    return production - degradation - usage

def phosphoglucoisomerase(t, y):
    production = phosphoglucoisomerase_vmax
    degradation = phosphoglucoisomerase_gamma * y[3]
    return production - degradation

def phosphoglucoisomerase2(t, y):
    production = phosphoglucoisomerase2_vmax
    degradation = phosphoglucoisomerase2_gamma * y[4]
    return production - degradation

def f6p(t, y):
    production = (phosphoglucoisomerase_rate * y[3] * y[2]) / (y[2] + phosphoglucoisomerase_km)
    degradation = f6p_gamma * y[5]
    usage = (phosphoglucoisomerase2_rate * y[4] * y[5]) / (y[5] + phosphoglucoisomerase2_km) + (phosphofructokinase_rate * y[6] * y[5]) / (y[5] + phosphofructokinase_km)
    return production - degradation - usage

def phosphofructokinase(t, y):
    production = phosphofructokinase_vmax
    degradation = phosphofructokinase_gamma * y[6]
    return production - degradation

def f16p(t, y):
    production = (phosphofructokinase_rate * y[6] * y[5]) / (y[5] + phosphofructokinase_km) + (fructosephosphatealdolase2_rate * y[9] * y[10]) / (y[10] + fructosephosphatealdolase_km) + (fructosephosphatealdolase2_rate * y[9] * y[11]) / (y[11] + fructosephosphatealdolase_km)
    degradation = f16p_gamma * y[7]
    usage = (fructosephosphatealdolase_rate * y[8] * y[7]) / (y[7] + fructosephosphatealdolase_km)
    return production - degradation - usage

def fructosephosphatealdolase(t, y):
    production = fructosephosphatealdolase_vmax
    degradation = fructosephosphatealdolase_gamma * y[8]
    return production - degradation

def fructosephosphatealdolase2(t, y):
    production = fructosephosphatealdolase2_vmax
    degradation = fructosephosphatealdolase2_gamma * y[9]
    return production - degradation

def g3p(t, y):
    production = (fructosephosphatealdolase_rate * y[8] * y[7]) / (y[7] + fructosephosphatealdolase_km)
    degradation = g3p_gamma * y[10]
    usage = (fructosephosphatealdolase2_rate * y[9] * y[10]) / (y[10] + fructosephosphatealdolase_km)
    return production - degradation - usage

def dhp(t, y):
    production = (fructosephosphatealdolase_rate * y[8] * y[7]) / (y[7] + fructosephosphatealdolase_km)
    degradation = dhp_gamma * y[11]
    usage = (fructosephosphatealdolase2_rate * y[9] * y[11]) / (y[11] + fructosephosphatealdolase_km)
    return production - degradation - usage

def g6pdehydrogenase(t, y):
    production = g6pdehydrogenase_vmax
    degradation = g6pdehydrogenase_gamma * y[12]
    return production - degradation

def pg6(t, y):
    production = (g6pdehydrogenase_rate * y[12] * y[2]) / (y[2] + g6pdehydrogenase_km)
    degradation = pg6_gamma * y[13]
    return production - degradation


circuitODE = range(14)
circuitODE[0] = hexokinase
circuitODE[1] = glucose
circuitODE[2] = g6p
circuitODE[3] = phosphoglucoisomerase
circuitODE[4] = phosphoglucoisomerase2
circuitODE[5] = f6p
circuitODE[6] = phosphofructokinase
circuitODE[7] = f16p
circuitODE[8] = fructosephosphatealdolase
circuitODE[9] = fructosephosphatealdolase2
circuitODE[10] = g3p
circuitODE[11] = dhp
circuitODE[12] = g6pdehydrogenase
circuitODE[13] = pg6


t0 = 0.0
tmax = 1800.0
dt = 0.1
outfile = 'Glycosis.csv'
f = open(outfile, 'w')
header = ['time', 'hexokinase', 'glucose', 'g6p', 'phosphoglucoisomerase', 'phosphoglucoisomerase2', 'f6p', 'phosphofructokinase', 'f16p', 
          'fructosephosphatealdolase', 'fructosephosphatealdolase2', 'g3p', 'dhp', 'g6pdehydrogenase', 'pg6' ]
f.write(','.join(header) + '\n')
for x in ode.multirk4(circuitODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\n')
f.close()
