'''
ME Modeller: Generation of Ordinary Differential Equation Based
Python Simulation from Metabolic Pathway Definition.
 
Date created: 3rd September 2015
'''
 
import sys
import os
import itertools
 
def get_kinetics(data):
    kinetics = {}
    data = data.split('\n')
    data = [x.split(',') for x in data if len(x) > 0]
    for x in data:
        kinetics[x[0].strip()] = [k.strip() for k in x[1:]]
    return kinetics
     
def get_enzymes(data):
    data = data.split('\n')
    data = [x.split(',')[0].strip() for x in data if len(x) > 0]
    data = [x for x in data]
    data = list(set(data))
    data.sort()
    return data
 
def get_compounds(data):
    data = data.split('\n')
    data = [x.split(',')[1].strip() for x in data if len(x) > 0]
    reactants = [x.split('>')[0].strip() for x in data if len(x) > 0]
    reactants = [[r.strip() for r in x.split('+')] for x in reactants]
    reactants = list(itertools.chain(*reactants))
    products = [x.split('>')[1].strip() for x in data if len(x) > 0]
    products = [[r.strip() for r in x.split('+')] for x in products]
    products = list(itertools.chain(*products))
    data = reactants + products
    data = [x for x in data if len(x) > 0]
    data = list(set(data))
    data.sort()
    return data
 
def get_initial_conditions(initials, labeling):
    initial_conditions = {}
    initials = initials.split('\n')
    initials = [x for x in initials if len(x) > 0]
    initials = [[x.split('=')[0].strip(), x.split('=')[1].strip()]
                for x in initials]
    for x in initials:
        initial_conditions[x[0]] = x[1]
    for k in labeling.keys():
        if k not in initial_conditions:
            initial_conditions[k] = 0.0
    return initial_conditions
 
def get_compound_gammas(gammas, labeling):
    gamma_rates = {}
    gammas = gammas.split('\n')
    gammas = [x for x in gammas if len(x) > 0]
    gammas = [(x.split('=')[0].strip(), x.split('=')[1].strip())
                for x in gammas]
    for x in gammas:
        gamma_rates[x[0]] = x[1]
    for k in labeling.keys():
        if k not in gamma_rates:
            gamma_rates[k] = 0.0
    return gamma_rates
     
def generate_vectors(data):
    enzymes = get_enzymes(data)
    compounds = get_compounds(data)
    y = [0.0 for x in range(len(enzymes + compounds))]
    ODE = range(len(enzymes + compounds))
    labeling = {}
    index = 0
    for x in enzymes + compounds:
        labeling[x] = str(index)
        ODE[index] = x
        index = index + 1
    return (y, ODE, labeling)
 
def print_initials(y, initials, ODE, labeling, ofile):
    rlabel = {}
    for k in labeling.keys():
        rlabel[labeling[k]] = k
    k = rlabel.keys()
    k.sort()
    ofile.write('y = range(%s) \n' % len(y))
    for x in k:
        ofile.write('y[%s] = %s    # %s \n' % \
                    (x,
                     initials[rlabel[x]],
                     '_'.join(rlabel[x].split('.'))))
    ofile.write('\n')
    return ofile
 
def print_ODE_assignments(y, ODE, labeling, ofile):
    rlabel = {}
    for k in labeling.keys():
        rlabel[labeling[k]] = k
    k = rlabel.keys()
    k.sort()
    ofile.write('ODE = range(%s) \n' % len(ODE))
    for x in k:
        ofile.write('ODE[%s] = %s \n' % (x, rlabel[x]))
    ofile.write('\n')
    return ofile
 
def get_formation(data):
    # key in form is product; values are reactants
    # {'A': ['E1', 'B']} means 'A' is formed from 'B'
    # using enzyme E1
    form = {}
    data = data.split('\n')
    data = [x.split(',') for x in data if len(x) > 0]
    for x in data:
        enzyme = x[0].strip()
        reactants = x[1].strip()
        reactants = reactants.split('>')[0].strip()
        reactants = [r.strip() for r in reactants.split('+')]
        products = x[1].strip()
        products = products.split('>')[1].strip()
        products = [r.strip() for r in products.split('+')]
        for p in products:
            if p in form:
                t = form[p]
                t.append([enzyme] + reactants)
                form[p] = t
            else:
                form[p] = [[enzyme] + reactants]
    # print 'Form:', form
    return form
 
def get_usage(data):
    # key in form is product; values are reactants
    # {'A': ['E1', 'B']} means 'A' is used by 'B'
    # using enzyme E1
    use = {}
    data = data.split('\n')
    data = [x.split(',') for x in data if len(x) > 0]
    for x in data:
        enzyme = x[0].strip()
        reactants = x[1].strip()
        reactants = reactants.split('>')[0].strip()
        reactants = [r.strip() for r in reactants.split('+')]
        products = x[1].strip()
        products = products.split('>')[1].strip()
        products = [r.strip() for r in products.split('+')]
        for r in reactants:
            if r in use:
                t = use[r]
                t.append([enzyme] + products)
                use[r] = t
            else:
                use[r] = [[enzyme] + products]
    # print 'Use:', use
    return use
 
#def print_enzyme_
 
def print_enzyme_ODEs(data, kinetics, labeling, ofile):
    for enzyme in get_enzymes(data):
        ofile.write('''
def %s(t, y):
    production = %s
    degradation = %s * y[%s]
    return production - degradation \n\n''' \
                    % (enzyme,
                       0.0, #kinetics[enzyme][1],
                       kinetics[enzyme][0],
                       labeling[enzyme]))
    return ofile
     
def print_compound_ODEs(data, kinetics, gammas, labeling, ofile):
    formation = get_formation(data)
    usage = get_usage(data)
    for cpd in get_compounds(data):
        form_table = []
        use_table = []
        ofile.write('def %s(t, y): \n' % cpd)
        # Step 1: equations for formation
        if cpd in formation:
            for form in formation[cpd]:
                form_table.append('p' + form[1])
                ofile.write('\tp%s = (%s * y[%s] * (%s)) / ((%s) + %s) \n' \
                % (form[1],
                   kinetics[form[0]][2],
                   labeling[form[0]],
                   '*'.join(['y[%s]' % labeling[r] for r in form[1:]]),
                   '*'.join(['y[%s]' % labeling[r] for r in form[1:]]),
                   str(float(kinetics[form[0]][1])/float(kinetics[form[0]][3]))))
        # Step 2: equations for degradation
        ofile.write('\td%s = %s * y[%s] \n' \
              % (cpd,
                 gammas[cpd],
                 labeling[cpd]))
        # Step 3: equations for usage
        if cpd in usage:
            for use in usage[cpd]:
                use_table.append('u' + use[1])
                ofile.write('\tu%s = (%s * y[%s] * y[%s]) / (y[%s] + %s) \n' \
                % (use[1],
                   kinetics[use[0]][2],
                   labeling[use[0]],
                   labeling[cpd],
                   labeling[cpd],
                   str(float(kinetics[use[0]][1])/float(kinetics[use[0]][3]))))
        # Step 4: returning
        form_list = ' + '.join(form_table)
        use_list = ' + '.join(use_table)
        if (len(form_list) > 0) and (len(use_list) > 0):
            ofile.write('\treturn (%s) - d%s - (%s) \n\n' \
                  % (' + '.join(form_table),
                     cpd,
                     ' + '.join(use_table)))
        elif (len(form_list) > 0) and (len(use_list) == 0):
            ofile.write('\treturn (%s) - d%s \n\n' \
                  % (' + '.join(form_table),
                     cpd))
        elif (len(form_list) == 0) and (len(use_list) > 0):
            ofile.write('\treturn d%s - (%s) \n\n' \
                  % (cpd,
                     ' + '.join(use_table)))
    return ofile
 
def print_header(ofile):
    ofile.write('''
# This ODE simulation script is generated by ME Modeller (Generation of
# Ordinary Differential Equation Based Python Simulation from Metabolic
# Pathway Definition).
 
# --------- Start of ODE solver codes -------------------------------
def boundary_checker(y, boundary, type):
    for k in boundary.keys():
        if y[int(k)] < boundary[k][0] and type == 'lower':
            y[int(k)] = boundary[k][1]
        if y[int(k)] > boundary[k][0] and type == 'higher':
            y[int(k)] = boundary[k][1]
    return y
 
def RK4(funcs, x0, y0, step, xmax, 
             lower_bound=None, upper_bound=None):
    """
    Integrates a system of ODEs, y' = f(x, y), using fourth
    order Runge-Kutta method.
 
    @param funcs: system of differential equations
    @type funcs: list
    @param x0: initial value of x-axis, which is usually starting time
    @type x0: float
    @param y0: initial values for variables
    @type y0: list
    @param step: step size on the x-axis (also known as step in calculus)
    @type step: float
    @param xmax: maximum value of x-axis, which is usually ending time
    @type xmax: float
    """
    n = len(funcs)
    yield [x0] + y0
    f1, f2, f3, f4 = [0]*n, [0]*n, [0]*n, [0]*n
    max = 1e100
    while x0 < xmax:
        y1 = [0]*n
        for i in range(n):
            try: f1[i] = funcs[i](x0, y0)
            except TypeError: pass
            except ZeroDivisionError: f1[i] = max
            except OverflowError: f1[i] = max
        for j in range(n):
            y1[j] = y0[j] + (0.5*step*f1[j])
        for i in range(n):
            try: f2[i] = funcs[i]((x0+(0.5*step)), y1)
            except TypeError: pass
            except ZeroDivisionError: f2[i] = max
            except OverflowError: f2[i] = max
        for j in range(n):
            y1[j] = y0[j] + (0.5*step*f2[j])
        for i in range(n):
            try: f3[i] = funcs[i]((x0+(0.5*step)), y1)
            except TypeError: pass
            except ZeroDivisionError: f3[i] = max
            except OverflowError: f3[i] = max
        for j in range(n):
            y1[j] = y0[j] + (step*f3[j])
        for i in range(n):
            try: f4[i] = funcs[i]((x0+step), y1)
            except TypeError: pass
            except ZeroDivisionError: f4[i] = max
            except OverflowError: f4[i] = max
        x0 = x0 + step
        for i in range(n):
            y1[i] = y0[i] + (step * (f1[i] + (2.0*f2[i]) + (2.0*f3[i]) + f4[i]) / 6.0)
        if lower_bound: 
            y1 = boundary_checker(y1, lower_bound, 'lower')
        if upper_bound: 
            y1 = boundary_checker(y1, upper_bound, 'upper')
        y0 = y1
        yield [x0] + y1
# ----------- End of ODE solver codes -------------------------------
\n\n''')
    return ofile
     
def print_footer(start_time, end_time, timestep,
                 resultsfile, labeling, ofile):
    rlabel = {}
    for k in labeling.keys():
        rlabel[labeling[k]] = k
    k = rlabel.keys()
    k.sort()
    header = ['time'] + [rlabel[i] for i in k]
    ofile.write('''
t0 = %s
tmax = %s
dt = %s
outfile = '%s'
f = open(outfile, 'w')
header = %s
f.write(','.join(header) + '\\n')
for x in RK4(ODE, t0, y, dt, tmax):
    f.write(','.join([str(item) for item in x]) + '\\n')
f.close() \n''' % (start_time,
                   end_time,
                   timestep,
                   resultsfile,
                   str(header)))
    return ofile
 
def generate_model(modelfile, outputfile):
    splashscreen()
    print('Execution steps ...... \n')
    print('1. Read metabolic pathway definition file (%s)\n' % str(modelfile))
    exec('import %s as model' % str(modelfile))
    ofile = open(outputfile, 'w')
    print('2. Generate vector of enzymes and compounds')
    (y, ODE, labeling) = generate_vectors(model.reactions)
    print('   ... %s enzymes and compounds found in total\n' % str(len(y)))
    print('3. Processing enzymatic kinetics')
    kinetics = get_kinetics(model.enzyme_kinetics)
    print('   ... %s enzymes found in total\n' % str(len(kinetics)))
    print('4. Get initial conitions of enzymes and compounds\n')
    initials = get_initial_conditions(model.initial_conditions,
                                      labeling)
    print('5. Get degradation rates (gamma) of enzymes and compounds\n')
    gammas = get_compound_gammas(model.compound_gammas,
                                 labeling)
    print('6. Write ODE solver to generated script file (%s)\n' % str(outputfile))
    ofile = print_header(ofile)
    print('7. Write initial condition vector to generated script file (%s)\n' % str(outputfile))
    ofile = print_initials(y, initials, ODE, labeling, ofile)
    print('8. Write ODEs for enzymes to generated script file (%s)\n' % str(outputfile))
    ofile = print_enzyme_ODEs(model.reactions,
                              kinetics,
                              labeling,
                              ofile)
    print('9. Write ODEs for compounds to generated script file (%s)\n' % str(outputfile))
    ofile = print_compound_ODEs(model.reactions,
                                kinetics,
                                gammas,
                                labeling,
                                ofile)
    print('10. Write ODEs assignments to generated script file (%s)\n' % str(outputfile))
    ofile = print_ODE_assignments(y, ODE, labeling, ofile)
    print('11. Write ODEs executor to generated script file (%s)\n' % str(outputfile))
    ofile = print_footer(model.start_time,
                         model.end_time,
                         model.timestep,
                         model.resultsfile,
                         labeling,
                         ofile)
    print('Model generation complete. \n')
    splashscreen()
    ofile.close()
 
def splashscreen():
    print('''
 
ME Modeller: Generation of Ordinary Differential Equation Based
             Python Simulation from Metabolic Pathway Definition
 
Copyright (c) 2015, Maurice HT Ling (mauriceling@acm.org)
All Rights Reserved
 
Please type 'python me_modeller.py help' for further help.
 
''')
     
def usage():
    print('''
Usage: python me_modeller.py <model definition file> <model output file>
 
where <model definition file> is the metabolic or pathway definition to be
used as input to generate the <model output file> as an executable Python
module.
 
For example, "python me_modeller.py riboflavin.py ribo_script.py" will
take the metabolic pathway definition from riboflavin.py file to generate
the corresponding Python code representation of the ODEs as ribo_script.py
file.
 
How to write the metabolic pathway definition file?
<See when I have enough energy for this>
''')
     
if __name__=='__main__':
    if len(sys.argv) < 3:
        usage()
    elif sys.argv[1] == 'help':
        usage()
    else:
        generate_model(sys.argv[1].split('.')[0].strip(), sys.argv[2])