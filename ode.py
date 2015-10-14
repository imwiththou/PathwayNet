import sys

def rk4(func, x, y, step, xmax):
    """
    Integrates y'=f(x,y) using 4th step-order Runge-Kutta.
    
    @param func: a differential equation
    @type func: list
    @param x: initial value of x-axis, which is usually starting time
    @type x: float
    @param y: initial value for y-axis
    @type y: float
    @param step: step size on the x-axis (also known as step in calculus)
    @type step: float
    @param xmax: maximum value of x-axis, which is usually ending time
    @type xmax: float
    """
    yield [x, y]
    while x < xmax:
        f1 = func(x, y)
        f2 = func(x+0.5*step, y+0.5*step*f1)
        f3 = func(x+0.5*step, y+0.5*step*f2)
        f4 = func(x+step, y+step*f3)
        x = x + step
        y = y + step*(f1+2.0*f2+2.0*f3+f4)/6.0
        yield [x, y]

def boundary_checker(y, boundary, type):
    for k in boundary.keys():
        if y[int(k)] < boundary[k][0] and type == 'lower':
            y[int(k)] = boundary[k][1]
        if y[int(k)] > boundary[k][0] and type == 'higher':
            y[int(k)] = boundary[k][1]
    return y

def multirk4(funcs, x0, y0, step, xmax, 
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
            y1[i] = y0[i] + (step * \
                (f1[i] + (2.0*f2[i]) + (2.0*f3[i]) + f4[i]) / 6.0)
        if lower_bound: 
            y1 = boundary_checker(y1, lower_bound, 'lower')
        if upper_bound: 
            y1 = boundary_checker(y1, upper_bound, 'upper')
        y0 = y1
        yield [x0] + y1
