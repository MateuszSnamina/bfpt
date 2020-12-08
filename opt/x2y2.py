#!/usr/bin/env python
'''
Simple resut plotter for bfpt-app

Usage:
  x2y2.py [-X <x_0>] [-Y <y_0>] -x <x> -y <y> 

Options:
  -X <x_0>, --x0 <x_0>
  -Y <y_0>, --y0 <y_0>
  -x <x>
  -y <y>      
'''

from docopt import docopt
#import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
#import numpy as np

po = docopt(__doc__)

x = po['-x']
y = po['-y']
x_0 = po['--x0']
y_0 = po['--y0']

try:
  x = float(x)
  y = float(y)
  x_0 = float(x_0) if x_0 else 0.0
  y_0 = float(y_0) if y_0 else 0.0
except Exception as e:
  message = "[OPTIMIZATION] [ERROR] Problem with parsing arguments: x_0, y_0, x, y. " + f"Exception details: {e}."
  print(message)
  exit(1)

print(f"{x=}")
print(f"{y=}")
print(f"{x_0=}")
print(f"{x_0=}")

result = (x - x_0)**2 + (y - y_0)**2

print(f"[OPTIMIZATION] [OK] {result}")