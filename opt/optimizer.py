#!/usr/bin/env python

'''
Simple resut plotter for bfpt-app

Usage:
  optimizer.py <path_to_f> -- [<extra_args>...]
'''

from docopt import docopt
from scipy.optimize import minimize
import subprocess
import math

po = docopt(__doc__)

path_to_f = po['<path_to_f>']
extra_args = po['<extra_args>']
print(f"[PO] {path_to_f=}")
print(f"[PO] {extra_args=}")

def f(eps, phi):
    print(f"f: {eps=:+14.8f}, {phi=:+14.8f} | ", end="")
    args = [f"-x {eps}", f"-y {phi}"]
    completed_process = subprocess.run(
        #['python3', 'x2y2.py', *args, *extra_args],
        ['path_to_f', *args, *extra_args],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    if completed_process.returncode != 0:
        print('ERROR as RC != 0|')
        return math.nan
    stdout_lines = completed_process.stdout.decode("utf-8").split('\n')
    optimizatoin_lines = [line for line in stdout_lines if line.startswith('[OPTIMIZATION]')]
    if len(optimizatoin_lines) != 1:
        print('ERROR as len(optimizatoin_lines) != 1|')
        return math.nan
    optimizatoin_line = optimizatoin_lines[0]
    optimizatoin_line_tokens = optimizatoin_line.split()
    if len(optimizatoin_line_tokens) < 2:
        print('ERROR as len(optimizatoin_line_tokens) < 2|')
        return math.nan
    if optimizatoin_line_tokens[1] != '[OK]':
        print('ERROR as optimizatoin_line_tokens[1] != [OK]|')
        return math.nan
    if len(optimizatoin_line_tokens) != 3:
        print('ERROR as len(optimizatoin_line_tokens) != 3|')
        return math.nan
    optimizatoin_line_value_token = optimizatoin_line_tokens[2]
    try:
        value = float(optimizatoin_line_value_token)
    except:
        print('ERROR value is not the float|')
        return math.nan
    print(f"OK {value=:+14.8f} |")
    return value

def f1(x):
    eps, phi = x
    return f(eps, phi)

x0 = (0.0, 0.0)
res = minimize(f1, x0, method='nelder-mead', options={'xatol': 1e-6, 'disp': True})
print(res.x)