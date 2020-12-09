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

f_global_id_counter = 0
def f(eps, phi):
    global f_global_id_counter
    f_global_id_counter = f_global_id_counter + 1
    log_tag = f"[FUNCTION-CALL] [#{f_global_id_counter:05d}] "
    #print(log_tag + f"[ARGS] {eps=:+14.8f}, {phi=:+14.8f}")
    args = ["--hamiltonian_eps", f"{eps}", "--hamiltonian_phi", f"{phi}"]
    cli_list = [path_to_f, *args, *extra_args]
    #print(log_tag + f"[CLI] {' '.join(cli_list)}")
    completed_process = subprocess.run(cli_list, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #print(log_tag + "[ERROR] as RC != 0")
    #print(log_tag + f"[NOTE] {completed_process.returncode=}")
    #stdout_lines = completed_process.stdout.decode("utf-8").split('\n')
    #stderr_lines = completed_process.stderr.decode("utf-8").split('\n')
    #stdout_lines = [log_tag + "[NOTE] [STDOUT] " + line for line in stdout_lines]
    #stderr_lines = [log_tag + "[NOTE] [STDERR] " + line for line in stderr_lines]
    #print('\n'.join(stdout_lines))
    #print('\n'.join(stderr_lines))
    if completed_process.returncode != 0:
        print(log_tag + "[ERROR] as RC != 0")
        print(log_tag + f"[NOTE] [DETAILS] {completed_process.returncode=}")
        print(log_tag + f"[NOTE] [ARGS] {eps=:+14.8f}, {phi=:+14.8f}")
        print(log_tag + f"[NOTE] [CLI] {' '.join(cli_list)}")
        stdout_lines = completed_process.stdout.decode("utf-8").split('\n')
        stderr_lines = completed_process.stderr.decode("utf-8").split('\n')
        stdout_lines = [log_tag + "[NOTE] [STDOUT] " + line for line in stdout_lines]
        stderr_lines = [log_tag + "[NOTE] [STDERR] " + line for line in stderr_lines]
        print('\n'.join(stdout_lines))
        print('\n'.join(stderr_lines))
        raise RuntimeError("Error as RC != 0.")
        #return math.nan
    stdout_lines = completed_process.stdout.decode("utf-8").split('\n')
    optimizatoin_lines = [line for line in stdout_lines if line.startswith('[OPTIMIZATION]')]
    if len(optimizatoin_lines) != 1:
        print(log_tag + "[ERROR] as len(optimizatoin_lines) != 1")
        print(log_tag + f"[NOTE] [DETAILS] {len(optimizatoin_lines)=}")
        print(log_tag + f"[NOTE] [ARGS] {eps=:+14.8f}, {phi=:+14.8f}")
        print(log_tag + f"[NOTE] [CLI] {' '.join(cli_list)}")
        raise RuntimeError("Error as len(optimizatoin_lines) != 1.")
        #return math.nan
    optimizatoin_line = optimizatoin_lines[0]
    optimizatoin_line_tokens = optimizatoin_line.split()
    if len(optimizatoin_line_tokens) < 2:
        print(log_tag + "[ERROR] as len(optimizatoin_line_tokens) < 2")
        print(log_tag + f"[NOTE] [DETAILS] {len(optimizatoin_line_tokens)=}")
        print(log_tag + f"[NOTE] [ARGS] {eps=:+14.8f}, {phi=:+14.8f}")
        print(log_tag + f"[NOTE] [CLI] {' '.join(cli_list)}")
        raise RuntimeError("Error as len(optimizatoin_line_tokens) < 2.")
        #return math.nan
    if optimizatoin_line_tokens[1] != '[OK]':
        print(log_tag + "[ERROR] as optimizatoin_line_tokens[1] != [OK]")
        print(log_tag + f"[NOTE] [DETAILS] {optimizatoin_line_tokens[1]=}")
        print(log_tag + f"[NOTE] [ARGS] {eps=:+14.8f}, {phi=:+14.8f}")
        print(log_tag + f"[NOTE] [CLI] {' '.join(cli_list)}")
        raise RuntimeError("Error as optimizatoin_line_tokens[1] != [OK].")
        #return math.nan
    if len(optimizatoin_line_tokens) != 3:
        print(log_tag + "[ERROR] as len(optimizatoin_line_tokens) != 3")
        print(log_tag + f"[NOTE] [DETAILS] {len(optimizatoin_line_tokens)=}")
        print(log_tag + f"[NOTE] [ARGS] {eps=:+14.8f}, {phi=:+14.8f}")
        print(log_tag + f"[NOTE] [CLI] {' '.join(cli_list)}")
        raise RuntimeError("Error as len(optimizatoin_line_tokens) != 3.")
        #return math.nan
    try:
        value = float(optimizatoin_line_tokens[2])
    except:
        print(log_tag + "[ERROR] value string does not represent a float number")
        print(log_tag + f"[NOTE] [DETAILS] {optimizatoin_line_tokens[2]=}")
        print(log_tag + f"[NOTE] [ARGS] {eps=:+14.8f}, {phi=:+14.8f}")
        print(log_tag + f"[NOTE] [CLI] {' '.join(cli_list)}")
        raise RuntimeError("Value string does not represent a float number.")
        #return math.nan
    print(log_tag + f"[OK] function args: {eps=:+14.8f}, {phi=:+14.8f} | function value: {value=:+14.8f}")
    return value

def f1(x):
    eps, phi = x
    return f(eps, phi)

x0 = (0.0, 0.0)
res = minimize(f1, x0, method='nelder-mead', options={'xatol': 1e-6, 'disp': True})
print(res.x)