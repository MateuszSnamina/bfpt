#!/usr/bin/env python
'''
Simple resut plotter for bfpt-app

Usage:
  plot_results.py [-4] [-2] [-p] [-s] af [-n <n_sites>] [-J <coupling-value>] <log_file>...
  plot_results.py [-4] [-2] [-p] [-s] fm [-n <n_sites>] [-J <coupling-value>] <log_file>...
  plot_results.py [-4] [-2] [-p] [-s] <log_file>...

Options:
  -s, --save                             Save the figures as png files in cwd
  -n <n_sites>, --n_sites <n_sites>      n_sites of the chain [default: 8]
  -J <value>, --coupling-value <value>   Ising coupling value [default: 1.0]
  -4, --reconstruct4                     Uses data from [0, Pi/2] to infer data for (Pi/2, Pi]
  -2, --reconstruct2                     Uses data from [0, Pi] to infer data for (Pi, 2Pi]
  -p, --patch_two_pi                     Fill value for 2Pi as it is for 0
'''

import yaml
import math

from docopt import docopt
import matplotlib.pyplot as plt
from utility.get_value_if_prefix_matches import get_value_if_prefix_matches
from utility.multiple_formatter import multiple_formatter

###############################################
## reconstruct data functions                ##
###############################################

def is_near(v1, v2, epsilon=0.00001):
  return abs(v1 - v2) < epsilon

def reconstruct_pivot(pivot, X, Y1, Y2):
  assert(len(X) == len(Y1) == len(Y2))
  assert(X == sorted(X))
  assert(len(X) >= 1)
  if is_near(X[-1], pivot) and len(X) >= 2:
    X = X + [2 * pivot - x for x in X[-2::-1]]
    Y1 = Y1 + Y1[-2::-1]
    Y2 = Y2 + Y2[-2::-1]
  else:
    X = X + [2 * pivot - x for x in X[::-1]]
    Y1 = Y1 + Y1[::-1]
    Y2 = Y2 + Y2[::-1]
  return X, Y1, Y2

def reconstruct2(X, Y1, Y2):
  return reconstruct_pivot(math.pi, X, Y1, Y2)

def reconstruct4(X, Y1, Y2):
  return reconstruct_pivot(math.pi / 2, X, Y1, Y2)

def patch_two_pi(X, Y1, Y2):
  idx_of_zero = X.index(0.0)
  X.append(2 * math.pi)
  Y1.append(Y1[idx_of_zero])
  Y2.append(Y2[idx_of_zero])
  return X, Y1, Y2


###############################################
## build data functions                      ##
###############################################

def grep_log_data(log_file_path):
    with open(log_file_path) as log_file:
        for line in log_file:
            line = line.strip()
            prefix_gs_energy = "[RESULT] [POST] gs_energy: "
            prefix_domain = "[RESULT] [POST] domain: "            
            prefix_es_absolute_energies = "[RESULT] [POST] es_absolute_energies: "
            prefix_es_excitation_energies = "[RESULT] [POST] es_excitation_energies: " 
            if (_ := get_value_if_prefix_matches(line, prefix_gs_energy)):
                gs_energy = yaml.safe_load(_)
            if (_ := get_value_if_prefix_matches(line, prefix_domain)):
                domain = yaml.safe_load(_)
            if (_ := get_value_if_prefix_matches(line, prefix_es_absolute_energies)):
                es_absolute_energies = yaml.safe_load(_)
            if (_ := get_value_if_prefix_matches(line, prefix_es_excitation_energies)):
                es_excitation_energies = yaml.safe_load(_)
    if po['--reconstruct4']:
       domain, es_absolute_energies, es_excitation_energies = reconstruct4(domain, es_absolute_energies, es_excitation_energies)
    if po['--reconstruct2']:
       domain, es_absolute_energies, es_excitation_energies = reconstruct2(domain, es_absolute_energies, es_excitation_energies)
    if po['--patch_two_pi']:
       domain, es_absolute_energies, es_excitation_energies = patch_two_pi(domain, es_absolute_energies, es_excitation_energies)
    data = {
        'gs_energy': gs_energy,
        'domain': domain,
        'es_excitation_energies' : es_excitation_energies,
        'es_absolute_energies' : es_absolute_energies}
    return data

def prepare_af_reference_data(n_sites, J):
    GRID = 100
    gs_energy = - J * n_sites * (math.log(2) - 0.25)
    domain = [2 * math.pi * x / GRID for x in range(0, GRID)]
    es_excitation_energies =  [J * math.pi / 2 * abs(math.sin(k)) for k in domain]
    es_absolute_energies =  [J * gs_energy + es_exciation_enery for es_exciation_enery in es_excitation_energies]
    domain, es_absolute_energies, es_excitation_energies = patch_two_pi(domain, es_absolute_energies, es_excitation_energies)
    data = {
        'gs_energy': gs_energy,
        'domain': domain,
        'es_excitation_energies' : es_excitation_energies,
        'es_absolute_energies' : es_absolute_energies}
    return data

def prepare_fm_reference_data(n_sites, J):
    assert(False) # TODO implement
    return None

###############################################
## put on plot functions                     ##
###############################################

def put_data_on_absolute_energy_plot(ax, data, scatter=False, **style):
  ax.axhline(y = data['gs_energy'], **style)
  ax.plot(data['domain'], data['es_absolute_energies'], **style)
  if scatter:
      ax.scatter(data['domain'], data['es_absolute_energies'], marker='*', **style)

def put_data_on_excitation_energy_plot(ax, data, scatter=False, **style):
  ax.plot(data['domain'], data['es_excitation_energies'], **style)
  if scatter:
      ax.scatter(data['domain'], data['es_excitation_energies'], marker='*', **style)

def index_to_color(idx, length):
    if length > 1:
        λ = index / (length - 1)
        color = (1 - λ, 0, λ, 0.7)
    else:
        color = 'blue' 
    return color

#################################################
## Prepare data                                ##
#################################################

po = docopt(__doc__)

if po['af']:
  n_sites = int(po['--n_sites'])
  J = float(po['--coupling-value'])
  reference_data = prepare_af_reference_data(n_sites, J)
elif po['fm']:
  n_sites = int(po['--n_sites'])
  J = float(po['--coupling-value'])
  reference_data = prepare_fm_reference_data(n_sites, J)
else:
  reference_data = None

data_list = []
for log_file_path in po['<log_file>']:
  print(f"Grep {log_file_path=} ...")
  data = grep_log_data(log_file_path)
  data_list.append(data)


#################################################
## Build the plot                              ##
#################################################

fig = plt.figure()
ax1 = fig.add_subplot(111)
if reference_data:
    put_data_on_absolute_energy_plot(ax1, reference_data, color='gray', linewidth=4)
for (index , data) in enumerate(data_list):
    color = index_to_color(index, len(data_list))  
    put_data_on_absolute_energy_plot(ax1, data, scatter=True, color=color, linewidth=2)
ax1.grid()
ax1.set_xlim(0, 2 * math.pi)
ax1.xaxis.set_major_locator(plt.MultipleLocator(math.pi / 2))
ax1.xaxis.set_minor_locator(plt.MultipleLocator(math.pi / 12))
ax1.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
ax1.set_xlabel('k', color='black')
ax1.set_ylabel('absolute enery', color='black')
if po['--save']:
  plt.savefig('absolute_enery.png')
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(111)
if reference_data:
    put_data_on_excitation_energy_plot(ax1, reference_data, color='gray', linewidth=4)
for (index , data) in enumerate(data_list):
    color = index_to_color(index, len(data_list))  
    put_data_on_excitation_energy_plot(ax1, data, scatter=True, color=color, linewidth=2)
ax1.grid()
ax1.set_xlim(0, 2 * math.pi)
ax1.set_ylim(0, None)
ax1.xaxis.set_major_locator(plt.MultipleLocator(math.pi / 2))
ax1.xaxis.set_minor_locator(plt.MultipleLocator(math.pi / 12))
ax1.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
ax1.set_xlabel('k', color='black')
ax1.set_ylabel('excitation enery', color='black')
if po['--save']:
  plt.savefig('excitation_enery.png')
plt.show()
