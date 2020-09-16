#!/usr/bin/env python
import sys
import yaml
import math

import matplotlib.pyplot as plt
import utility.multiple_formatter

n_sites=20

def get_value_if_prefix_matches(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix):]
    return None

def patch_data_for_two_pi(data):
    idx_of_zero = data['domain'].index(0.0)
    data['domain'].append(2 * math.pi)
    data['es_excitation_energies'].append(data['es_excitation_energies'][idx_of_zero])
    data['es_absolute_energies'].append(data['es_absolute_energies'][idx_of_zero])

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
    data = {
        'gs_energy': gs_energy,
        'domain': domain,
        'es_excitation_energies' : es_excitation_energies,
        'es_absolute_energies' : es_absolute_energies}
    patch_data_for_two_pi(data)
    return data

def prepare_reference_data(n_sites):
    GRID = 100
    gs_energy = - n_sites * (math.log(2) - 0.25)
    domain = [2 * math.pi * x / GRID for x in range(0, GRID)]
    es_excitation_energies = [math.pi / 2 * abs(math.sin(k)) for k in domain]
    es_absolute_energies = [gs_energy + es_exciation_enery for es_exciation_enery in es_excitation_energies]
    data = {
        'gs_energy': gs_energy,
        'domain': domain,
        'es_excitation_energies' : es_excitation_energies,
        'es_absolute_energies' : es_absolute_energies}
    patch_data_for_two_pi(data)
    return data

def put_data_on_absolute_energy_plot(ax, data, **style):
  ax.axhline(y = data['gs_energy'], **style)
  ax.plot(data['domain'], data['es_absolute_energies'], **style)

def put_data_on_excitation_energy_plot(ax, data, **style):
  ax.plot(data['domain'], data['es_excitation_energies'], **style)

def index_to_color(idx, length):
    if len(data_list) > 1:
        λ = index / (len(data_list) - 1)
        color = (1 - λ, 0, λ, 0.7)
    else:
        color = 'blue' 
    return color

#################################################
## Prepare data                                ##
#################################################

data_list = []
print(f"{sys.argv=}")
for log_file_path in  sys.argv[1:]:
  print(f"{log_file_path=}")
  data = grep_log_data(log_file_path)
  data_list.append(data)

reference_data = prepare_reference_data(n_sites)

#################################################
## Build the plot                              ##
#################################################

fig = plt.figure()
ax1 = fig.add_subplot(111)
put_data_on_absolute_energy_plot(ax1, reference_data, color='gray', linewidth=4)
for (index , data) in enumerate(data_list):
    color = index_to_color(index, len(data_list))  
    put_data_on_absolute_energy_plot(ax1, data, color=color, linewidth=2)
ax1.grid()
ax1.set_xlim(0, 2 * math.pi)
ax1.xaxis.set_major_locator(plt.MultipleLocator(math.pi / 2))
ax1.xaxis.set_minor_locator(plt.MultipleLocator(math.pi / 12))
ax1.xaxis.set_major_formatter(plt.FuncFormatter(utility.multiple_formatter.multiple_formatter()))
ax1.set_xlabel('k', color='black')
ax1.set_ylabel('absolute enery', color='black')
plt.savefig('absolute_enery.png')
plt.show()

fig = plt.figure()
ax1 = fig.add_subplot(111)
put_data_on_excitation_energy_plot(ax1, reference_data, color='gray', linewidth=4)
for (index , data) in enumerate(data_list):
    color = index_to_color(index, len(data_list))  
    put_data_on_excitation_energy_plot(ax1, data, color=color, linewidth=2)
ax1.grid()
ax1.set_xlim(0, 2 * math.pi)
ax1.xaxis.set_major_locator(plt.MultipleLocator(math.pi / 2))
ax1.xaxis.set_minor_locator(plt.MultipleLocator(math.pi / 12))
ax1.xaxis.set_major_formatter(plt.FuncFormatter(utility.multiple_formatter.multiple_formatter()))
ax1.set_xlabel('k', color='black')
ax1.set_ylabel('excitation enery', color='black')
plt.savefig('excitation_enery.png')
plt.show()
