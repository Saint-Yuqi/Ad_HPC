#!/usr/bin/env python
import numpy as np
from sys import argv,exit

header_type = np.dtype([
  ('time', '>f8'),('N', '>i4'), ('Dims', '>i4'),
  ('Ngas', '>i4'), ('Ndark', '>i4'), ('Nstar', '>i4'), 
  ('pad', '>i4')])
gas_type  = np.dtype([
  ('mass','>f4'), 
  ('r', '>f4', (3,)),  ('v', '>f4', (3,)),
  ('rho','>f4'), ('temp','>f4'), ('hsmooth','>f4'), 
  ('metals','>f4'), ('phi','>f4')])

dark_type = np.dtype([
  ('mass', '>f4'), 
  ('r', '>f4', (3,)), ('v', '>f4', (3,)),
  ('eps', '>f4'), ('phi', '>f4')
])

star_type  = np.dtype([
  ('mass','>f4'), 
  ('r', '>f4', (3,)),  ('v', '>f4', (3,)),
  ('metals','>f4'), ('tform','>f4'), 
  ('eps','>f4'), ('phi','>f4')])

with open(argv[1],'rb') as tipsy:
  header = np.fromfile(tipsy,dtype=header_type,count=1)
  header=header[0] # there is only a single row
  gas  = np.fromfile(tipsy,dtype=gas_type,count=header['Ngas'])
  dark = np.fromfile(tipsy,dtype=dark_type,count=header['Ndark'])
  star = np.fromfile(tipsy,dtype=star_type,count=header['Nstar'])

# Extract positions
x_coords = dark['r'][:, 0]
y_coords = dark['r'][:, 1]
z_coords = dark['r'][:, 2]

# Compute min and max for each coordinate
x_min, x_max = np.min(x_coords), np.max(x_coords)
y_min, y_max = np.min(y_coords), np.max(y_coords)
z_min, z_max = np.min(z_coords), np.max(z_coords)

# Compute total mass
total_mass = np.sum(dark['mass'])

# Print results
print(f"Minimum x, y, z coordinates: ({x_min}, {y_min}, {z_min})")
print(f"Maximum x, y, z coordinates: ({x_max}, {y_max}, {z_max})")
print(f"Total mass in the box: {total_mass}")
print(header)
print(dark)
print(dark['r'])
print(dark['r'][0])
print(dark['r'][:,0])