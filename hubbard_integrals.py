#!/usr/local/bin/python3
# Given a file in xyz format, U, and t, this script creates the HUBBARD integral file
# LV 2020

import re
import sys
import numpy
import math

fin = open(sys.argv[1], "r")

cgeom = []
for line in fin:
	if re.search("C", line):
		aux = line.split()
		cgeom.append(aux[1:4])
fin.close()

nsites = len(cgeom)

# Hubbard Hamiltonian parameters
tt = float(sys.argv[2])
uu = float(sys.argv[3])
t = numpy.zeros((nsites, nsites)) # hopping terms

# secondly t
ii = 0
for i in cgeom:
	x1 = float(i[0])
	y1 = float(i[1])
	z1 = float(i[2])

	jj = 0
	for j in cgeom:
		x2 = float(j[0])
		y2 = float(j[1])
		z2 = float(j[2])

		r12 = math.sqrt((x2  - x1) * (x2 - x1) +  (y2  - y1) * (y2 - y1) + (z2  - z1) * (z2 - z1))

		if r12 < 1.5 and ii != jj: # adjacent atoms
			t[ii, jj] = -tt

		jj = jj + 1
	ii = ii + 1

# create HUBBARD file
original_stdout = sys.stdout
fout = open("HUBBARD", "w")
sys.stdout = fout

# hopping terms
# alpha
for i in range(0, nsites):
	for j in range(i, nsites):
		if(abs(t[i, j]) > 1e-12):
			print("{:4.12e}   {:d}   {:d}   {:d}   {:d}".format(t[i, j], i+1, j+1, 0, 0))

print("{:4.12e}   {:d}   {:d}   {:d}   {:d}".format(0, 0, 0, 0, 0))

# beta
for i in range(0, nsites):
	for j in range(i, nsites):
		if(abs(t[i, j]) > 1e-12):
			print("{:4.12e}   {:d}   {:d}   {:d}   {:d}".format(t[i, j], i+1, j+1, 0, 0))

print("{:4.12e}   {:d}   {:d}   {:d}   {:d}".format(0, 0, 0, 0, 0))

# on-site repulsion
for i in range(0, nsites):
	print("{:4.12e}   {:d}   {:d}   {:d}   {:d}".format(0.5 * uu, i+1, i+1, i+1, i+1))

# nuclear repulsion plus core energy
print("{:4.12e}   {:d}   {:d}   {:d}   {:d}".format(0, 0, 0, 0, 0))

fin.close()
sys.stdout = original_stdout

