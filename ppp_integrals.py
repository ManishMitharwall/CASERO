#!/usr/local/bin/python3
# Given a file in xyz format, this script creates a FCIDUMP with integrals corresponding to the
# Parisier-Parr-Pople (PPP) Hamiltonian according to the parametrization from JCP 134, 024114 (2011)
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

# PPP Hamiltonian parameters
t = numpy.zeros((nsites, nsites)) # one-electron resonance integrals
g = numpy.zeros((nsites, nsites)) # two-electron repulsion integrals

# firstly gamma, as it is needed for t
ii = 0
for i in cgeom:
	x1 = float(i[0])
	y1 = float(i[1])
	z1 = float(i[2])

	g[ii, ii] = 11.56718

	jj = 0
	for j in cgeom:
		x2 = float(j[0])
		y2 = float(j[1])
		z2 = float(j[2])

		r12 = math.sqrt((x2  - x1) * (x2 - x1) +  (y2  - y1) * (y2 - y1) + (z2  - z1) * (z2 - z1))

		if ii != jj:
			g[ii, jj] = 21.88221 / (math.sqrt(r12 * r12 + 7.05909)) 

		if r12 < 1.5 and ii != jj: # adjacent atoms
			g[ii, ii] = g[ii, ii] - 0.16640 / ((r12 - 2.30448) *(r12 - 2.30448)) 

		jj = jj + 1
	ii = ii + 1

# secondly t
ii = 0
for i in cgeom:
	x1 = float(i[0])
	y1 = float(i[1])
	z1 = float(i[2])

	t[ii, ii] = -9.13597

	jj = 0
	for j in cgeom:
		x2 = float(j[0])
		y2 = float(j[1])
		z2 = float(j[2])

		r12 = math.sqrt((x2  - x1) * (x2 - x1) +  (y2  - y1) * (y2 - y1) + (z2  - z1) * (z2 - z1))

		if ii != jj:
			t[ii, ii] = t[ii, ii] - g[ii,jj] * r12 / (math.sqrt(r12 * r12 - 0.15076))

		if r12 < 1.5 and ii != jj: # adjacent atoms
			t[ii, jj] = -28.07749 * math.exp(-1.65878 * r12)

		jj = jj + 1
	ii = ii + 1

# create FCIDUMP file
original_stdout = sys.stdout
fout = open("FCIDUMP", "w")
sys.stdout = fout

# header
print("&FCI NORB=  {:d},NELEC=  {:d},MS2= 0,".format(nsites, nsites))
aux = ""
for i in range(0, nsites):
	aux = aux + "1,"
print(" ORBSYM={:s}".format(aux))
print("  ISYM=1,")
print(" /")

# two-electron integrals
for i in range(0, nsites):
	for j in range(i, nsites):
		if(abs(g[i, j]) > 1e-12):
			print("{:4.12e}   {:d}   {:d}   {:d}   {:d}".format(g[i, j], i+1, i+1, j+1, j+1))

# one-electron integrals
for i in range(0, nsites):
	for j in range(i, nsites):
		if(abs(t[i, j]) > 1e-12):
			print("{:4.12e}   {:d}   {:d}   {:d}   {:d}".format(t[i, j], i+1, j+1, 0, 0))

# nuclear repulsion plus core energy
print("{:4.12e}   {:d}   {:d}   {:d}   {:d}".format(0, 0, 0, 0, 0))

fin.close()
sys.stdout = original_stdout

