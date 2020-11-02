"""

This file loads the data generated from the simulations and 
generates the engram-related figures.

usage:
python make_figs.py


"""
import numpy as np
from numpy.linalg import norm

import matplotlib.pyplot as plt
import random
import pandas as pd
import sys
import scipy.stats as stats

import matplotlib as mpl
mpl.rcParams["errorbar.capsize"] = 2
mpl.rcParams["lines.linewidth"] = 1
mpl.rcParams['pdf.fonttype'] = 42



def trevrolls(frates):
	r2 = 0.
	rs = 0.
	n = float(frates.size)
	for i in range(frates.size): # np.nditer(frates):
		r2 += (frates[i]**2)/n
		rs += frates[i]/n
	return 1. - ((rs**2)/r2)





np.set_printoptions( threshold=999999999999999)


NPYRS = 400
NINH = 100
NRUNS=10
CLUSTERED=0

def scos(a,b):
	 return np.dot(a, b)/(norm(a) * norm(b))




for SUFF in ('5', '40'):

	plt.figure()

	pts = np.loadtxt("pyr_active%sC.txt"%(SUFF));
	plt.plot(pts/4., label="caap");

	pts = np.loadtxt("pyr_active%sN.txt"%(SUFF));
	plt.plot(pts/4., label="norm");

	plt.legend()

	plt.xlabel('#memory')
	plt.ylabel('% Active Neurons (Engram size)')


	pts = np.loadtxt("patterns%sN.txt"%(SUFF));
	sims = [];
	pairs = [];

	simX  = [];
	simYC  = [];
	simYN  = [];

	spC  = [];
	spN  = [];

	plt.figure()
	plt.imshow(pts)

	spikesC = np.loadtxt("test_spikes%sC.txt"%(SUFF));
	spikesN = np.loadtxt("test_spikes%sN.txt"%(SUFF));

	plt.figure();
	plt.imshow(spikesC, interpolation='none', aspect='auto');
	plt.title('Spikes Caap %s'%(SUFF));
	plt.colorbar()

	plt.figure();
	plt.imshow(spikesN, interpolation='none', aspect='auto');
	plt.title('Spikes Norm %s'%(SUFF));
	plt.colorbar()



	for i in range(0, len(pts)):

		spC.append(trevrolls(spikesC[i, 0:400]))
		spN.append(trevrolls(spikesN[i, 0:400]))

		for j in range(0, i):
			print(i,j);
			#pairs.append(i,j);
			
			simX.append( scos( pts[i], pts[j] ))

			a = spikesC[i, 0:400]
			b = spikesC[j, 0:400]

			simYC.append( scos( a, b ))

			a = spikesN[i, 0:400]
			b = spikesN[j, 0:400]

			simYN.append( scos( a, b ))


	plt.figure()
	plt.plot( spC);
	plt.plot( spN);
	plt.ylabel('Sparsity ')
	plt.title('Sparsity (Treves-Rolls) %s'%(SUFF));
	plt.xlabel('# memory');

	plt.figure()
	plt.hist( simX, 100)
	plt.title('Input-output similarity');
	plt.ylabel('similarity');
	plt.xlabel('# pairs');

	plt.figure();
	plt.scatter(simX, simYC, label='caap');
	plt.scatter(simX, simYN, label='norm');
	plt.legend()
	


plt.show();



