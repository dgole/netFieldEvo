#!/usr/bin/python
from __future__ import unicode_literals
import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import os
import math
from scipy import fftpack as fft
import sys
import netFieldEvoAnalysis as reader
import resource
from matplotlib.backends.backend_pdf import PdfPages
import time
import matplotlib.animation as animation
m.rcParams['text.usetex'] = True
m.rcParams['text.latex.unicode'] = True

######################
#self.header = [r"$\alpha$", r"$h/r$", r"$T_c^4$", r"$T_disk^4$", r"$c_s$", r"$\rho$", r"$\kappa_R$", r"$\nu$", r"$\tau$", r"$mdot$", r"$v_{adv}$", r"$v_{diff}$", r"$\beta$", r"$B_z$", r"$B_{rs}$", r"$\psi$", r"$\Sigma$"  		
		
# 0  alpha
# 1  h/r
# 2  Tc4
# 3  Tdisk4
# 4  cs
# 5  rho
# 6  kappaR
# 7  nu
# 8  tau
# 9  mdot
# 10 vadv
# 11 vdiff
# 12 Bz
# 13 Brs
# 14 Psi
# 15 Sigma
# 16 beta
# 17 Bz/Brs
# 18 Pr eff
# 19 tadv 
# 20 tdiff
# 21 dFlux

# avaliable routines
# profile(self, col, time, vLineCoords=[1.0, 50.0], hLineCoords=[], logOption=0, save=None, savePath=None, legendLabel=None, ymin=None, ymax=None)
# multiProfile(self, col, nProfiles, spacing="log", vLineCoords1=[1.0,50.0], hLineCoords1=[], logOption=0, save=None, savePath=None, ymin=None, ymax=None)
# makeMultiAnim(self, timeCutFactor=10, lengthInSeconds=20, savePath=None, show=None, save=None)
# timeEvo(self, col, r, logOption=0, save=None, savePath=None, legendLabel=None, ymin=None, ymax=None)
# stPlot(self, col, cmapType="viridis", logOption=0, save=None, savePath=None, clim1=None, clim2=None)

for idNum in sys.argv[1:]:
	# make data object
	do = reader.Data("../output/run" + str(idNum) + "/")

	for i in range(do.nr):
		print(i, do.r[i], do.data[20][0,i], do.data[19][0,i])

	# all time scales on the same plot
	n=0
	plt.loglog(do.r, do.tOrbit, label="tOrbit")
	plt.loglog(do.r, do.data[19][n], label="tAdv")
	plt.loglog(do.r, do.data[20][n], label="tDiff")
	plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
	plt.ylim(1.e-4,1.e5)
	plt.ylabel("timescale (code units)")
	plt.xlabel("r (code units)")
	plt.savefig(do.pdfName, format='pdf', bbox_inches='tight')
	plt.savefig(do.savePath + "timeScales.png", bbox_inches='tight')
	plt.clf()

	# space time plot
	do.stPlot(12, cmapType="coolwarm", logOption=2, save="pdf")
	do.stPlot(0, cmapType="viridis", logOption=1, save="pdf")
	do.stPlot(9, cmapType="viridis", logOption=1, save="pdf", clim1=1.e8, clim2=1.e10)

	do.stPlot(12, cmapType="coolwarm", logOption=2, save="png", xLabel=0)
	do.stPlot(0, cmapType="viridis", logOption=1, save="png", xLabel=1)
	do.stPlot(9, cmapType="viridis", logOption=1, save="png", clim1=1.e8, clim2=1.e10, xLabel=1)

	# sigma multi profile
	do.multiProfile(15, 10, logOption=1, save="pdf", spacing="log", vLineCoords1=[1.0,10.0,100.0])
	# alpha multi profile
	do.multiProfile(0, 10, logOption=1, save="pdf", spacing="log", vLineCoords1=[1.0,10.0,100.0])
	# h/r multi profile
	do.multiProfile(1, 10, logOption=1, ymin=1.e-2, ymax=1.e0, save="pdf", spacing="log", vLineCoords1=[1.0,10.0,100.0])
	# nu multi profile
	do.multiProfile(7, 10, logOption=1, save="pdf", spacing="log", vLineCoords1=[1.0,10.0,100.0])
	# mdot multi profile
	do.multiProfile(9, 10, logOption=1, save="pdf", spacing="log", vLineCoords1=[1.0,10.0,100.0])
	# bz multi profile
	do.multiProfile(12, 10, logOption=0, save="pdf", spacing="log", vLineCoords1=[1.0,10.0,100.0])
	do.multiProfile(12, 10, logOption=1, save="pdf", vLineCoords1=[1.0,10.0,100.0])
	# beta multi profile
	do.multiProfile(16, 10, logOption=1, save="pdf", vLineCoords1=[1.0,10.0,100.0])

	# brs multi profile
	#do.multiProfile(6, 10, logOption=0, save="pdf", ymin=-1.0, ymax=1.0)
	#do.multiProfile(6, 10, logOption=1, save="pdf", ymin=1.e-3, ymax=1.e3)
	# psi multi profile
	do.multiProfile(14, 10, logOption=1, spacing="log", save="pdf", vLineCoords1=[1.0,10.0,100.0])
	# brs/bz multi profile
	do.multiProfile(17, 10, logOption=0, save="pdf", ymin=-2.0, ymax=2.0, vLineCoords1=[1.0,10.0,100.0])




	do.profile(0, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
	do.profile(1, 0.01, vLineCoords=[10.0, 100.0], logOption=1, ymin=1.e-2, ymax=1.e0, save='png')
	do.profile(9, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
	do.profile(12, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
	do.profile(15, 0.01, vLineCoords=[10.0, 100.0], logOption=1, save='png')
	do.stPlot(21, cmapType="viridis", logOption=2, save="png", xLabel=1)
	do.stPlot(12, cmapType="coolwarm", logOption=2, save="png", xLabel=1)
	do.stPlot(0, cmapType="viridis", logOption=1, save="png", xLabel=1)
	do.stPlot(9, cmapType="viridis", logOption=1, save="png", clim1=1.e8, clim2=1.e10, xLabel=1)
	plt.figure(figsize=(10, 4))
	plt.semilogy(do.t, do.lum)
	plt.ylabel("Integrated Luminosity (code units)")
	plt.xlabel("Time (code units)")
	plt.axhline(y=1.e0, color='k', linestyle='--')
	plt.savefig(do.pdfName, format='pdf', bbox_inches='tight')
	plt.savefig(do.savePath + "lightCurve.png", bbox_inches='tight')
	plt.clf()
	


	do.pdfName.close()









'''
if str(sys.argv[2])=="save": 
	t0 = time.time()
	do.makeMultiAnim(save="yes", timeCutFactor=int(sys.argv[3]), lengthInSeconds=int(60))
	print(time.time() - t0)
if str(sys.argv[2])=="show": do.makeMultiAnim(show="yes")
'''


























