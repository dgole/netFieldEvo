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
from matplotlib.backends.backend_pdf import PdfPages
import resource
import time
import numpy.random
import matplotlib.animation as animation
from matplotlib.pylab import *
from mpl_toolkits.axes_grid1 import host_subplot
m.rcParams['text.usetex'] = True
m.rcParams['text.latex.unicode'] = True

font = {'family' : 'normal', 'weight' : 'bold', 'size'   : 18}
m.rc('font', **font)

class Data:
	def getrindex(self, r1):
		return (np.abs(self.r-r1)).argmin()

	def gettindex(self, t1):
		return (np.abs(self.t-t1)).argmin()

	def __init__(self, path, savePath=None):
		print "initializing data structure from " + str(path)
		self.path= path
		sgrid = np.load(self.path+"sgrid.npy")
		dgrid = np.load(self.path+"dgrid.npy")
		state = np.load(self.path+"state.npy")
		time  = np.load(self.path+"time.npy")
		self.data = []
		for col in range(dgrid.shape[0]): self.data.append(dgrid[col])
		for col in range(state.shape[0]): self.data.append(state[col])
		if savePath is None:
			self.savePath = self.path
		if not os.path.exists(self.savePath): os.makedirs(self.savePath)
		ntRightNow = self.data[0].shape[0]	
		print "I want the time to be " + str(len(time))
		print "but I have to make it " + str(ntRightNow)
		ntRightNow = ntRightNow - 10
		print "but to be even safer I'm making it " + str(ntRightNow)
		self.r = sgrid[0]
		self.dr = sgrid[1]
		self.Omega = sgrid[2]
		self.tOrbit = (2.0*3.14159)/self.Omega
		self.vKep=self.Omega*self.r
		self.t = time[:ntRightNow]
		self.dt = self.t-np.roll(self.t,1); self.dt[0]=self.dt[1]
		self.dt = self.dt[:ntRightNow]
		self.header = [r"$\alpha$", r"$h/r$", r"$T_c^4$", r"$T_disk^4$", r"$c_s$", r"$\rho$", r"$\kappa_R$", r"$\nu$", r"$\tau$", r"$\dot{M}$", r"$v_{adv}$", r"$v_{diff}$", r"$B_z$", r"$B_{rs}$", r"$\psi$", r"$\Sigma$", r"$\beta$"]  		
		self.pdfName = PdfPages(self.savePath + "/plots.pdf")
		self.tmax = self.t.max()
		self.rmax = self.r.max()
		self.rmin = self.r.min()
		self.nr = self.r.shape[0]
		self.nt = self.t.shape[0]
		self.data.append(self.data[13]/self.data[12]); self.header.append(r"$B_{rs}/B_z$")	
		self.data.append(self.data[10]/self.data[11]); self.header.append(r"$Pr_{eff}$")
		self.data.append(-self.r/self.data[10]); self.header.append(r"$t_{adv}$")		
		self.data.append(self.r/self.data[11]); self.header.append(r"$v_{adv}$")
		self.data.append((3.0/4.0)*self.data[9]*np.square(self.Omega)*self.r*self.dr); self.header.append("dFlux")
		self.lum = np.sum(self.data[21], axis=1)/np.sum(self.data[21][self.gettindex(10)]); 		
			
	def getrindex(self, r1):
		return (np.abs(self.r-r1)).argmin()

	def gettindex(self, t1):
		return (np.abs(self.t-t1)).argmin()

	def getLogSpacedTimeList(self, nProfiles):
		dtLog = np.log10(self.tmax+1.0)/(nProfiles-1.0)
		logtArray = np.asarray([n*dtLog for n in range(nProfiles)])
		tArray = np.power(10,logtArray)-1.0 
		return tArray

	def getLinSpacedTimeList(self, nProfiles):
		dt = self.tmax/(nProfiles)
		tArray = np.asarray([n*dt for n in range(nProfiles+1)])
		return tArray

	def profile(self, col, time, vLineCoords=[1.0, 50.0], hLineCoords=[], logOption=0, save=None, savePath=None, legendLabel=None, ymin=None, ymax=None):
		n=self.gettindex(time)
		if savePath is None:
			savePath=self.savePath
		plotData = self.data[col][n,...]; title = self.header[col];
		if logOption==0 and save is not None:
			plt.semilogx(self.r, plotData);	plt.ylabel(self.header[col]);
		if logOption==1 and save is not None:
			plt.loglog(self.r, np.absolute(plotData));	plt.ylabel(self.header[col]);
		if logOption==0 and save is None:
			plt.semilogx(self.r, plotData, label=legendLabel);	plt.ylabel(self.header[col]);
		if logOption==1 and save is None:
			plt.loglog(self.r, np.absolute(plotData), label=legendLabel);	plt.ylabel(self.header[col]);
		if ymin is None and ymax is None:
			plt.ylim(np.abs(plotData).min()*0.8, np.abs(plotData).max()*1.1)
		if ymin is not None and ymax is None:
			plt.ylim(ymin, plotData.max()*1.1)
		if ymin is None and ymax is not None:
			plt.ylim(plotData.min()*0.8, ymax)
		if ymin is not None and ymax is not None:
			plt.ylim(ymin, ymax)
		plt.xlabel("r (code units)");
		for xc in vLineCoords: plt.axvline(x=xc, color='k', linestyle='--')
		for yc in hLineCoords: plt.axhline(y=yc, color='k', linestyle='--')
		plt.tight_layout()
		if save=="png":
			plt.savefig(savePath + "profile_" + str(col) + ".png", bbox_inches='tight')
			print "saved profile plot for column " + str(col) + " to png"
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
			print "saved profile plot for column " + str(col) + " to pdf"
		if save is None:
			plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
			#plt.legend()
		if save is not None:
			plt.clf()

	def multiProfile(self, col, nProfiles, spacing="log", vLineCoords1=[1.0,50.0], hLineCoords1=[], logOption=0, save=None, savePath=None, ymin=None, ymax=None):
		if spacing=="lin":
			tList = self.getLinSpacedTimeList(nProfiles)
		if spacing=="log":
			tList = self.getLogSpacedTimeList(nProfiles)
		for n in range(len(tList)): 
			self.profile(col, tList[n], logOption=logOption, legendLabel="t = " + str(np.round(tList[n],1)), ymin=ymin, ymax=ymax, vLineCoords=vLineCoords1, hLineCoords=hLineCoords1)
		if save=="png":
			plt.savefig(savePath + "multiProfile_" + str(col) + ".png", bbox_inches='tight'); plt.clf()
			print "saved multi profile plot for column " + str(col) + " to png"
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight'); plt.clf()
			print "saved multi profile plot for column " + str(col) + " to pdf"

	def makeMultiAnim(self, timeCutFactor=10, lengthInSeconds=20, savePath=None, show=None, save=None):
		# Setup figure and subplots
		f0 = figure(num = 0, figsize = (16, 12))
	
		colList = [12, 0, 16, 15, 1, 9]		
		
		ax = []
		ax.append(subplot2grid((4, 2), (0, 0)))
		ax.append(subplot2grid((4, 2), (0, 1)))
		ax.append(subplot2grid((4, 2), (1, 0)))
		ax.append(subplot2grid((4, 2), (1, 1)))
		ax.append(subplot2grid((4, 2), (2, 0)))
		ax.append(subplot2grid((4, 2), (2, 1)))
		ax.append(subplot2grid((4, 2), (3, 0), colspan=2))
		
		i=0; ax[i].set_ylim(1.1*np.min(self.data[colList[i]]), 1.1*np.max(self.data[colList[i]]))
		i=1; ax[i].set_ylim(0.7*np.min(np.abs(self.data[colList[i]])), 1.5*np.max(np.abs(self.data[colList[i]])))
		i=2; ax[i].set_ylim(0.7*np.min(np.abs(self.data[colList[i]])), min(1.5*np.max(np.abs(self.data[colList[i]])),1.e18))
		i=3; ax[i].set_ylim(0.7*np.min(np.abs(self.data[colList[i]])), 1.5*np.max(np.abs(self.data[colList[i]])))
		i=4; ax[i].set_ylim(0.7*np.min(np.abs(self.data[colList[i]])), 1.5*np.max(np.abs(self.data[colList[i]])))
		i=5; ax[i].set_ylim(0.7*np.min(np.abs(self.data[colList[i]])), 1.5*np.max(np.abs(self.data[colList[i]])))
		i=6; ax[i].set_ylim(0.7*np.min(self.lum), 1.5*np.max(self.lum))
		
		for axObject in ax[0:6] : axObject.set_xlim(self.rmin, self.rmax)
		ax[6].set_xlim(0.0, self.tmax)

		i=0; ax[i].set_ylabel(self.header[colList[i]])
		i=1; ax[i].set_ylabel(self.header[colList[i]])
		i=2; ax[i].set_ylabel(self.header[colList[i]])
		i=3; ax[i].set_ylabel(self.header[colList[i]])
		i=4; ax[i].set_ylabel(self.header[colList[i]])
		i=5; ax[i].set_ylabel(self.header[colList[i]])
		i=6; ax[i].set_ylabel(r"$L_{BB}/L_{BB,0}$")
		
		ax[6].set_xlabel("t")

		# Data Placeholders
		yp1=zeros(0)
		t=zeros(0)

		# set plots
		p0, = ax[0].semilogx(t,yp1)
		p1, = ax[1].loglog(t,yp1)
		p2, = ax[2].loglog(t,yp1)
		p3, = ax[3].loglog(t,yp1)
		p4, = ax[4].loglog(t,yp1)
		p5, = ax[5].loglog(t,yp1)
		p6, = ax[6].semilogy(t,yp1)

		def updateData(n):
			n1 = n*timeCutFactor
			i=0; p0.set_data(self.r, self.data[colList[i]][n1])
			i=1; p1.set_data(self.r, np.abs(self.data[colList[i]][n1]))
			i=2; p2.set_data(self.r, np.abs(self.data[colList[i]][n1]))
			i=3; p3.set_data(self.r, np.abs(self.data[colList[i]][n1]))
			i=4; p4.set_data(self.r, np.abs(self.data[colList[i]][n1]))
			i=5; p5.set_data(self.r, np.abs(self.data[colList[i]][n1]))
			p6.set_data(self.t[0:n1], np.abs(self.lum[0:n1]))
			f0.suptitle("t = " + str(np.round(self.t[n1],1)), fontsize=16)
			return p0, p1, p2, p3, p4, p5, p6

		# interval: draw new frame every 'interval' ms
		# frames: number of frames to draw
		nFrames = self.nt/timeCutFactor
		framesPerSecond = nFrames/lengthInSeconds
		simulation = animation.FuncAnimation(f0, updateData, blit=False, frames=nFrames, interval=10, repeat=False)
		# Uncomment the next line if you want to save the animation
		#simulation.save(filename='sim.mp4',fps=30,dpi=600)
		if show is not None:
			plt.tight_layout()
			plt.show()
		if save is not None:
			Writer = animation.writers['ffmpeg']
			writer = Writer(fps=framesPerSecond, metadata=dict(artist='Me'), bitrate=1800)
			simulation.save(self.savePath+'anim1.mp4', writer=writer)

	def timeEvo(self, col, r, logOption=0, save=None, savePath=None, legendLabel=None, ymin=None, ymax=None):
		i=self.getrindex(r) 
		if savePath is None:
			savePath=self.savePath
		plotData = self.data[col][...,i]; title = self.header[col];
		if logOption==0 and save is not None:
			plt.plot(self.t, plotData);	plt.ylabel(self.header[col]);
		if logOption==1 and save is not None:
			plt.semilogy(self.t, np.absolute(plotData));	plt.ylabel(self.header[col]);
		if logOption==0 and save is None:
			plt.plot(self.t, plotData, label=legendLabel);	plt.ylabel(self.header[col]);
		if logOption==1 and save is None:
			plt.semilogy(self.t, np.absolute(plotData), label=legendLabel);	plt.ylabel(self.header[col]);
		if ymin is not None and ymax is None:
			plt.ylim(ymin, plotData.max()*1.1)
		if ymin is None and ymax is not None:
			plt.ylim(plotData.min()*0.8, ymax)
		if ymin is not None and ymax is not None:
			plt.ylim(ymin, ymax)
		plt.xlabel("t (code units)");
		plt.tight_layout()
		if save=="png":
			plt.savefig(savePath + "timeEvo_" + str(col) + ".png", bbox_inches='tight')
			print "saved time evo plot for column " + str(col) + " to png"
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
			print "saved time evo plot for column " + str(col) + " to pdf"
		if save is None:
			plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
			#plt.legend()
		if save is not None:
			plt.clf()




'''
	



	def multiTimeEvo(self, col, rList, logOption=0, save=None, savePath=None, ymin=None, ymax=None):
		for r in rList: 
			self.timeEvo(col, r, logOption=logOption, legendLabel="r = " + str(np.round(r,1)), ymin=ymin, ymax=ymax)
		if save=="png":
			plt.savefig(savePath + "multiTimeEvo_" + str(col) + ".png", bbox_inches='tight')
			print "saved multi time evo plot for column " + str(col) + " to png"
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
			print "saved multi time evo plot for column " + str(col) + " to pdf"
		plt.clf()

	def makeAnim(self, animCol, logOption=0, savePath=None, show=None, save=None):
		if savePath is None:
			savePath=self.savePath
		fig = plt.figure()
		if logOption==1:
			ylim1=np.amin(np.abs(self.data[animCol]))*0.8
			ylim2=np.amax(np.abs(self.data[animCol]))*1.1
		if logOption==0:
			ylim1=np.amin(self.data[animCol])
			ylim2=np.amax(self.data[animCol])
			if ylim1<0:
				ylim1=ylim1*1.1
			else:
				ylim1=ylim1*0.8
			if ylim2<0:
				ylim2=ylim2*1.1
			else:
				ylim2=ylim2*0.8
		ax = plt.axes(xlim=(self.r[0], self.r[self.r.shape[0]-1]), ylim=(ylim1, ylim2))
		if logOption==1:
			line, = ax.loglog([], [], lw=2)
		if logOption==0:
			line, = ax.semilogx([], [], lw=2)
		def init():
			line.set_data([], [])
			return line, 
		if logOption==1:
			def animate(i):
				x = self.r
				y = np.abs(self.data[animCol][i,...])
				line.set_data(x, y)
				return line,
		if logOption==0:
			def animate(i):
				x = self.r
				y = self.data[animCol][i,...]
				line.set_data(x, y)
				return line,
		anim = animation.FuncAnimation(fig, animate, init_func=init, frames=self.nt, interval=20, blit=True)
		plt.xlabel("r (code units)")
		plt.ylabel(self.header[animCol])
		if show is not None:
			plt.tight_layout()
			plt.show()
		if save is not None:
			Writer = animation.writers['ffmpeg']
			writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
			anim.save(savePath+'anim.mp4', writer=writer)

	def stPlot(self, col, cmapType="viridis", logOption=0, save=None, savePath=None, clim1=None, clim2=None):
		print self.path + ": making ST plot for column " + str(col)
		if savePath is None:
			savePath=self.path
		if logOption==0:
			plotData = self.data[col]; title = self.header[col];
		if logOption==1:
			plotData = np.log10(np.absolute(self.data[col])); title = "log " + self.header[col];
		plt.imshow(np.transpose(np.fliplr(plotData)), extent=[0,self.tmax,self.rmin,self.rmax], aspect=(0.2*self.tmax/(self.rmax-self.rmin)), cmap=plt.get_cmap(cmapType))
		plt.title(title); plt.xlabel("Time (code units)"); plt.ylabel("r (code units)");
		plt.colorbar(shrink=0.5)
		if (clim1 is not None and clim2 is not None):
			plt.clim(clim1,clim2)
		if (clim1 is not None and clim2 is None):
			plt.clim(clim1,np.amax(plotData))
		if (clim1 is None and clim2 is not None):
			plt.clim(np.aminx(plotData),clim2)
		plt.tight_layout()
		if save=="png":
			plt.savefig(savePath + "ST_" + str(col) + ".png", bbox_inches='tight')
			print "saved ST plot for column " + str(col) + " to png"
		if save=="pdf":
			plt.savefig(self.pdfName, format='pdf', bbox_inches='tight');
			print "saved ST plot for column " + str(col) + " to pdf"
		plt.clf()

	def addCol(self, funcName, label, *args, **kwargs):
		print self.path + ": adding column named " + label
		self.data.append(funcName(self.data))
		self.header.append(label)


		
'''











