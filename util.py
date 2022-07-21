import sys
import os
import getopt
import shutil
import subprocess

from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd 


# needed for 2D plots
import pyPLUTO as pypl
import pyPLUTO.pload as pp
import pyPLUTO.Image as img
import pyPLUTO.Tools as tl
# end 

# Conversion factors

mp=1.6726e-24
b_fac_pluto=0.0458505
r_fac_pluto=2.149e+02
temp_fac_pluto=1
rho_fac_pluto=mp
time_fac_pluto=1.49597871e+08/3600 


def find_output_dir_name(pluto_file):

	""" read pluto.ini file and find the name of the output_dir """

	file1 = open(pluto_file, "r")

	string1 = 'output_dir'
	
	# setting flag  to 0
	flag = 0

	# Loop through the file line by line 

	for line in file1:

		# checking string is present 

		if string1 in line:

			flag = 1
			break


	if (flag == 0): 
		print(' ')
		print(' Error: Missing output_dir name in pluto input file ')
		print(' Code will Exit')
		print(' ')
		quit()

	index = line.index('/')

	wdir = line[(index+1):-1]

	return wdir


def find_nearest(array, value):

	""" Given an array and a value find index and value of nearest element """
	array = np.array(array)
	idx = (np.abs(array - value)).argmin()
	return idx, array[idx]


def read_cme_params(pluto_ini_file):

	""" Read pluto.ini file fo run and retrieve start time 
	of CME and duration. Return both in code units"""


	cme_params_dict = {'CME_START_TIME': 0.0, 'CME_DURATION': 0.0, 'OBS_LOC': 1.0, 'X1-grid': 1.0}

	with open(pluto_ini_file, 'r') as fp:
		line = fp.readline()
		while line:
			newline = line;
			for sub in cme_params_dict.keys():
				if (line.find(sub) != -1):
					line2 = line.strip().split(' ')
					last = len(line2) - 1
					cme_params_dict[sub] = line2[last]
			line = fp.readline()

	start = float(cme_params_dict['CME_START_TIME']) * time_fac_pluto
	dur   = float(cme_params_dict['CME_DURATION']) * time_fac_pluto
	r1    = float(cme_params_dict['OBS_LOC']) * r_fac_pluto
	rmax  = float(cme_params_dict['X1-grid']) * r_fac_pluto


	return start, dur, rmax, r1



def read_obs_dat(obs_file):

	""" read observer file and convert units to hours nT"""

	df = pd.read_table(obs_file,delim_whitespace=True, skiprows=1)
	df['time'] = df['time'] * time_fac_pluto
	df['Bp']   = df['Bp'] * b_fac_pluto

	return (df)

def read_tracer_dat(tracer_file):

	""" read tracers file """

	df = pd.read_table(tracer_file,delim_whitespace=True, skiprows=0)
	df['time'] = df['time'] * time_fac_pluto
	df['tracer1'] = df['tracer1'] * r_fac_pluto
	df['tracer2'] = df['tracer2'] * r_fac_pluto

	return (df)

def plot_vars_at_bc(df, cme_start_time, cme_duration, filename):

	""" Plot Perturbations at inner boundary """

	time = df['time']
	bpr0 = df['Bp']
	vrr0 = df['vr']
	npr0 = df['rho']
	tr0 = df['temp']

	cme_end_time = cme_start_time + cme_duration


	#########################################################################
	#
	# Find the time range we should use for the plot at the inner boundary
	#


	time_max = cme_start_time + cme_duration / 2.0 #time[imax]
	
	wwindow = round(cme_duration)
	tmin = round(time_max) - wwindow
	tmax = round(time_max) + wwindow

	#########################################################################
	# 
	# Plot Vr, Bp, and np perturbations at the inner boundary 
	#

	f1 = figure(figsize=[7,8], num=1)

	ax1 = f1.add_subplot(311)
	time = np.array(time)
	ydata = np.array(vrr0) 
	ax1.plot(time, ydata, color = 'red')
	ax1.axvline(x = cme_start_time, color = 'grey', linestyle = '--')
	ax1.axvline(x = cme_end_time, color = 'grey', linestyle = '--')
	xlim([tmin,tmax])

	ylabel(r'V [km/s]')
	title("Initial Perturbation at Inner Boundary")
	ax1.grid(which='major',axis='y')
	ax1.axes.xaxis.set_ticklabels([])

	ax1 = f1.add_subplot(312)
	ydata = np.array(bpr0) 
	ax1.plot(time, ydata, color = 'blue')
	ax1.axvline(x = cme_start_time, color = 'grey', linestyle = '--')
	ax1.axvline(x = cme_end_time, color = 'grey', linestyle = '--')
	xlim([tmin,tmax])
	#xlabel(r'Time [hours]')
	ylabel(r'B$_{\rm \phi}$ [nT]')

	ax1.grid(which='major',axis='y')
	ax1.axes.xaxis.set_ticklabels([])

	ax1 = f1.add_subplot(313)
	ydata = np.array(npr0) 
	ax1.plot(time, ydata, color = 'green')
	ax1.axvline(x = cme_start_time, color = 'grey', linestyle = '--')
	ax1.axvline(x = cme_end_time, color = 'grey', linestyle = '--')
	xlim([tmin,tmax])
	xlabel(r'Time [hours]')
	ylabel(r'n [$cm^{-3}$]')

	ax1.grid(which='major',axis='y')
	plt.savefig(filename)

	plt.clf()
	plt.cla()
	plt.close()
	return


def plot_vars_at_r1(df, cme_start_time, cme_duration, filename):

	""" Plot Time Series of variables at observer point """

	time  = df['time'] 
	bp1au = df['Bp']
	vr1au = df['vr']
	np1au = df['rho']
	t1au  = df['temp']

	cme_end_time = cme_start_time + cme_duration
	#########################################################################
	#
	# Find the time range we should use for the plot
	#

	imax = np.argmax(vr1au)
	
	time_max = time[imax] # This is the time that Vr has a maximum

	# # find index of CME start time 
	# itmin, time_min = find_nearest(time, cme_start_time)

	itmax = len(time)-1   


	# set tmin to 24 hours before the maximum in velocity of the time-series
	tmax = round(time[itmax])
	tmin = time_max - 24.0

	itmin, time_min = find_nearest(time, tmin)
	#########################################################################
	#
	# Plot the variables at 1AU
	#

	f1 = figure(figsize=[15,8], num=1)

	ax1 = f1.add_subplot(411)
	time_1au = np.array(time)
	ydata = np.array(vr1au) 
	ax1.plot(time_1au, ydata, color = 'red')
	ymin=min(ydata[itmin:itmax])
	ymax=max(ydata[itmin:itmax])
	xlim([tmin,tmax])
	ylim([0.9*ymin, 1.1*ymax])
	xlabel(r'Time [hours]')
	ylabel(r'v$_{r}$ [km/s]')

	ax1 = f1.add_subplot(412)
	ydata = np.array(np1au) 
	ax1.plot(time_1au, ydata, color = 'green')
	ymin=min(ydata[itmin:itmax])
	ymax=max(ydata[itmin:itmax])
	yscale('log')

	xlabel(r'Time [hours]')
	ylabel(r'n [$cm^{-3}$]')
	xlim([tmin,tmax])
	ylim([0.9*ymin, 1.1*ymax])

	ax1 = f1.add_subplot(413)
	ydata = np.array(bp1au) 
	ax1.plot(time_1au, ydata, color = 'blue')

	ymin=min(ydata[itmin:itmax])
	ymax=max(ydata[itmin:itmax])         
	xlabel(r'Time [hours]')
	ylabel(r'B$_{\rm \phi}$ [nT]')
	xlim([tmin,tmax])
	ylim([0.9*ymin, 1.1*ymax])

	ax1 = f1.add_subplot(414)
	ydata = np.array(t1au) 
	ax1.plot(time_1au, ydata, color = 'magenta')

	yscale('log')
	ymin=min(ydata[itmin:itmax])
	ymax=max(ydata[itmin:itmax])     
	xlabel(r'Time [hours]')
	ylabel(r'T [K]')
	xlim([tmin,tmax])
	ylim([0.9*ymin, 1.1*ymax])
	plt.savefig(filename)
	plt.clf()
	plt.cla()
	plt.close()

	#plt.show()

	return


def plot_cme_width(df, rmax, filename):

	""" Plot CME width as a function of distance"""

	time  = df['time'] 
	tracer1 = df['tracer1']
	tracer2 = df['tracer2']

	npts = len(time)

	# find when the CME ends - this is when tracer2 begins to move
	istart, ifirst = 0, 0
	for ii in range(1, npts):
		if (tracer2[ii] > tracer2[(ii-1)] and ifirst == 0):
			istart = ii
			ifirst = 1


	# find when the first tracer exits the domain
	iend, rval = find_nearest(tracer1, rmax)
	
	npts = iend - istart 

	pl_tr1, pl_tr2, pl_time = np.zeros(npts), np.zeros(npts), np.zeros(npts)

	icount =0
	for ii in range(istart, iend):
		pl_tr1[icount] = tracer1[ii]
		pl_tr2[icount] = tracer2[ii]
		pl_time[icount] = time[ii]
		icount = icount +1

	r_mid, w_cme = np.zeros(npts), np.zeros(npts)

	for ii in range(npts):

		r_mid[ii] = (pl_tr1[ii] + pl_tr2[ii])/2.0 
		w_cme[ii] = pl_tr1[ii] - pl_tr2[ii]
	
	
	f1 = figure(figsize=[5,5], num=1)

	fig, ax1 = plt.subplots(figsize = [15,5], num = 1)
	ax1.plot(r_mid, w_cme, 'g', linewidth =4)
	ax1.grid(True)
	xlabel(r'R [Rs]')
	ylabel(r'CME Width [Rs]')	
	ylim([min(w_cme)*0.95, max(w_cme)*1.05])

	plt.savefig(filename)

	plt.clf()
	plt.cla()
	plt.close()

	return

def read_pluto_output(pluto_dir):

	""" Read 2-D PLUTO output """

	nlinf = pypl.nlast_info(w_dir=pluto_dir)
	D=pp.pload(nlinf['nlast'],w_dir=pluto_dir,datatype='dbl')

	time2d = pypl.read_time(w_dir=pluto_dir,datatype='dbl')

	time2d = np.array(time2d, dtype=np.float32) 

	r=D.x1 * r_fac_pluto
	time2d = time2d * time_fac_pluto
	nr=D.n1
	nt = nlinf['nlast']
	nt1=nt+1
	vr2d = np.zeros((nr,nt1))
	np2d = np.zeros((nr,nt1))
	bp2d = np.zeros((nr,nt1))
	t2d = np.zeros((nr,nt1))

	for nframe in range(0, nt1):
		D=pp.pload(nframe,w_dir=pluto_dir,datatype='dbl')
		vr2d[:,nframe] = D.vx1
		np2d[:,nframe] = D.rho
		bp2d[:,nframe] = D.Bx3
		t2d[:, nframe] = D.T


	return(vr2d, np2d, bp2d, t2d, time2d, r)


def plot_pluto_2d(vr2d, np2d, bp2d, t2d, time2d, r, r1, cme_start_time, cme_duration, filename):

	""" Plot 2-D PLUTO output """

	shrink=0.8

	tmin=cme_start_time - 24

	time2d = time2d - tmin
	tmin = 0.0
	tmax=max(time2d)

	tmin2d = tmin 
	tmax2d = tmax 

	I=img.Image()
	mfig = figure(8,figsize=[15,10])
	Nrows,Ncols=2,2

	# Vr
	icount = 1
	data = vr2d
	mfig.add_subplot(Nrows,Ncols,icount)
	plt.pcolormesh(time2d, r, data,vmin=None,vmax=None,shading='auto',cmap='rainbow')
	plt.ylabel('R [R$_s$]')
	plt.axis([tmin2d,tmax2d,r.min(),r.max()])
	colorbar(orientation='vertical',shrink=shrink)
	mytitle = 'Vr [km/s] '
	title(mytitle)
	plt.axhline(y = r1, color = 'black', linestyle = '--')


	# Np
	icount = icount + 1
	data = np2d
	mfig.add_subplot(Nrows,Ncols,icount)
	plt.pcolormesh(time2d, r, data,norm=LogNorm(vmin=data.min(), vmax=data.max()), shading='auto')
	plt.ylabel('R [R$_s$]')
	plt.axis([tmin2d,tmax2d,r.min(),r.max()])
	colorbar(orientation='vertical',shrink=shrink)
	mytitle = 'n [$cm^{-3}$]'
	title(mytitle)
	plt.axhline(y = r1, color = 'black', linestyle = '--')

	# Bp
	icount = icount + 1
	data = bp2d
	mfig.add_subplot(Nrows,Ncols,icount)
	plt.pcolormesh(time2d, r, data* b_fac_pluto, vmin = -50, vmax = 50, shading='auto',cmap='RdBu_r')
	plt.ylabel('R [R$_s$]')
	plt.xlabel('$Time [hours]$')
	plt.axis([tmin2d,tmax2d,r.min(),r.max()])
	colorbar(orientation='vertical',shrink=shrink)
	mytitle = 'B$\\phi$'
	title(mytitle)
	plt.axhline(y = r1, color = 'black', linestyle = '--')


	# Temperature
	icount = icount + 1
	data = t2d
	mfig.add_subplot(Nrows,Ncols,icount)
	plt.pcolormesh(time2d, r, data, norm=LogNorm(vmin=data.min(), vmax=data.max()), shading='auto',cmap='rainbow')
	plt.ylabel('R [R$_s$]')
	plt.xlabel('$Time [hours]$')
	plt.axis([tmin2d,tmax2d,r.min(),r.max()])
	colorbar(orientation='vertical',shrink=shrink)
	mytitle = 'T [K]'
	title(mytitle)
	plt.axhline(y = r1, color = 'black', linestyle = '--')


	plt.savefig(filename)

	print('Saving Figure ', filename)
	plt.clf()
	plt.cla()
	plt.close()


	return
