#!/Users/michal/miniconda3/envs/psi/bin/python

import sys
import os
import getopt
import shutil
import subprocess
from util import *


#
# Write the welcome banner 
#
print('\n################################################')
print('################################################\n')
print(' ')
print(' Welcome to sunRunner1D 1.0 (07/12/22) ')
print(' This is a 1D MHD calculation with a simple CME ')
print(' You will need to provide a run number and have a pluto_X.ini file ')
print(' for your run with X being your run number ')
print(' To reproduce events 1-4 of the manuscript please use 1,2,3 or 4 ')
print(' For your own event please use the pluto_5.ini file as a starting ')
print(' point and edit the background and CME parameters,and possibly the grid range etc.')
print(' The output of the run can be found in the runs/event_X/output_dir ')
print(' directory. output_dir is defined in the pluto_X.ini file (typically output) ')
print(' The script determines this directory name and creates it ')
print(' Note that making the 2D plots requires the pyPLUTO package ')
print(' For more on that please see the PLUTO manual ')
print(' If you can not or do not want to install pyPLUTO comment the 2-D ')
print(' section of this script (line 170 onward) and the import commands in util.py ')
print('\n################################################')
print('################################################\n\n')



run_name = input("Please enter event number [1-4] or your own event name ")

print("\nYou entered: " + run_name)


##########################################################################
#
# Check that pluto_name.ini exits and retrieve name of sub-directory 
# for output of run 


template_file='pluto_'+str(run_name)+'.ini'

file_exists = os.path.exists(template_file)

if(file_exists == False):
	print(' ')
	print(' Error: Missing '+template_file+ ' File ')
	print(' Code Will Exit. ')
	print(' ')
	#quit()


wdir = find_output_dir_name(template_file)

run_dir = 'event_'+run_name

##########################################################################
#
# set up the working directory for the PLUTO run
#


if not os.path.exists('runs'):
	os.mkdir('runs')
if not os.path.exists('runs/'+run_dir):
	os.mkdir('runs/'+run_dir)
if not os.path.exists('runs/'+run_dir+'/'+wdir):
	os.mkdir('runs/'+run_dir+'/'+wdir)


print('\nAll output will be found in: runs' + '/' + run_dir)

##########################################################################
#
# copy the pluto run specific input file and the pluto executable to wdir
#


########################################################################
#
# copy the pluto_X.ini and and pluto executable 
# to PLUTO run directory 


src_path = os.getcwd()+'/src/'
dst_path = os.getcwd() + '/runs/'+run_dir

pluto_ini_org = os.getcwd() + '/' + template_file
pluto_exe_org = src_path+'pluto'

pluto_ini_dst = dst_path+'/pluto.ini'
pluto_exe_dst = dst_path+'/pluto'	


shutil.copyfile(    pluto_ini_org, pluto_ini_dst)
shutil.copyfile(pluto_exe_org, pluto_exe_dst)

########################################################################
#
# Run Pluto
#
  
os.chdir('runs/')
os.system('chmod -R 777 *')
os.chdir(run_dir)

print('\n\nRunning Pluto...For Progress see:\n\n')
print(os.getcwd()+'/out.txt')


with open('out.txt','w+') as fout:
	with open('err.txt','w+') as ferr:
		out=subprocess.call(["./pluto"],stdout=fout,stderr=ferr)
   
print('\n\n Run Completed Successfully \n\n')
print(' Output Files saved at: \n')


#########################################################################
#
# make  1-D plots 
#

cme_start_time, cme_duration, rmax, r1 = read_cme_params(pluto_ini_file = pluto_ini_dst)


# plots at inner boundary

df_bc = read_obs_dat(obs_file = dst_path+'/'+wdir+'/obs_r0.dat')

plot_file = dst_path+'/event_'+run_name+'_bc.png'

print('\n Saving Plot of Variables at Inner boundary: ', plot_file)
	
plot_vars_at_bc(df_bc, cme_start_time, cme_duration, plot_file)

	
# Time series plots at R1

df_r1 = read_obs_dat(obs_file = dst_path+'/'+wdir+'/obs_r1.dat')


plot_file = dst_path+'/event_'+run_name+'_ts.png'

print('\nSaving Time-Series Plot of Variables at Observer Point: ', plot_file)

plot_vars_at_r1(df_r1, cme_start_time, cme_duration, plot_file)


# CME width plot - requires reading the tracer file

tracer_file = dst_path+'/'+wdir+'/tracers.dat'

df_tracers = read_tracer_dat(tracer_file)

plot_filename = dst_path+'/event_'+run_name+'_width.png'

print('Saving CME width Plot: ', plot_filename, '\n')

plot_cme_width(df_tracers, rmax, plot_filename)


#########################################################################
#
# make  2-D plots - requires pyPLUTO 
#

pluto_dir = dst_path + '/' + wdir + '/'

vr2d, np2d, bp2d, t2d, time2d, r = read_pluto_output(pluto_dir)
  
plot_2d_filename = dst_path+'/event_'+run_name+'_2d.png'

err = plot_pluto_2d(vr2d, np2d, bp2d, t2d, time2d, r, r1, cme_start_time, cme_duration,plot_2d_filename)


#########################################################################



