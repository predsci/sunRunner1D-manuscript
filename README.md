# sunRunner1D
The 1D version of the sunRunner tool

##########################################################

PLUTO INSTALLATION

The sunRunner1D tool is based on the PLUTO code. 

Hence PLUTO must first be downloaded and installed from here:

http://plutocode.ph.unito.it/

If you are not familiar with the PLUTO code we recommend that you spend time reading the manual and running the examples. If you would like to understand how the boundary conditions (BC) can be defined please see chapter 5 of the PLUTO manual.  

The 2D plotting function of sunRunner1D uses the pyPLUTO Python library.  If you would like to have this capability please create a miniconda, or Conda, environment and install the pyPLUTO package. Details on how to install the package are provided in the Doc/pyPLUTO.html file (of your downloaded PLUTO package). 


COMPILING PLUTO FOR sunRunner1D

We have made modifications to only two of the PLUTO routines: init.c and userdef_output.c These routines along with the required definitions.h file are in the 'src' sub-directory.  To compile the PLUTO code for sunRunner1D please go to the src sub-directory:

% cd src

Invoke the PLUTO setup.py script (see Section 1.3 in the PLUTO manual)

%  python $PLUTO_DIR/setup.py

Once you are done with this configuration setup you will see a makefile and sysconf.out file that are specific to your computer platform (we recommend selecting the gcc compiler in setup.py) and you can now compile PLUTO:

%make 

Clean the directory:

%make clean

And go back to the directory with the sunRunner1D.py script

%cd ..

If you have installed the pyPLUTO library, and would like to use the 2D capabilities of sunRunner1D, please be sure  to activate the Conda environment of pyPLUTO.  Before running sunRunner1D.py please change the first line in this script to point to your python installation and if you have NOT installed pyPLUTO you will need to comment some lines in sunRunner1D.py and util.py


RUNNING sunRunner1D.py

To run the script use:
%./sunRunner1D.py 

You will be prompted to enter an event number [1-4] or your own event number.  We recommend that you start with events 1-4 of the manuscript.   Then you can use the provided pluto_5.ini input file and either run it as is (event number will be 5 in this case) or edit it and then run.  The logic of the sunRunner1D.py script is that it expects a pluto_X.ini file when running event number 'X'.

The output of a run is saved in an event specific sub-directory: runs/event_X/ (where X is your event number).

The following plots will be generated:

event_X_ts.png - A Time series plot at the observer point (defined by the OBS_LOC parameter in pluto_X.ini file) 

event_X_bc.png - A Time series plot of the perturbation at the inner boundary

event_X_width.ong - The width of the CME as a function of radial distance 

event_X_2d.png - 2D plots (distance/time) of the event 


PARAMETERS in pluto.ini

Upstream background values:

T_0   - upstream Temperature
RHO_0 - upstream density
V_0   - upstream velocity
BP_0  - upstream magnetic field

Perturbation parameters;

RHO_PERT - density perturbation
V_PERT   - velocity perturbation
BP_PERT  - magnetic field perturbation

CME_START      - start time of CME
CME_DURATION   - duration of CME

OBS_LOC        - location of observer (used by init.c to produce data for the time-series plot event_X_ts.png) 

You can modify any of these parameters but we recommend that you do not change the order of the parameters.  If you do, you must also edit the definitions.h file (in src sub-directory) and recompile the code.  


DATA

The observations for events 1-4 are in the data/eventX_obs.csv files [X=1,2,3,4]


FIGS-PDF

The figs-pdf directory includes pdf files of events 1-5 with the following name convention:

eventX_ts.pdf - Time series comparison between the model and observation (e.g, figure 1 in manuscript)

eventX_bc.pdf - Boundary conditions for the event (e.g., figure 2 in manuscript)

eventX_r.pdf - Radial speed, radial density, transverse magnetic field and tempearture as a function of heliospheric
distance for three snapshots in time (e.g., figure 3 in manuscript)

eventX_bc.pdf - The width of the ICME as function of heliospheric distance (e.g., figure 4 in manuscript)

eventX_2d.pdf - Radial speed, radial density, transverse magnetic field and tempearture as a function of time (x-axis) and
heliospheric distance (y-axis, e.g, figure 5 in manuscript)


HELP and SUPPORT

For questions and/or help please email us: pete@predsci.com, mbennun@predsci.com



