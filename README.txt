DESCRIPTION

This folder contains the  network simulator that simulates memory engram formation 
in a population consisting of excitatory and inhibitory neurons with independent
dendritic subunits. 

Details of the simulator model are available in: 

Kastellakis, G., Silva, A. J., & Poirazi, P. (2016).
Linking memories across time via neuronal and dendritic overlaps in model neurons
with active dendrites. Cell reports, 17(6), 1491-1504.


Directory layout:

data/ : Contains simulation output data. used to generate figures

src/  : Contains the implementation of the simulator. Specifically:

	src/lamodel.cpp:    Main simulator entry point with command line option parsing
	src/constructs.h:   Data structure definitions
	src/constructs.cpp: Implementation file of simulation dynamics / connectivity and plasticity
	src/tests.cpp: 	    Unit tests 

figs/ : Figure-generating python scripts (requires the simulator output data)
	figs/engrams.py: generates main text figure and the .txt files for the figure data
	figs/supl.py: generates supplemental figure

.:
	run_simulations.sh:  Script to run all simulations serially
	submit_lamodel.sh:   Submission script used to run the simulations in a PBS compatible cluster  (not used by default)


REQUIREMENTS 

gcc 4.4.7 
python 2.7
GNU Make 3.81


RUNNING

To compile the simulator and generate  data:

> make -C src clean all
> sh run_simulations.sh


To generate figures 

> cd figs
> python  engrams.py
> python  supl.py


To run unit tests:

> ./tests



LICENSE
GPLv2 (included)
