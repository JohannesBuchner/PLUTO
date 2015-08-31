NAME : pyPLUTO 

TASK : Quick Tool for Visualization of PLUTO 4 data [Mignone 2007]

AUTHOR : Bhargav Vaidya [University of Torino, Italy.]

CONTRIBUTION : A. Strugarek (Dept. of Physics, University of Montreal)
	       D. Stepanovs (MPI Astronomy, Heidelberg)

DESCRIPTION:

The code is completely written using the Python Language. 
Further the GUI is developed with the Tkinter Interface.

Features of this code: 

1. Completely based on Python and easy to work without need of any licenses
like in IDL. 
2. The GUI environment provides a tool for quick-check of data during the
simulations runs. 
3. The code is user friendly and allows the user to even do further plotting
of contours, velocity vectors
on the surface plot.
4. Also the code can read the saved user defined variables. 


CHANGES from 4-1.0 to 4-2.0:

1. The Data reader is made 10x faster using the array.fromstring module 
replacing the previously used struct module.

2. The reader can also now read VTK data files generated using PLUTO code.
Python/VTK wrapping library is NOT required for the same.

3. HDF5 Reader is also provided with this version which can analyse and visualize
PLUTO-Chombo AMR data files.    

MANUAL :

Further details of the Installation and Getting Started can be found in 
the html documentation - PLUTO/Doc/pyPLUTO.html 
