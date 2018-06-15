# -*- Mode: Makefile -*- 

# trace the chain of included makefiles
makefiles += releasedExamples_EBMenagerie_exec

# This is an embedded boundary example
USE_EB=TRUE

## Define the variables needed by Make.example

# the base name(s) of the application(s) in this directory
ebase = plane sphere cylinder ellipsoid torus threeSpheres cappedCylinder \
        polynomial ribbonTorus mobiusTorus borromeanRings nozzle          \
        twoSpheres channelWithSpheres datafile digitalElevationModel      \
        sphereFlake gaussian mollifiedBox helix helicoil spiral           \
        fuelElement swirl artery channel packedChannel turbine

# the location of the Chombo "lib" directory
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRTools Workshop EBTools AMRTools BoxTools

# the locations of the source code directories
base_dir = .
src_dirs = ../src

# input file for 'run' target
## NOTE: this relies on the 'foreach' loop used in the 'run-only' target in "Make.rules"
INPUT = 

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets

run-only: datafile$(DIM)d.ASCII digitalElevationModel.sfbay.dat


datafile: datafile$(DIM)d.ASCII

datafile$(DIM)d.ASCII: datafile$(DIM)d.ASCII.gz
	gzip -dc $@.gz > $@

digitalElevationModel: digitalElevationModel.sfbay.dat

digitalElevationModel.sfbay.dat: digitalElevationModel.sfbay.dat.gz
	gzip -dc $@.gz > $@
