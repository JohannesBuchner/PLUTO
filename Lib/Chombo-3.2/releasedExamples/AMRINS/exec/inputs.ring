main.max_step = 10000  #max number of timesteps to compute
main.max_step = 10    #max number of timesteps to compute
main.max_time = 200.0  #stop time
#main.max_time = 16.0  #stop time
#main.num_cells = 16 16  #base level domain
#main.num_cells = 8 8 16  #base level domain
#main.num_cells = 32 32 64  #base level domain
#main.num_cells = 128 64 8  #base level domain
#main.num_cells =  32 32 96  #base level domain
main.num_cells =  32 32 48  #base level domain
#main.num_cells =  64 64 96  #base level domain
#main.max_level = 4
#main.max_level = 3
main.max_level = 2
#main.max_level = 1
#main.max_level = 0
main.ref_ratio = 2 2 2 2 2 2
#main.ref_ratio = 4 4 4 4 4 4
main.regrid_interval = 100000 10000
main.regrid_interval = 4 4 4 4 4 4
main.regrid_interval = 2 2 2 2 2 2 
main.block_factor = 8
#main.max_grid_size = 1024
#main.max_grid_size = 28
main.max_grid_size = 16
#main.max_base_grid_size = 1024 
main.max_base_grid_size = 32
main.fill_ratio = 0.8
#main.fill_ratio = 0.5
main.grid_buffer_size = 1
main.checkpoint_interval = 100
#main.checkpoint_interval = -1
main.plot_interval = 2
#main.plot_interval = 50
main.plotPrefix = plt.
main.cfl = 0.5
main.limitSolverCoarsening = 0 # solver parameter
main.verbosity = 1  #higher number means more verbose
#main.gridfile = grids.dat.128

main.is_periodic = 0 0  1

#ns.domainLength = 1.0 1.0 0.25
ns.init_shrink = 1.0
ns.tag_vorticity = 1
ns.vorticity_tagging_factor = 0.00250
ns.vorticity_tagging_factor = 0.00125
ns.project_initial_vel = 1
ns.init_pressures = 1
ns.num_init_passes = 1
ns.tags_grow = 1

#initial grids
ns.specifyInitialGrids = 0
#ns.initialGridFile = grids.init
ns.initVelFromVorticity = 1
#ns.backgroundVelocity = 1.0
ns.backgroundVelocity = 0.0

ns.viscosity = 0.000001
#ns.viscosity = 0.00

ns.num_scalars = 1
#ns.scal_diffusion_coeffs = 0.00 0.00
ns.scal_diffusion_coeffs = 0.01 0.02

#multigrid solver parameters
#ns.viscous_num_smooth_up = 3
#ns.viscous_num_smooth_down = 3

#ns.max_dt = 1.0e8
ns.reflux_momentum = 1

#if this is 0, freestream preservation correction is not applied
ns.applyFreestreamCorrection = 1

# extra plotfile variables

ns.writeDivergence = 1
ns.writeLambda = 1
ns.writeTimeDerivatives = 0
ns.writeVorticity = 1
ns.writeScalars = 1
ns.writeDscalDt = 0
ns.writeStrains = 0
ns.writeGradELambda = 0
ns.writeProcIDs = 0

projection.doSyncProjection = 1
projection.applyFreestreamCorrection = 1
projection.eta = 0.9

# multigrid solver parameters
#projection.numSmoothUp = 3
#projection.numSmoothDown = 3


# this is physical BC info
# 0 = solidWall, 1=inflow, 2=outflow, 3=symmetry, 4=noShear
physBC.lo = 4 4 4
physBC.hi = 4 4 4	
physBC.maxInflowVel = 1.0
