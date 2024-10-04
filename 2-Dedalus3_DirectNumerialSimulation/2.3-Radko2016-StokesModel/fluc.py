'''
Dedalus script to simulate thermohaline-shear equations in 2D periodic domain.
um, Tm, Sm: <u>, <T>, <S> (<.> is horizontally averaged variable)
up, Tp, Sp: u', T', S'

Parameters: 
Pr (Prandtl), 
tau (diffusivity ratio), 
Rp (R_rho, diffusive density ratio), 
Pe (Peclet)

To run, restart, and plot using e.g. 4 processes:
    $ mpiexec -n 24 python3 fluc.py
    $ mpiexec -n 24 python3 fluc.py --restart
'''
import sys
import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
import math
import h5py
import logging
logger = logging.getLogger(__name__)

# Allow restarting via command line
restart = (len(sys.argv) > 1 and sys.argv[1] == '--restart')

dealias = 3/2 
pi = np.pi
Rp, Ri, Pe, tau = 2., 10., 1e2, 0.01 # figure 4
Pr = 10.  # Prandtl number
Lx, Lz = 64., 1.
Nx, Nz = 384, 192
stop_sim_time = 600 + 300*restart # Stopping criteria
# Bases
coords = d3.CartesianCoordinates('x','z')
dist = d3.Distributor(coords, dtype=np.float64)
# define the coordinate system
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
zbasis = d3.RealFourier(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)

# define fields
# um, Tm, Sm: <u>, <T>, <S> (<.> is horizontally averaged variable)
# up, Tp, Sp:  u',  T',  S'
p = dist.Field(name='p', bases=(xbasis,zbasis)) # pressure 
um = dist.VectorField(coords,name='um', bases=(zbasis)) # horizontally averaged u(z)
up = dist.VectorField(coords,name='up', bases=(xbasis,zbasis)) # fluctuation u'(x,z)
Sm = dist.Field(name='Sm', bases=(zbasis)) # horizontally averaged S(z)
Sp = dist.Field(name='Sp', bases=(xbasis,zbasis)) # fluctuation S'(x,z)
Tm = dist.Field(name='Tm', bases=(zbasis)) # horizontally averaged T(z)
Tp = dist.Field(name='Tp', bases=(xbasis,zbasis)) # fluctuation T'(x,z)
Ubg = dist.Field(bases=(zbasis)) #background velocity of basic state

# Substitutions
x, z = dist.local_grids(xbasis, zbasis) # get coordinate arrays in horizontal and vertical directions
ex, ez = coords.unit_vector_fields(dist) # get unit vectors in horizontal and vertical directions
# define velocity components
U = um + up # full velocity
T = Tm + Tp # full temperature
S = Sm + Sp # full salinity
umx = um@ex
umz = um@ez
upx = up@ex
upz = up@ez

# create constant sub-field for incompressible flow condition's equation
tau_p = dist.Field(name='tau_p') 
tau_um = dist.Field(name='tau_um') 
tau_up = dist.Field(name='tau_up') 
# because this term is only a contant added to the equation, we don't need to instantiate it for bases system

# First-order form: "div(A)" becomes "trace(grad_A)"
lap = lambda A: d3.div(d3.grad(A)) # First-order form: "lap(f)" becomes "div(grad_f)"
dx = lambda A: d3.Differentiate(A, coords['x']) 
dz = lambda A: d3.Differentiate(A, coords['z'])
vol_avg = lambda A: d3.Integrate(A)/(Nx*Nz)
# vol_avg = lambda A: d3.Average(A) # Horizontal mean of A 
# h_mean = lambda A: d3.Average(A,'x') # Horizontal mean of A <------- this function will get ERROR for computing 'grad(h_mean(p))'
h_mean = lambda A: d3.Integrate(A,'x')/Nx # <----- highly recommended for computing 'grad(h_mean(p))' instead
Ubg['g'] = np.sin(2*pi*z)

# Problem
problem = d3.IVP([p,tau_p,tau_um,tau_up,
                  um,up,
                  Tm,Tp, 
                  Sm,Sp], namespace=locals())

problem.add_equation("trace(grad(U))+tau_p = 0")
problem.add_equation("integ(p) = 0") # Pressure gauge

problem.add_equation("grad(h_mean(p)) - (Pr/Pe)*lap(um) - (4*pi*pi*Ri/(Rp-1))*(Tm-Sm)*ez + tau_um= 0")
problem.add_equation("grad(p-h_mean(p)) - (Pr/Pe)*lap(up) - (4*pi*pi*Ri/(Rp-1))*(Tp-Sp)*ez + tau_up= 0")

problem.add_equation("integ(U) = 0") # Velocity gauge

problem.add_equation("dt(Tm)              - (1./Pe)*dz(dz(Tm)) - umz = - dz(h_mean(upz*Tp))")
problem.add_equation("dt(Tp) + Ubg*dx(Tp) - (1./Pe)*lap(Tp) - upz = - up@grad(Tm)-um@grad(Tp)-up@grad(Tp)+dz(h_mean(upz*Tp))")

problem.add_equation("dt(Sm)             - (tau/Pe)*dz(dz(Sm)) - Rp*(umz) = - dz(h_mean((upz)*Sp))")
problem.add_equation("dt(Sp) + Ubg*dx(Sp) - (tau/Pe)*lap(Sp) - Rp*(upz) = - up@grad(Sm)-um@grad(Sp)-up@grad(Sp)+dz(h_mean((upz)*Sp))")


# timestepper = d3.RK443
timestepper = d3.RK222

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Initial conditions
if not restart:
    p.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    up.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    Tp.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    Sp.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    file_handler_mode = 'overwrite'
    initial_timestep = 0.02
    max_timestep = 0.02
else:
    write, initial_timestep = solver.load_state('/home/jms24002/reproduce_radko_results/fig5/checkpoints/checkpoints_s90.h5')
    max_timestep = 0.02
    file_handler_mode = 'append'
    logger.info('Imported last-step data successfully')

# store instantaneous data for analysis later
# first, we save horizontal averaged fields (1D array) following small output time-step
timehistory = solver.evaluator.add_file_handler('timehistory', sim_dt=1.0, max_writes=10000, mode=file_handler_mode)
timehistory.add_task(umx, name='um') # horizontal averaged x-axis velocity
timehistory.add_task(umz, name='wm') # horizontal averaged z-axis velocity, optional, becasue wm=0
timehistory.add_task(Tm, name='Tm') # horizontal averaged temperature
timehistory.add_task(Sm, name='Sm') # horizontal averaged sanility
timehistory.add_task(dz(Sm), name='S_grad') # horizontally averaged profiles of scalar gradients
timehistory.add_task(dz(Tm), name='T_grad')
timehistory.add_task(h_mean(U@ez*S), name='S_conv_flux') # horizontally averaged profiles of scalar convective fluxes
timehistory.add_task(h_mean(U@ez*T), name='T_conv_flux')
timehistory.add_task(h_mean(U@ez*S)/dz(Sm), name='S_tur_diff') # horizontally averaged profiles of turbulent diffusivity
timehistory.add_task(h_mean(U@ez*T)/dz(Tm), name='T_tur_diff')
# timehistory.add_task(h_mean(U@ez*Rho)/dz(Rhom), name='den_tur_diff')
timehistory.add_task(vol_avg(np.sqrt(upx*upx)), name='urms')
timehistory.add_task(vol_avg(np.sqrt(upz*upz)), name='wrms')
timehistory.add_task(vol_avg(np.sqrt(Tp*Tp)), name='Trms')
timehistory.add_task(vol_avg(np.sqrt(Sp*Sp)), name='Srms')

# timehistory.add_task(vol_avg((U@ez)*S)/(KS*dT/Lz), name='nusselt_s')
# timehistory.add_task(vol_avg((U@ez)*S)/(KT*dT/Lz), name='nusselt_t')
# timehistory.add_task(up@ex*up@ex*Lz/nu, name='Re_x')
# timehistory.add_task(up@ez*up@ez*Lz/nu, name='Re_z')
# timehistory.add_task(-vol_avg(U@ez*Rho), name='conv_den_flux')

# second, we save horizontal averaged fields (1D array) and coressponding pluctuation fields (2D array) following higher output time-step  
# this can help us to construct full fields for preocessing later
snapshot = solver.evaluator.add_file_handler('snapshots', sim_dt=50.0, max_writes=10000, mode=file_handler_mode)
snapshot.add_task(umx, name='um')
snapshot.add_task(umz, name='wm') #optional, becasue wm=0
snapshot.add_task(upx, name='up')
snapshot.add_task(upz, name='wp')
snapshot.add_task(Tm, name='Tm')
snapshot.add_task(Tp, name='Tp')
snapshot.add_task(Sm, name='Sm')
snapshot.add_task(Sp, name='Sp')
snapshot.add_task(p, name='p')
snapshot.add_task(-d3.div(d3.skew(U)), name='vorticity')

# store data to restart later
checkpoints = solver.evaluator.add_file_handler('checkpoints', sim_dt=100, max_writes=1, mode=file_handler_mode)
checkpoints.add_tasks(solver.state)

# CFL
CFL = d3.CFL(solver, initial_timestep, cadence=10,
             max_dt=max_timestep,min_dt = 1e-6,
             safety=0.2, threshold=0.1,
             max_change=1.5, min_change=0.5
             )
CFL.add_velocity(U)

# xg = xbasis.global_grid(dist, scale=dealias)
# zg = zbasis.global_grid(dist, scale=dealias)
# Tmg = Tm.allgather_data('g')
# umg = um.allgather_data('g')
# Main loop
try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        solver.step(timestep)  
        if (solver.iteration-1) % 5000 == 0:
            logger.info('Completed iteration {}, time={:.3f}, dt={:.10f}'.format(solver.iteration, solver.sim_time, timestep))        
            ########################### <--- plot instantaneous temperature distribution
            # Tg = Tm.allgather_data('g')+Tp.allgather_data('g')
            # Sg = Sm.allgather_data('g')+Sp.allgather_data('g')
            # umg = um.allgather_data('g')
            # if dist.comm.rank == 0:
            #     plt.figure(figsize=(3,3))
            #     plt.plot(umg[0].T,np.linspace(0,1.,int(Nz*dealias)))
            #     plt.savefig('fluc/um_time={:010.3f}.png'.format(solver.sim_time),dpi=200, bbox_inches='tight')
            #     plt.close()
            #     # plot temperature distribution
            #     plt.figure(figsize=(10,3))
            #     plt.pcolormesh(xg.ravel(),zg.ravel(),Tg.transpose(),cmap='jet',shading='gouraud')
            #     plt.colorbar() 
            #     plt.xticks([0,Lx])
            #     plt.yticks([0,Lz])
            #     plt.xlabel(r'$x$')
            #     plt.ylabel(r'$z$')
            #     plt.title("t = {:.3f}".format(solver.sim_time))
            #     # plt.show()
            #     plt.savefig('fluc/T_time={:010.3f}.png'.format(solver.sim_time),dpi=200, bbox_inches='tight')
            #     plt.close()
            #     # plot temperature distribution
            #     plt.figure(figsize=(10,3))
            #     plt.pcolormesh(xg.ravel(),zg.ravel(),Sg.transpose(),cmap='jet',shading='gouraud')
            #     plt.colorbar() 
            #     plt.xticks([0,Lx])
            #     plt.yticks([0,Lz])
            #     plt.xlabel(r'$x$')
            #     plt.ylabel(r'$z$')
            #     plt.title("t = {:.3f}".format(solver.sim_time))
            #     # plt.show()
            #     plt.savefig('fluc/S_time={:010.3f}.png'.format(solver.sim_time),dpi=200, bbox_inches='tight')
            #     plt.close()
            #     # plot horizontaly averaged density profiles
            #     meanT = np.mean(Tg,axis=0,keepdims=True)
            #     meanS = np.mean(Sg,axis=0,keepdims=True)
            #     meanDensity = (meanS-meanT)+(1-Rp)*np.linspace(0,1.,int(Nz*dealias))
            #     plt.figure(figsize=(5,5))
            #     plt.plot(np.copy(meanDensity[0]),np.linspace(0,1.,int(Nz*dealias)))
            #     plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)# Hide xticks
            #     plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
            #     plt.title("t = {:.3f}".format(solver.sim_time))
            #     # plt.show()
            #     plt.savefig('fluc/meanDensity_time={:010.3f}.png'.format(solver.sim_time),dpi=200, bbox_inches='tight')
            #     plt.close()
            ###########################
            if math.isnan(np.max(Tg)):
                logger.error('NaN values')
                break
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()
