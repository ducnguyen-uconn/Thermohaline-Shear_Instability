'''
Dedalus script to simulate thermohaline-shear equations in 2D periodic domain.
um, wm, Tm, Sm: <u>, <w>, <T>, <S> (<.> is horizontally averaged variable)
up, wp, Tp, Sp: u', w', T', S'

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
stop_sim_time = 100 + 300*restart # Stopping criteria
# Bases
coords = d3.CartesianCoordinates('x','z')
dist = d3.Distributor(coords, dtype=np.float64)
# define the coordinate system
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
zbasis = d3.RealFourier(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)

# define fields
# um, wm, pm, Tm, Sm: <u>, <w>, <p>, <T>, <S> (<.> is horizontally averaged variable)
# up, wp, pp, Tp, Sp:  u',  w',  p',  T',  S'
pm = dist.Field(name='pm', bases=(zbasis)) # pressure
pp = dist.Field(name='pp', bases=(xbasis,zbasis)) # pressure
um = dist.VectorField(coords,name='um', bases=(zbasis)) # mean u(z)
up = dist.VectorField(coords,name='up', bases=(xbasis,zbasis)) # fluctuation u'(x,z)
# wm = dist.Field(name='wm', bases=(zbasis)) # mean w(z)
# wp = dist.Field(name='wp', bases=(xbasis,zbasis)) # fluctuation w'(x,z)
Sm = dist.Field(name='Sm', bases=(zbasis)) # mean S(z)
Sp = dist.Field(name='Sp', bases=(xbasis,zbasis)) # fluctuation S'(x,z)
Tm = dist.Field(name='Tm', bases=(zbasis)) # mean T(z)
Tp = dist.Field(name='Tp', bases=(xbasis,zbasis)) # fluctuation T'(x,z)
Ubg = dist.Field(bases=(zbasis)) #background velocity of basic state

# Substitutions
x, z = dist.local_grids(xbasis, zbasis) # get coordinate arrays in horizontal and vertical directions
ex, ez = coords.unit_vector_fields(dist) # get unit vectors in horizontal and vertical directions
# define velocity components
U = um + up # full velocity
T = Tm + Tp # full temperature
S = Sm + Sp # full salinity

# create constant sub-field for incompressible flow condition's equation
tau_pm = dist.Field(name='tau_pm') 
tau_pp = dist.Field(name='tau_pp') 
# tau_wm = dist.Field(name='tau_wm') 
# because this term is only a contant added to the equation, we don't need to instantiate it for bases system

# grad_te = d3.grad(te) # First-order reduction
# grad_sa = d3.grad(sa) # First-order reduction
# grad_u = d3.grad(u) # First-order reduction

# lap_u = d3.div(grad_u)
# lap_te = d3.div(grad_te)
# lap_sa = d3.div(grad_sa)
# First-order form: "div(A)" becomes "trace(grad_A)"

# grad = lambda A: d3.grad(A)
lap = lambda A: d3.div(d3.grad(A)) # First-order form: "lap(f)" becomes "div(grad_f)"
dx = lambda A: d3.Differentiate(A, coords['x']) 
dz = lambda A: d3.Differentiate(A, coords['z']) 
h_mean = lambda A: d3.Integrate(A,'x')/Lx         # Horizontal mean of A

Ubg['g'] = np.sin(2*pi*z)

# Problem
problem = d3.IVP([pm,pp,tau_pm, tau_pp,
                  um,up,
                  Tm,Tp, 
                  Sm,Sp], namespace=locals())

problem.add_equation("trace(grad(um))+tau_pm = 0")
problem.add_equation("trace(grad(up))+tau_pp = 0")
problem.add_equation("integ(pm) = 0") # Pressure gauge
problem.add_equation("integ(pp) = 0") # Pressure gauge

problem.add_equation("dt(um) + Ubg*dx(um) + (um@ez)*dz(Ubg)*ex + grad(pm) - (Pr/Pe)*lap(um) - (4*pi*pi*Ri/(Rp-1))*(Tm-Sm)*ez = - dz(h_mean((up@ez)*up))")
problem.add_equation("dt(up) + Ubg*dx(up) + (up@ez)*dz(Ubg)*ex + grad(pp) - (Pr/Pe)*lap(up) - (4*pi*pi*Ri/(Rp-1))*(Tp-Sp)*ez = - up@grad(um)-um@grad(up)-up@grad(up)+dz(h_mean((up@ez)*up))")

problem.add_equation("dt(Tm) - (1./Pe)*lap(Tm) - um@ez = - dz(h_mean((up@ez)*Tp))")
problem.add_equation("dt(Tp) - (1./Pe)*lap(Tp) - up@ez = - up@grad(Tm)-um@grad(Tp)-up@grad(Tp)+dz(h_mean((up@ez)*Tp))")

problem.add_equation("dt(Sm) - (tau/Pe)*lap(Sm) - Rp*um@ez = - dz(h_mean((up@ez)*Sp))")
problem.add_equation("dt(Sp) - (tau/Pe)*lap(Sp) - Rp*up@ez = - up@grad(Sm)-um@grad(Sp)-up@grad(Sp)+dz(h_mean((up@ez)*Sp))")

# timestepper = d3.RK443
timestepper = d3.RK222

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Initial conditions
if not restart:
    pm.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    pp.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    um.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    up.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    Tm.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    Tp.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    Sm.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    Sp.fill_random('g', seed=42, distribution='normal', scale=1e-4) # Random noise
    file_handler_mode = 'overwrite'
    initial_timestep = 0.02
    max_timestep = 0.02
else:
    write, initial_timestep = solver.load_state('/home/jms24002/reproduce_radko_results/fig5/checkpoints/checkpoints_s90.h5')
    max_timestep = 0.02
    file_handler_mode = 'append'
    logger.info('Imported last-step data successfully')

# store data for analysis later
# dataset = solver.evaluator.add_file_handler('snapshots', sim_dt=50.0, max_writes=10000, mode=file_handler_mode)
# dataset.add_task(umx, name='umx')
# dataset.add_task(umz, name='umz')
# dataset.add_task(Tm, name='Tm')
# dataset.add_task(Tp, name='Tp')
# dataset.add_task(Sm, name='Sm')
# dataset.add_task(Sp, name='Sp')
# dataset.add_task(P, name='p')
# dataset.add_task(-d3.div(d3.skew(um+up)), name='vorticity')
# store data to restart later
# checkpoints = solver.evaluator.add_file_handler('checkpoints', sim_dt=100, max_writes=1, mode=file_handler_mode)
# checkpoints.add_tasks(solver.state)

# CFL
CFL = d3.CFL(solver, initial_timestep, cadence=10,
             max_dt=max_timestep,min_dt = 1e-6,
             safety=0.2, threshold=0.1,
             max_change=1.5, min_change=0.5
             )
CFL.add_velocity(U)

xg = xbasis.global_grid(dist, scale=dealias)
zg = zbasis.global_grid(dist, scale=dealias)
# Tmg = Tm.allgather_data('g')
umg = um.allgather_data('g')
# if dist.comm.rank == 0:
#     print(np.copy(umg[0]))
# Main loop
oldtime = 0.
try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        solver.step(timestep)  
        if (solver.iteration-1) % 1000 == 0:
            logger.info('Completed iteration {}, time={:.3f}, dt={:.10f}'.format(solver.iteration, solver.sim_time, timestep))        
            ########################### <--- plot instantaneous temperature distribution
            Tg = Tm.allgather_data('g')+Tp.allgather_data('g')
            Sg = Sm.allgather_data('g')+Sp.allgather_data('g')
            umg = um.allgather_data('g')
            if dist.comm.rank == 0:
                plt.figure(figsize=(3,3))
                plt.plot(umg[0].T,np.linspace(0,1.,int(Nz*dealias)))
                plt.savefig('fluc/um_time={:010.3f}.png'.format(solver.sim_time), bbox_inches='tight')
                plt.close()
                # plot temperature distribution
                plt.figure(figsize=(10,3))
                plt.pcolormesh(xg.ravel(),zg.ravel(),Tg.transpose(),cmap='jet')
                plt.colorbar() 
                plt.xticks([0,Lx])
                plt.yticks([0,Lz])
                plt.xlabel(r'$x$')
                plt.ylabel(r'$z$')
                plt.title("t = {:.3f}".format(solver.sim_time))
                # plt.show()
                plt.savefig('fluc/T_time={:010.3f}.png'.format(solver.sim_time), bbox_inches='tight')
                plt.close()
                # plot temperature distribution
                plt.figure(figsize=(10,3))
                plt.pcolormesh(xg.ravel(),zg.ravel(),Sg.transpose(),cmap='jet')
                plt.colorbar() 
                plt.xticks([0,Lx])
                plt.yticks([0,Lz])
                plt.xlabel(r'$x$')
                plt.ylabel(r'$z$')
                plt.title("t = {:.3f}".format(solver.sim_time))
                # plt.show()
                plt.savefig('fluc/S_time={:010.3f}.png'.format(solver.sim_time), bbox_inches='tight')
                plt.close()
                # plot horizontaly averaged density profiles
                meanT = np.mean(Tg,axis=0,keepdims=True)
                meanS = np.mean(Sg,axis=0,keepdims=True)
                meanDensity = (meanS-meanT)+(1-Rp)*np.linspace(0,1.,int(Nz*dealias))
                plt.figure(figsize=(5,5))
                plt.plot(np.copy(meanDensity[0]),np.linspace(0,1.,int(Nz*dealias)))
                plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)# Hide xticks
                plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
                plt.title("t = {:.3f}".format(solver.sim_time))
                # plt.show()
                plt.savefig('fluc/meanDensity_time={:010.3f}.png'.format(solver.sim_time), bbox_inches='tight')
                plt.close()
            ###########################
            if math.isnan(np.max(Tg)):
                logger.error('NaN values')
                break
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()
