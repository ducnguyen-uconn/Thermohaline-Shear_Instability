"""
Plot 2D cartesian snapshots.
mpiexec -n 1 python3 plot_timehistory.py ./*.h5
Usage:
    plot_timehistory.py <files>... [--output=<dir>]

Options:
    --output=<dir>  Output directory [default: ./frames]

"""

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dedalus.extras import plot_tools



def main(filename, start, count, output):
    """Save plot of specified tasks for given range of analysis writes."""

    # Plot settings
    scale = 1.5
    dpi = 200
    title_func = lambda sim_time: 't = {:.3f}'.format(sim_time)
    savename_func = lambda write: 'write_{:06}.png'.format(write)

    Rp = 2
    Nx, Nz = 384, 192
    kx = 8
    Lx, Lz = 64, 1
    x = np.linspace(0,Lx,Nx)
    z = np.linspace(0,Lz,Nz)
    

    # Plot writes
    with h5py.File(filename, mode='r') as file:
        t_list = np.array(file['scales/sim_time'])
        # Create the meshgrid for plotting 
        X, Z = np.meshgrid(t_list, z)

        um_list = np.array(file['tasks']['um'])[:,0,:].T # horizontal averaged x-axis velocity
        Tm_list = np.array(file['tasks']['Tm'])[:,0,:].T # horizontal averaged temperature
        Sm_list = np.array(file['tasks']['Sm'])[:,0,:].T # horizontal averaged sanility
        Rhom_list = (Sm_list-Tm_list)+(1-Rp)*Z # horizontal averaged density
        dzSm_list = np.array(file['tasks']['T_grad'])[:,0,:].T # horizontally averaged profiles of sanility gradients
        dzTm_list = np.array(file['tasks']['S_grad'])[:,0,:].T # horizontally averaged profiles of temperature gradients
        S_conv_flux_list = np.array(file['tasks']['T_conv_flux'])[:,0,:].T # horizontally averaged profiles of convective sanility fluxes
        T_conv_flux_list = np.array(file['tasks']['S_conv_flux'])[:,0,:].T # horizontally averaged profiles of convective temperature fluxes
        S_tur_diff_list = T_conv_flux_list/dzTm_list # horizontally averaged profiles of turbulent sanility diffusivity
        T_tur_diff_list = S_conv_flux_list/dzSm_list # horizontally averaged profiles of turbulent temperature diffusivity

        urms_list = np.array(file['tasks']['urms'])[:,0,:]
        wrms_list = np.array(file['tasks']['wrms'])[:,0,:]
        Trms_list = np.array(file['tasks']['Trms'])[:,0,:]
        Srms_list = np.array(file['tasks']['Srms'])[:,0,:]

        # plot horizontal averaged density
        plt.figure(figsize=(10,3))
        plt.pcolormesh(t_list,z,Rhom_list,shading='gouraud')
        plt.set_cmap('Purples')
        plt.colorbar()
        plt.xlabel(r'$t$')
        plt.ylabel(r'$z$')
        plt.xticks([0,max(t_list)])
        plt.yticks([0,Lz])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/Rhom.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # plot horizontal averaged x-axis velocity
        plt.figure(figsize=(10,3))
        plt.pcolormesh(t_list,z,um_list,shading='gouraud')
        plt.set_cmap('PRGn')
        plt.colorbar()
        plt.xlabel(r'$t$')
        plt.ylabel(r'$z$')
        plt.xticks([0,max(t_list)])
        plt.yticks([0,Lz])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/um.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # plot horizontal averaged temperature
        plt.figure(figsize=(10,3))
        plt.pcolormesh(t_list,z,Tm_list,shading='gouraud')
        plt.set_cmap('PRGn')
        plt.colorbar()
        plt.xlabel(r'$t$')
        plt.ylabel(r'$z$')
        plt.xticks([0,max(t_list)])
        plt.yticks([0,Lz])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/Tm.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # plot horizontal averaged sanility
        plt.figure(figsize=(10,3))
        plt.pcolormesh(t_list,z,Sm_list,shading='gouraud')
        plt.set_cmap('PRGn')
        plt.colorbar()
        plt.xlabel(r'$t$')
        plt.ylabel(r'$z$')
        plt.xticks([0,max(t_list)])
        plt.yticks([0,Lz])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/Sm.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # plot horizontally averaged profiles of temperature gradients
        plt.figure(figsize=(10,3))
        plt.pcolormesh(t_list,z,dzTm_list,shading='gouraud')
        plt.set_cmap('PRGn')
        plt.colorbar()
        plt.xlabel(r'$t$')
        plt.ylabel(r'$z$')
        plt.xticks([0,max(t_list)])
        plt.yticks([0,Lz])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/dzTm.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # plot horizontally averaged profiles of sanility gradients
        plt.figure(figsize=(10,3))
        plt.pcolormesh(t_list,z,dzSm_list,shading='gouraud')
        plt.set_cmap('PRGn')
        plt.colorbar()
        plt.xlabel(r'$t$')
        plt.ylabel(r'$z$')
        plt.xticks([0,max(t_list)])
        plt.yticks([0,Lz])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/dzSm.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # plot horizontally averaged profiles of convective sanility fluxes
        plt.figure(figsize=(10,3))
        plt.pcolormesh(t_list,z,S_conv_flux_list,shading='gouraud')
        plt.set_cmap('Purples')
        plt.colorbar()
        plt.xlabel(r'$t$')
        plt.ylabel(r'$z$')
        plt.xticks([0,max(t_list)])
        plt.yticks([0,Lz])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/S_conv_flux.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # plot horizontally averaged profiles of convective temperature fluxes
        plt.figure(figsize=(10,3))
        plt.pcolormesh(t_list,z,T_conv_flux_list,shading='gouraud')
        plt.set_cmap('Reds')
        plt.colorbar()
        plt.xlabel(r'$t$')
        plt.ylabel(r'$z$')
        plt.xticks([0,max(t_list)])
        plt.yticks([0,Lz])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/T_conv_flux.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # plot horizontally averaged profiles of turbulent sanility diffusivity
        plt.figure(figsize=(10,3))
        plt.pcolormesh(t_list,z,S_tur_diff_list,shading='gouraud')
        plt.set_cmap('PRGn')
        plt.colorbar()
        plt.xlabel(r'$t$')
        plt.ylabel(r'$z$')
        plt.xticks([0,max(t_list)])
        plt.yticks([0,Lz])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/S_tur_diff.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        # plot horizontally averaged profiles of turbulent temperature diffusivity
        plt.figure(figsize=(10,3))
        plt.pcolormesh(t_list,z,T_tur_diff_list,shading='gouraud')
        plt.set_cmap('PRGn')
        plt.colorbar()
        plt.xlabel(r'$t$')
        plt.ylabel(r'$z$')
        plt.xticks([0,max(t_list)])
        plt.yticks([0,Lz])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/T_tur_diff.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(10,2))
        plt.plot(t_list,urms_list)
        plt.xlabel(r'$t$')
        plt.xticks([0,max(t_list)])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/urms.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(10,2))
        plt.plot(t_list,wrms_list)
        plt.xlabel(r'$t$')
        plt.xticks([0,max(t_list)])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/wrms.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(10,2))
        plt.plot(t_list,Trms_list)
        plt.xlabel(r'$t$')
        plt.xticks([0,max(t_list)])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/Trms.png', dpi=dpi, bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(10,2))
        plt.plot(t_list,Srms_list)
        plt.xlabel(r'$t$')
        plt.xticks([0,max(t_list)])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/Srms.png', dpi=dpi, bbox_inches='tight')
        plt.close()


if __name__ == "__main__":

    import pathlib
    from docopt import docopt
    from dedalus.tools import logging
    from dedalus.tools import post
    from dedalus.tools.parallel import Sync

    args = docopt(__doc__)

    output_path = pathlib.Path(args['--output']).absolute()
    # Create output directory if needed
    with Sync() as sync:
        if sync.comm.rank == 0:
            if not output_path.exists():
                output_path.mkdir()
    post.visit_writes(args['<files>'], main, output=output_path)