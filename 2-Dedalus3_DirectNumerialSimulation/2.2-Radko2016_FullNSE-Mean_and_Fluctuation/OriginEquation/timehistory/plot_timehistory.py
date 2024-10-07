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
    S = np.zeros((Nx,Nz))

    # Plot writes
    with h5py.File(filename, mode='r') as file:
        t_list = []
        um_list = [] # horizontal averaged x-axis velocity
        Tm_list = [] # horizontal averaged temperature
        Sm_list = [] # horizontal averaged sanility
        Rhom_list = [] # horizontal averaged density
        dzSm_list = [] # horizontally averaged profiles of sanility gradients
        dzTm_list = [] # horizontally averaged profiles of temperature gradients
        S_conv_flux_list = [] # horizontally averaged profiles of convective sanility fluxes
        T_conv_flux_list = [] # horizontally averaged profiles of convective temperature fluxes
        S_tur_diff_list = [] # horizontally averaged profiles of turbulent sanility diffusivity
        T_tur_diff_list = [] # horizontally averaged profiles of turbulent temperature diffusivity
        
        urms_list = []
        wrms_list = []
        Trms_list = []
        Srms_list = []

        for index in range(start, start+count):
            t_list.append(file['scales/sim_time'][index])

            um = file['tasks']['um'][index]
            Tm = file['tasks']['Tm'][index]
            Sm = file['tasks']['Sm'][index]
            Rhom = (Sm-Tm)+(1-Rp)*z
            dzTm = file['tasks']['T_grad'][index]
            dzSm = file['tasks']['S_grad'][index]
            T_conv_flux = file['tasks']['T_conv_flux'][index]
            S_conv_flux = file['tasks']['S_conv_flux'][index]
            T_tur_diff = file['tasks']['T_tur_diff'][index]
            S_tur_diff = file['tasks']['S_tur_diff'][index]

            urms = file['tasks']['urms'][index]
            wrms = file['tasks']['wrms'][index]
            Trms = file['tasks']['Trms'][index]
            Srms = file['tasks']['Srms'][index]

            um_list.append(um.T)
            Tm_list.append(Tm.T)
            Sm_list.append(Sm.T)
            Rhom_list.append(Rhom.T)
            dzTm_list.append(dzTm.T)
            dzSm_list.append(dzSm.T)
            T_conv_flux_list.append(T_conv_flux.T)
            S_conv_flux_list.append(S_conv_flux.T)
            T_tur_diff_list.append(T_tur_diff.T)
            S_tur_diff_list.append(S_tur_diff.T)
            
            urms_list.append(urms)
            wrms_list.append(wrms)
            Trms_list.append(Trms)
            Srms_list.append(Srms)

        # plot horizontal averaged density
        plt.figure(figsize=(10,3))
        plt.pcolormesh(np.array(t_list),z,np.array(Rhom_list)[:,:,0].T,shading='gouraud')
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
        plt.pcolormesh(np.array(t_list),z,np.array(um_list)[:,:,0].T,shading='gouraud')
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
        plt.pcolormesh(np.array(t_list),z,np.array(Tm_list)[:,:,0].T,shading='gouraud')
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
        plt.pcolormesh(np.array(t_list),z,np.array(Sm_list)[:,:,0].T,shading='gouraud')
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
        plt.pcolormesh(np.array(t_list),z,np.array(dzTm_list)[:,:,0].T,shading='gouraud')
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
        plt.pcolormesh(np.array(t_list),z,np.array(dzSm_list)[:,:,0].T,shading='gouraud')
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
        plt.pcolormesh(np.array(t_list),z,np.array(S_conv_flux_list)[:,:,0].T,shading='gouraud')
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
        plt.pcolormesh(np.array(t_list),z,np.array(T_conv_flux_list)[:,:,0].T,shading='gouraud')
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
        plt.pcolormesh(np.array(t_list),z,np.array(S_tur_diff_list)[:,:,0].T,shading='gouraud')
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
        plt.pcolormesh(np.array(t_list),z,np.array(T_tur_diff_list)[:,:,0].T,shading='gouraud')
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
        plt.pcolormesh(np.array(t_list),np.array(T_tur_diff_list)[:,0].T,shading='gouraud')
        plt.xlabel(r'$t$')
        plt.xticks([0,max(t_list)])
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
        plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
        plt.savefig('frames/urms.png', dpi=dpi, bbox_inches='tight')
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