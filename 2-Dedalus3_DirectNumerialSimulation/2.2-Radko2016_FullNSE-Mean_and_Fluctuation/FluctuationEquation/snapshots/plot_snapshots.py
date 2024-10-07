"""
Plot 2D cartesian snapshots.
mpiexec -n 24 python3 plot_snapshots.py ./*.h5
Usage:
    plot_snapshots.py <files>... [--output=<dir>]

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
        for index in range(start, start+count):
            Tm = file['tasks']['Tm'][index]
            Tp = file['tasks']['Tp'][index]
            Sm = file['tasks']['Sm'][index]
            Sp = file['tasks']['Sp'][index]
            
            Density = (Sm+Sp-Tm-Tp)+(1-Rp)*z
            meanDensity = (Sm-Tm)+(1-Rp)*z

            # plot Density
            plt.figure(figsize=(10,3))
            plt.pcolormesh(x,z,Density.T,shading='gouraud')
            plt.set_cmap('jet')
            plt.colorbar()
            plt.xlabel(r'$x$')
            plt.ylabel(r'$z$')
            time = title_func(file['scales/sim_time'][index])
            plt.title(time)
            plt.xticks([0,Lx])
            plt.yticks([0,Lz])
            plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
            plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
            plt.savefig('frames/Rho_time={:010.3f}.png'.format(file['scales/sim_time'][index]), dpi=dpi, bbox_inches='tight')
            plt.close()

            # plot T
            plt.figure(figsize=(10,3))
            plt.pcolormesh(x,z,(Tm+Tp).transpose(),shading='gouraud')
            plt.set_cmap('jet')
            plt.colorbar()
            plt.xlabel(r'$x$')
            plt.ylabel(r'$z$')
            time = title_func(file['scales/sim_time'][index])
            plt.title(time)
            plt.xticks([0,Lx])
            plt.yticks([0,Lz])
            plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
            plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
            plt.savefig('frames/T_time={:010.3f}.png'.format(file['scales/sim_time'][index]), dpi=dpi, bbox_inches='tight')
            plt.close()

            # plot S
            plt.figure(figsize=(10,3))
            plt.pcolormesh(x,z,(Sm+Sp).transpose(),shading='gouraud')
            plt.set_cmap('jet')
            plt.colorbar()
            plt.xlabel(r'$x$')
            plt.ylabel(r'$z$')
            plt.title(time)
            plt.xticks([0,Lx])
            plt.yticks([0,Lz])
            plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)# Hide xticks
            plt.tick_params(axis='y', which='both', left=True, right=False, labelleft=True)
            plt.savefig('frames/S_time={:010.3f}.png'.format(file['scales/sim_time'][index]), dpi=dpi, bbox_inches='tight')
            plt.close()

            plt.figure(figsize=(3,3))
            plt.plot(meanDensity[0],z)
            plt.title(time)
            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)# Hide xticks
            plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
            plt.savefig('frames/meanDensity_time={:010.3f}.png'.format(file['scales/sim_time'][index]), dpi=dpi, bbox_inches='tight')
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