"""
Plot 2D cartesian snapshots.

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
    Lz = 1
    x = np.linspace(0,64.,Nx)
    z = np.linspace(0,1.,Nz)
    S = np.zeros((Nx,Nz))

    # Plot writes
    with h5py.File(filename, mode='r') as file:
        for index in range(start, start+count):
            dset1 = file['tasks']['temperature']
            dset2 = file['tasks']['salinity']
            
            
            meanT = np.mean(dset1[index],axis=0,keepdims=True)
            meanS = np.mean(dset2[index],axis=0,keepdims=True)
            meanDensity = (meanS-meanT)+(1-Rp)*np.linspace(0,1.,Nz)

            # plot T
            plt.figure(figsize=(10,3))
            plt.pcolormesh(x,z,dset1[index].transpose(),shading='nearest')
            plt.set_cmap('turbo')
            plt.colorbar()
            plt.xlabel(r'$x$')
            plt.ylabel(r'$z$')
            title = title_func(file['scales/sim_time'][index])
            plt.title(title)
            # plt.xticks([])
            # plt.yticks([])
            # Save figure
            # savename = savename_func(file['scales/write_number'][index])
            # savepath = output.joinpath(savename)
            plt.savefig('frames/T_time={:.3f}.png'.format(file['scales/write_number'][index]), dpi=dpi)
            plt.close()

            plt.figure(figsize=(3,3))
            plt.plot(meanDensity[0])
            # plt.colorbar()
            # plt.xlabel(r'$x$')
            # plt.ylabel(r'$z$')
            title = title_func(file['scales/sim_time'][index])
            plt.title(title)
            plt.xticks(None)
            plt.yticks(None)
            plt.xlabel(None)
            plt.ylabel(None)
            # Save figure
            # savename = savename_func(file['scales/write_number'][index])
            # savepath = output.joinpath(savename)
            plt.savefig('frames/meanDensity_time={:.3f}.png'.format(file['scales/write_number'][index]), dpi=dpi)
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