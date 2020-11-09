import argparse
import vtk_reader
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

try:
    plt.style.use('mystyle')
except OSError:
    print('...')


def show_grid(fname, save, L=1):
    (u,v,w), (p), (x,y,z) = vtk_reader.read_vtr(fname)
    x=(x-0.5)/L; y=(y-0.5)/L; z=(z-0.5)/L
    if len(z)==1:
        fig, ax = plot_2D_Grid(x, y, p[:,:,0], L, every=2)
    else:
        fig, ax = plot_3D_Gird(x, y, z, p, L, every=2)
    plt.tight_layout()
    if save: plt.savefig("grid.png", dpi=300)
    plt.show()
    

def plot_XY_grid(ax, x, y, dist, L, every, xlab=r"x/L", ylab=r"y/L"):
    ax.vlines(x[::every], ymin=min(y), ymax=max(y), color='gray', lw=0.2)
    ax.hlines(y[::every], xmin=min(x), xmax=max(x), color='gray', lw=0.2)
    ax.contour(x, y, dist.T, colors='k', levels=[0.])
    ax.set_xlabel(xlab); ax.set_ylabel(ylab)
    ax.set_xlim(min(x-50/L),max(x+50/L))
    ax.set_ylim(min(y-50/L),max(y+50/L))
    ax.set_aspect('equal', adjustable='datalim', anchor='C')
    return ax


def plot_2D_Grid(x, y, dist, L, every):
    fig = plt.figure()
    ax = plt.subplot2grid((1,1), (0,0))
    ax = plot_XY_grid(ax, x, y, dist, L, every)
    return fig, ax


def plot_3D_Gird(x, y, z, dist, L, every):
    fig = plt.figure(figsize=(10,6))
    xz = plt.subplot2grid((2,3), (0,0), colspan=2)
    xy = plt.subplot2grid((2,3), (1,0), colspan=2)
    yz = plt.subplot2grid((2,3), (0,2))
    xy = plot_XY_grid(xy, x, y, dist[:,:,np.argmin(abs(z))], L, every)
    yz = plot_XY_grid(yz, y, z, dist[np.argmin(abs(x)),:,:], L, every, xlab=r"y/L", ylab=r"z/L")
    xz = plot_XY_grid(xz, x, z, dist[:,np.argmin(abs(y)),:], L, every, ylab=r"z/L")
    return fig, (xy, yz, xz)


# def plot_3D_Gird(x, y, z, dist, every=2):
#     fig = plt.figure()
#     ax = fig.gca(projection='3d')
#     for i in range(0,len(x),every):
#         ax.plot([x[i],x[i]], [min(y),max(y)], [min(z),min(z)], color='gray', lw=0.2)
#         ax.plot([x[i],x[i]], [max(y),max(y)], [min(z),max(z)], color='gray', lw=0.2)
#     for j in range(0,len(y),every):
#         ax.plot([min(x),min(x)], [y[j],y[j]], [min(z),max(z)], color='gray', lw=0.2)
#         ax.plot([min(x),max(x)], [y[j],y[j]], [min(z),min(z)], color='gray', lw=0.2)
#     for k in range(0,len(z),every):
#         ax.plot([min(x),min(x)], [min(y),max(y)], [z[k],z[k]], color='gray', lw=0.2)
#         ax.plot([min(x),max(x)], [max(y),max(y)], [z[k],z[k]], color='gray', lw=0.2)
#     ax.set_xlabel(r"x/L"); ax.set_ylabel(r"y/L"); ax.set_zlabel(r"z/L")
#     ax.grid(False)
#     return fig, ax


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Specify Project to use to plot grid.')
    parser.add_argument('-f','--filename', help='Name of the Lotus project to use')
    parser.add_argument('-L','--lengthscale', help='Length scale of the problem')
    args = parser.parse_args()
    if (not args.lengthscale):
        L = 1.
    else:
         L = float(args.lengthscale)
    show_grid(args.filename+"/datp/bodyF.-1.pvtr", True, L)
