import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import binned_statistic

plt.style.use('mystyle')

def get_data(fname):
    return np.genfromtxt(fname+'.txt', skip_header=5)

def stats(x, y, nbins=16):
    mean,xs,_  = binned_statistic(x, y, 'mean', bins=nbins)
    std,xs,_ = binned_statistic(x, y, 'std', bins=nbins)
    mn,xs,_ = binned_statistic(x, y, 'min', bins=nbins)
    mx,xs,_ = binned_statistic(x, y, 'max', bins=nbins)
    xs = (xs[1:] + xs[:-1])/2
    return xs, mean, std, mn, mx

def read_it(folder):
    df = pd.read_csv(folder+'/fort.9', delim_whitespace=True,
                     names=["time","CFL","pforce_x","pforce_y","pforce_z","vforce_x","vforce_y","vforce_z"])
    df.drop(df.index[:3], inplace=True)
    return df


def get_forces(a, string,taverage=10):
    b = np.empty((2, len(a)))
    for i, j in enumerate(a):
        folder = string+str(j)
        df = read_it(folder)
        df['force_x'] = df.pforce_x+df.vforce_x
        df['force_y'] = df.pforce_y+df.vforce_y
        b[0, i] = df.force_x[df.time>taverage].mean()
        b[1, i] = df.force_y[df.time>taverage].mean()
    return b


dat1 = get_data('Souppez')
dat2 = get_data('Velychko')
dat3 = get_data('Bot_Re_68k')
dat4 = get_data('Bot_Re_238k')
dat5 = get_data('Bot_Re_369k')
dat6 = get_data('Collie')

# numerical data
c3d = get_forces(['256_3D_3','384_3D_1'],'../../projects/arc/C',taverage=40)
c3d_10 = get_forces(['384_12_3D','384_12_3D'],'../../projects/arc/C',taverage=40)
c3d[:,1] /= 8.
c3d_10 /= 8.
alpha_n = [18, 18]

# low RE data
alpha = np.hstack((dat1[:,0],dat2[:15,0],dat3[:,0]))
cl = np.hstack((dat1[:,1],dat2[:15,1],dat3[:,1]))
alpha_2 = np.hstack((dat1[:,2],dat2[:15,2],dat3[:,2]))
cd = np.hstack((dat1[:,3],dat2[:15,3],dat3[:,3]))

fig, ax = plt.subplots(1,2,figsize=(10,6))

x_, mean, std, mn, mx = stats(alpha, cl)
ax[0].plot(alpha, cl, '.',  color='black', alpha=0.15)
ax[0].fill_between(x_, mx, mn, alpha=0.05, color='C0')
ax[0].plot(x_, mean, 'C0', lw=2)
ax[0].plot(alpha_n,c3d[1,:],'xr')
ax[0].plot([12,12],c3d_10[1,:],'xr')

x_, mean, std, mn, mx = stats(alpha_2, cd)
p1, = ax[1].plot(alpha_2, cd, '.',  color='black', alpha=0.15)
p2 = ax[1].fill_between(x_, mx, mn, alpha=0.05, color='C0')
p3, = ax[1].plot(x_, mean, 'C0', lw=2)
p4, = ax[1].plot([0], marker='None',linestyle='None')
pp, = ax[1].plot(alpha_n,c3d[0,:],'xr')
ax[1].plot([12,12],c3d_10[0,:],'xr')
# High Re data
alpha = np.hstack((dat4[:,0],dat5[:,0],dat6[:,0]))
cl = np.hstack((dat4[:,1],dat5[:,1],dat6[:,1]))
alpha_2 = np.hstack((dat4[:,2],dat5[:,2],dat6[:,2]))
cd = np.hstack((dat4[:,3],dat5[:,3],dat6[:,3]))

x_, mean, std, mn, mx = stats(alpha, cl)
ax[0].plot(alpha, cl, 'x', color='black', alpha=0.15)
ax[0].fill_between(x_, mx,  mn, alpha=0.05, color='C1')
ax[0].plot(x_, mean, '--C1', lw=2)

x_, mean, std, mn, mx = stats(alpha_2, cd)
p5, = ax[1].plot(alpha_2, cd, 'x', color='black', alpha=0.15)
p6 = ax[1].fill_between(x_, mx, mn, alpha=0.05, color='C1')
p7, = ax[1].plot(x_, mean, '--C1', lw=2)
p8, = ax[1].plot([0], marker='None',linestyle='None')

ax[0].set_xlabel(r'Angle of attack $\alpha$ ($^\circ$)')
ax[1].set_xlabel(r'Angle of attack $\alpha$ ($^\circ$)')
ax[0].set_ylabel(r'$C_L$'); ax[1].set_ylabel(r'$C_D$')
leg = ['Experimetal Data',r'$\pm$ Envelope','Mean']
ax[1].legend([p4,p1,p2,p3,pp,p8,p5,p6,p7],
             [r'$\underline{Re : 53-68k}$']+leg+['3D BDIM']+[r'$\underline{Re : 238-525k}$']+leg,
              ncol=2)
plt.tight_layout()
plt.savefig('experimemntal_data.png',dpi=300)
plt.show()
