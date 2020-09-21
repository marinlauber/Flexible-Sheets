#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
try:
    plt.style.use('mystyle')
except OSError:
    print('nope...')


def read_lotus(fname):
    df = pd.read_csv(fname,delim_whitespace = True,
            names=["time","CFL","pforce_x","pforce_y","pforce_z",
                "vforce_x","vforce_y","vforce_z",
                "phi","theta","alpha",
                "dphi","dtheta","dalpha"])
    df['Fx'] = df.pforce_x + df.vforce_x
    df['Fy'] = df.pforce_y + df.vforce_y
    df['Fz'] = df.pforce_z + df.vforce_z
    df['Lift'] = df.Fy
    df['Drag'] = df.Fx*np.cos(df.phi) - df.Fz*np.sin(df.phi)

    return df

if __name__ == "__main__":
    
    lab = [rf'Bose et al. (2013)',rf'Zheng et al. (2020)']
    m = ['-^','-s']
    ylab = [r'Lift Coefficient', r'Drag Coefficient']

    for k,force in enumerate(['CL','CD']):
        plt.figure(figsize=(6,4))
        for i,author in enumerate(['bose', 'Zheng']):
            dat = np.genfromtxt('Dat/'+force+'_'+author+'.txt', delimiter=',')
            plt.plot(dat[:,0], dat[:,1], m[i], mfc='None', label=lab[i])
        for i,f in enumerate(['new_kernel_fort.9']):
            df = read_lotus(f)
            labels = 'BDIM'
            if i==0:
                labels = 'BDIM new kernel'
            if force=='CL':
                plt.plot(df.time, df.Lift, '-o', mfc='None', markevery=10, label=labels)
            else:
                plt.plot(df.time, df.Drag, '-o', mfc='None', markevery=10, label=labels)
        plt.xlim(0,1)
        plt.xlabel(r'$t/T$')
        plt.ylabel(ylab[k])
        if k==1:
            plt.legend()
        plt.tight_layout()
        plt.savefig(force+'_wing.png', dpi=900)
        plt.show()

