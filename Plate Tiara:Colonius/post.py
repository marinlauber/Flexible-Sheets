#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

try:
    plt.style.use('mystyle')
except OSError:
    print('OSError: using default style')

if __name__ == "__main__":
    
    lab = [rf'avg $C_D$', rf'avg $C_L$']
    cols = ['-ok', '-sk', '-^k']
    lgd = ['AR : 1','AR : 2','AR : 4',]

    for k,file in enumerate(['CD', 'CL']):
        dat = np.genfromtxt('plate_'+file+'.csv', delimiter=',')
        plt.figure(figsize=(4,4))
        for i in range(1,4):
            plt.plot(dat[:,0], dat[:,i], cols[i-1], lw=0.5, mfc='white', label=lgd[i-1])
        if k==0:
            plt.plot(30,0.6684, 'sr', mfc='None')
        else:
            plt.plot(30,0.5881, 'sr', mfc='None')
        plt.xlabel(rf'$\alpha$ $^\circ$')
        plt.ylabel(lab[k])
        plt.xlim(0,60); plt.ylim(0,1.5)
        if k==0:
            plt.legend()
        plt.savefig(file+'_plate.png', dpi=300)
        plt.show()
