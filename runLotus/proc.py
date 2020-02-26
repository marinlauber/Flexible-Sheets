#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def read_it(folder):
    df = pd.read_csv(folder+'/fort.9',delim_whitespace = True,
        names=["time","CFL","drag","tke","s","U","a"])
    df.drop(df.index[:3], inplace=True)
    return df

for accel in [0.5,2.0]:
    for reynolds in [2,3,4]:
        folder = 'acc{:02d}_rey5e{:01d}'.format(int(10*accel),reynolds)
        df = read_it(folder)
        plt.plot(df.s,df.drag,label=folder)
plt.xlabel('s/D')
plt.ylabel('drag')
plt.legend()
plt.savefig('compare.png')
