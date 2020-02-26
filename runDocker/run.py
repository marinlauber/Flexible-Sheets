#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import lotus
import subprocess

def runcase(folder,dic,resume=None,n=2):
    print(dic)
    lotus.replace('template.f90',dic)
    lotus.runDocker(n,folder,resume)


# low RE for validation
Re = 500

for resolution in [128]:
    for thick in [1,2,3,4,6]:
    # folder = 'arclow_{:01d}'.format(int(resolution))
        folder = 'ARC_{:01d}_{:01d}'.format(int(resolution),int(thick))
        runcase(folder,{'RES':str(resolution),'DR':str(thick),'RE':str(Re)})
print('Run all python files for postprocessing')
subprocess.call('python proc.py', shell=True)
print('Finished')
