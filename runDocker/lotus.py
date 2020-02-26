#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess
import os
import shutil


def runDocker(n_proc=0,run_folder='test',read_folder=None):
    "setup and run Lotus using lotus.f90 and the files in postproc using the Docker implementation"

    print('Number of proccessors :{}'.format(n_proc))
    print('Running               :{}'.format(run_folder))

    if read_folder:
        print('Read folder           :{}'.format(read_folder))
    else:
        print('No read folder')

    print('Running in Docker ')
    if read_folder:
        subprocess.call('runDocker '+str(n_proc)+' '+str(run_folder)+' '+str(read_folder), shell=True)
    else:
        subprocess.call('runDocker '+str(n_proc)+' '+str(run_folder), shell=True)


def replace(template,dic):
    """
    Write lotus.f90 file by replacing the dic on the template
    and return a potential folder name
    """
    f1 = open(template,'r')
    f2 = open('lotus.f90','w')
    for text in f1:
        for i, j in dic.items():
            text = text.replace(i, j)
        f2.write(text)
    f1.close()
    f2.close()
    return '_'.join([i+' '+j for i, j in dic.items()])
