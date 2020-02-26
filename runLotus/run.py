#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import lotus

def runcase(folder,dic,resume=None,n=1):
    print(dic)
    lotus.replace('template.f90',dic)
    lotus.run(n,folder,resume)

if __name__=="__main__":
    lotus.run(2,"out")
