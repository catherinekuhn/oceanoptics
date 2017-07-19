#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 20:18:21 2017

@author: maritescanto
"""

import os, sys, scipy, glob;
from os.path import basename
import numpy as np


#Get list of all L2 files
flist = glob.glob('/Users/maritescanto/Documents/BenthicIrradiance/OceanOptics2017/Data/Lab2/UVA/*')

stnList = np.zeros((len(flist)))    

count = 0
for l2file in flist:
    print(l2file)
    os.path.splitext(basename(l2file))[0]
    