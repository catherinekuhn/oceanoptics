#!/usr/bin/python

import IOP as IOP;

station = 'stn5';
ac9file = station + '/' + station + '.ac9';
ctdfile = station + '/' + station + '.ctd';
bb9file = station + '/' + station + '.bb9';
hs6file = station + '/' + station + '.hs6';
bb9Dev = 'BB9-278.dev';

output = IOP.iopCor(ac9file,ctdfile,bb9file,hs6file,bb9Dev,1,station);