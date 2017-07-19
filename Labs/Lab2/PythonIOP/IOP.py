#===============================================================================
#			---------------------
#				IOP.py
#			---------------------
#
# Description:
# Python module written to process IOP data collected using Wetlabs ac9, bb9,and/or 
# HOBI Labs HS6.  Module also has a function to read in SeaBird Electronics 
# (SBE) 19plus CTD profiler. 
# 
# We are assuming here that the ac9 and bb9 data is being logged to an WET Labs
# data handler (dh4) which time stamps each data line.
# 
# Written: 19 April 2012
# by: Lachlan McKinna, Curtin University
#
# Last modified: 13 July 2015
#
# Disclaimer:
# This code is distributed freely. It may be used and modified for non-commercial 
# purposes (i.e. education and research). This code is provided "as is" WITHOUT 
# ANY EXPRESS OR IMPLIED WARRANTY whatsover. The code is not guaranteed to be fault 
# tolerant. The author accepts no responsibility for potential errors, faults or bugs 
# within this code. The user must assume the entire risk of using this code.
#
#*******************************************************************************
#Import necessary modules
import numpy as np;
import datetime as dt;
import os, sys, csv,shutil,scipy;
import scipy.signal as signal;
import matplotlib.pyplot as plt;
from scipy.optimize.minpack import leastsq;
#
#-------------------------------------------------------------------------------
# Constants
#-------------------------------------------------------------------------------
#Calibration temperature for ac9 meter in degrees celcius.  Taken from the header
#of ac9 file.
Tcal = 19.2;
#-------
# Average pure water offsets collected at AIMS - 28-06/-2011
# Data required for ac9 processing
#
#aOffset = np.array((0.0571285, 0.14513943, 0.12520539, 0.22843816, 0.14116505, \
#			0.1953659, 0.19019437, 0.12314558,  0.18347907),dtype = np.float64)
#
#cOffset = np.array((0.14149699, 0.4120653, 0.29074376, 0.67741285, 0.40469971, \
#			0.40665453, 0.47635536, 0.28058782, 0.34157937),dtype = np.float64)
#
###
# Average pure water offsets collected onboard RV Cape Ferguson - 05-04/-2012
#aOffset = np.array((-0.042453792, 0.1168135, 0.10477414, 0.223477931, 0.140792865,0.187638269, 0.190845965, 0.138743028, 0.205017956))
#
#cOffset = np.array((0.01038963, 0.403365349, 0.29839399, 0.734282851, 0.42554881,0.460359771, 0.52486882, 0.33413454, 0.409931995))
#
#Assume zero offsets - instrument just unpacked (returned) from Wetlabs after factory cal
aOffset = np.zeros((9,));
cOffset = np.zeros((9,));
#-----
#Temperature-salinity correction coefficients for ac9
#
#ac9 wavelength are -> 412.,440.,488.,510.,532.,555.,650.,676., and 712.
#
psiT = np.array((0.0, 0.0, 0.0, 0.0001, 0.0, 0.0, 0.0, -0.0001, 0.0032 ),\
		dtype=np.float64);
#
psiSc = np.array((-0.00002, -0.00004, -0.00004, -0.00004, -0.00004, -0.00004, \
		-0.00001, -0.00004, -0.00023),dtype=np.float64);
#
psiSa = np.array((0.00004, 0.00002, 0.00001,0.00001, 0.00001, 0.0, 0.00003, -0.00001,\
		-0.0002),dtype=np.float64);
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Function: datainterp
#
# Description: Interpolate input data to a new array length
#
def datainterp(xnew,xin,yin):
	#
	ynew = scipy.interp(xnew, xin, yin, left=None, right=None);
	#
	return ynew
#
#-------------------------------------------------------------------------------
#Function: ctdRead
# 
#Description: Reads in calibrated salinity, temperature and depth data from a
#SeaBird Electronics (SBE) 19plus CTD seacat profiler.
#NOTE: File format are .CNV converted files - not raw .HEX files
#
def ctdRead(filename):
	#
	#how many lines in file?
	incr = 0;
	numtotal = 0;
	#
	for line in open(filename):
	#
		rowwidth = len(np.fromstring(line,dtype=np.float64,sep=' '));
		#
		if '*END*' in line:
			startrow = incr + 1;
		incr = incr + 1;
		numtotal = numtotal + 1;
	#
	#
	FID = open(filename);
	lines = FID.readlines();
	lines = lines[startrow:numtotal];
	#
	data = np.zeros((rowwidth,));
	#
	for i in range(0,len(lines)):
		datarow = np.fromstring(lines[i],dtype=np.float64,sep=' ');
		data = np.asarray(np.vstack((data,datarow)),dtype=np.float64);
	#
	FID.close();
	#
	data = np.delete(data,0,0);
	#
	#NOTE: These columns of the 'data' block may vary on the CTD instrument you use 
	#		and what additional instruments you have attached. In this case, we have 
	#		CTD with a wetlabs Chl fluorometer attached. I have commented out the 
	#		the backscattering and trasmissometer instruments here.
	time = data[:,0];
	depth = data[:,1];
	temp = data[:,2];
	sal = data[:,3];
	fl = data[:,4];
	#bb = data[:,5];
	#cstar = data[:,6];
	#
	return time, depth, temp, sal, fl
#
#-------------------------------------------------------------------------------
#Function: ctdDh4Read
#
#Description: Reads CTD data extracted from a WetLabs DH4 logger.
#
def ctdDh4Read(filename):
	#
	#DH4 extracted CTD files have a single header line.
	#Time(ms)	P(dbars)		Temp(C)		Conduct		Sal(PSU)		Spare(V)		Spare(V)		Spare	(V)		Spare(V)
	#
	#How many lines in ctd file?
	numtotal = 0
	#
	for line in open(filename):
		#
		rowwidth = len(np.fromstring(line,dtype=np.float64,sep=' '));
		#
		numtotal = numtotal + 1;
		##
	FID = open(filename);
	lines = FID.readlines();
	lines = lines[1:numtotal];
	data = np.zeros((rowwidth,));
	#
	for i in range(0,len(lines)):
		datarow = np.fromstring(lines[i],dtype=np.float64,sep=' ');
		data = np.asarray(np.vstack((data,datarow)),dtype=np.float64);
		##
	#
	time = data[:,0];
	depth = data[:,1];
	temp = data[:,2];
	sal = data[:,4];
	#
	return time, depth, temp, sal
##
#-------------------------------------------------------------------------------
# Function: ac9read
# Description: Reads WETLabs ac9 data and applies pure water offsets.
#
def ac9read(filename):
	#
	#how many lines in file?
	incr = 0;
	numtotal = 0;
	#
	#Find the row number at the beginning of data (ie after header)
	for line in open(filename):
		#
		if '1		; aquisition bin size' in line:
			startrow = incr + 1;
		#
		incr = incr + 1;
		numtotal = numtotal + 1;
	#
	#
	FID = open(filename);
	lines = FID.readlines();
	#
	#Get lines in the file that occur after the header
	lines = lines[startrow:numtotal];
	#
	#create null array
	data = np.zeros((24,));
	#
	for i in range(0,len(lines)):
		datarow = np.fromstring(lines[i],dtype=np.float64,sep=' ');
		data = np.asarray(np.vstack((data,datarow)),dtype=np.float64);
	#end loop
	#
	FID.close();
	#Remove dummy array at beginning of matrix
	data = np.delete(data,0,0);
	#
	#strip out the desired channels
	time = data[:,0];
	#time = time - time[0]

	depth = data[:,21];
	ac9Temp = data[:,19];
	#
	#wavelengths
	lam = np.array((412,440,488,510,532,555,650,676,715));
	#
	#absorption channels
	a412 = data[:,7] - aOffset[0];
	a440 = data[:,8] - aOffset[1];
	a488 = data[:,9] - aOffset[2];
	a510 = data[:,13] - aOffset[3];
	a532 = data[:,14] - aOffset[4];
	a555 = data[:,15] - aOffset[5];
	a650 = data[:,1] - aOffset[6];
	a676 = data[:,2] - aOffset[7];
	a712 = data[:,3] - aOffset[8];
	#
	#Attenuation channels
	c412 = data[:,16] - cOffset[0];
	c440 = data[:,17] - cOffset[1];
	c488 = data[:,18] - cOffset[2];
	c510 = data[:,4] - cOffset[3];
	c532 = data[:,5] - cOffset[4];
	c555 = data[:,6] - cOffset[5];
	c650 = data[:,10] - cOffset[6];
	c676 = data[:,11] - cOffset[7];
	c712 = data[:,12] - cOffset[8];
	#
	#Smooth data using a median average, moving window filter.  This will 
	#minimise any noisy spikes in the data
	a412= signal.medfilt(a412,15);
	a440= signal.medfilt(a440,15);
	a488= signal.medfilt(a488,15);
	a510= signal.medfilt(a510,15);
	a532= signal.medfilt(a532,15);
	a555= signal.medfilt(a555,15);
	a650= signal.medfilt(a650,15);
	a676= signal.medfilt(a676,15);
	a712= signal.medfilt(a712,15);
	#
	c412= signal.medfilt(c412,15);
	c440= signal.medfilt(c440,15);
	c488= signal.medfilt(c488,15);
	c510= signal.medfilt(c510,15);
	c532= signal.medfilt(c532,15);
	c555= signal.medfilt(c555,15);
	c650= signal.medfilt(c650,15);
	c676= signal.medfilt(c676,15);
	c712= signal.medfilt(c712,15);
	#
	a = np.vstack((a412,a440,a488,a510,a532,a555,a650,a676,a712));
	c = np.vstack((c412,c440,c488,c510,c532,c555,c650,c676,c712));
	x = a.shape;
	#
	return time, depth, lam, a, c
	#
#
#-------------------------------------------------------------------------------
# Function: bb9devRead
#
# Description: Function reads in the device file for a WET Labs bb9 instrument.  
# The dev file contains the calibration coefficients and dark offsets.
#
def bb9devRead(filename):
	#
	count = 0
	#
	for line in open(filename):
	#
		if count == 9:
			dataline = line.split();
			dataline = np.array((dataline));
			dark1 = dataline[2];
			scale1 = dataline[1];
		#
		elif count == 11:
			dataline = line.split();
			dataline = np.array((dataline));
			dark2 = dataline[2];
			scale2 = dataline[1];
		#
		elif count == 13:
			dataline = line.split();
			dataline = np.array((dataline));
			dark3 = dataline[2];
			scale3 = dataline[1];
		#
		elif count == 15:
			dataline = line.split();
			dataline = np.array((dataline));
			dark4 = dataline[2];
			scale4 = dataline[1];
		#
		elif count == 17:
			dataline = line.split();
			dataline = np.array((dataline));
			dark5 = dataline[2];
			scale5 = dataline[1];
			#
		elif count == 19:
			dataline = line.split();
			dataline = np.array((dataline));
			dark6 = dataline[2];
			scale6 = dataline[1];
		#
		elif count == 21:
			dataline = line.split();
			dataline = np.array((dataline));
			dark7 = dataline[2];
			scale7 = dataline[1];
		#
		elif count == 23:
			dataline = line.split();
			dataline = np.array((dataline));
			dark8 = dataline[2];
			scale8 = dataline[1];
		#
		elif count == 25:
			dataline = line.split();
			dataline = np.array((dataline));
			dark9 = dataline[2];
			scale9 = dataline[1];
		#
		count = count +1;
	
		#
	offset = np.array((dark1,dark2,dark3,dark4,dark5,dark6,dark7,dark8,dark9),dtype=np.float64);
	scale = np.array((scale1,scale2,scale3,scale4,scale5,scale6,scale7,scale8,scale9),dtype=np.float64);
	#
	return scale, offset
#
#-------------------------------------------------------------------------------
# Function: bb9read
#
# Description: Reads in WET Labs bb9 data, aplies calibration coefficients, 
# subtracts dark offset.  Returns the calibrated fixed-angle VSF 
#
#
def bb9Read(filename,devfile):
	#
	#Read in calibration data from .dev file
	scale, offset = bb9devRead(devfile);
	#
	#bb9lam = np.array((405.,443.,493.,504.,523.,592.,655.,676.,717.));
	bb9lam = np.array((412.,440.,488.,510.,532.,595.,660.,676.,715.));
	#
	#Create null values
	betaCal = np.zeros((9,));
	#
	#how many lines in file?
	incr = 0;
	numtotal = 0;
	#
	time = 0.;
	#
	for line in open(filename):
		#
		dataline = line.split();
		dataline = np.array((dataline));
		#
		time = np.append(time,dataline[0]);
		#
		beta= np.array((dataline[6],dataline[8],dataline[10],dataline[12],\
				dataline[14], dataline[16], dataline[18], dataline[20],\
				dataline[22]),dtype=np.float64);
		#
		betaNew = scale*(beta - offset);
		#
		betaCal = np.vstack((betaCal,betaNew));
		#
		#
	#
	#Assign the time array as a float64
	time = np.array(time,dtype=np.float64);
	#
	#Remove dummy values from beginning of arrays
	betaCal = np.delete(betaCal,0,0);
	time = np.delete(time,0,0);
	#
	#Subtract the time offset, so that the first time element is 0 ms
	time = time - time[0];
	#
	betaCal = np.transpose(betaCal);
	#
	return time, bb9lam, betaCal
	##
#
#-------------------------------------------------------------------------------
# Function: hs6read 
#
# Description: Function to read HydroScat6 ".dat" file.
#
def hs6read(filename):
	#
	#Hydrscat sensor-specific wavelengths
	#Note these may vary depending on the instrument used
	hs6lam = np.array((420,442,470,510,590,700));
	#
	#how many lines in file?
	incr = 0;
	numtotal = 0;
	#
	for line in open(filename):
		#
		rowwidth = len(np.fromstring(line,dtype=np.float64,sep=','));
		#
		if '[Data]' in line:
			startrow = incr + 1;
		#
		incr = incr + 1;
		numtotal = numtotal + 1;
	#
	#print startrow
	#print rowwidth
	#
	FID = open(filename);
	lines = FID.readlines();
	lines = lines[startrow:numtotal];
	#
	DATA = np.zeros((rowwidth,));
	#
	for i in range(0,len(lines)):
		DATAROW = np.fromstring(lines[i],dtype=np.float64,sep=',');
		DATA = np.asarray(np.vstack((DATA,DATAROW)),dtype=np.float64);
	#end loop
	#
	#Close file
	FID.close();
	#
	#Remove dummy array at beginning of matrix
	DATA = np.delete(DATA,0,0);
	#
	hs6time = DATA[:,0];
	hs6depth = DATA[:,1];
	#
	betabb420uncorr = DATA[:,26];
	betabb510uncorr = DATA[:,27];
	betabb442uncorr = DATA[:,28];
	betabb700uncorr = DATA[:,29];
	betabb470uncorr = DATA[:,30];
	betabb590uncorr = DATA[:,31];
	betafl510uncorr = DATA[:,32];
	betafl700uncorr = DATA[:,33];
	#
	betaUncorr = np.array((betabb420uncorr,betabb442uncorr,betabb470uncorr,\
	betabb510uncorr,betabb590uncorr,betabb700uncorr));
	#
	return hs6time, hs6depth, hs6lam, betaUncorr
#-------------------------------------------------------------------------------
#Function: tsCorMerged
#
#Description: Applied temperature-salinity corrections to ac9 data as described
#in the WET Labs ac-meter protocols.
#
#NOTE:
#In this function, we assume that the t,s,a,c have been depth merged already.
#tsCorMerged(aCrop,cCrop,tempCrop,salCrop)
#
def tsCorMerged(a,c,t,s):
	#
	#
	aCorr = np.zeros((9,));
	cCorr = np.zeros((9,));
	#
	for i in range(0,len(a.T)):
		#
		#Apply temperature-salinity corrections
		aCorNew = a[:,i] - (psiT*(t[i] - Tcal) + psiSa*s[i]);
		cCorNew = c[:,i] - (psiT*(t[i] -Tcal ) + psiSc*s[i]);
		#
		aCorr = np.vstack((aCorr,aCorNew));
		cCorr = np.vstack((cCorr,cCorNew));
	#
	#
	aCorr = np.delete(aCorr,0,0);
	cCorr = np.delete(cCorr,0,0);
	#
	return aCorr, cCorr
#
#-------------------------------------------------------------------------------
#Function betaSwHs6
#
#Descripton: computes back scattering of pure water and the fixed-angle volume 
#				scattering angle for pure wter. The centroid angle varies for
#				the fixed-angle VSF scattering meter used (e.g. bb9, h26)
#
def betaSwHs6(sensorLam,salinity,centroidAngle):
	#Calculate pure water scattering data using Morel functions
	#
	#define conststants in calculation
	S = salinity;
	delta = 0.09;
	theta = np.radians(centroidAngle);
	p90 = (1.0-delta)/(1.0+delta);
	#
	#Wavelength vector
	#
	longLam = np.arange(400,800,1);
	#
	#Morel 1974 model for pure water b and bb from Twardowski et al (2007)
	#bpwFull = ((3.50e-4)*((longLam/550.)**(-4.32)));
	bpwFull = 0.0029308*(longLam/500.)**(-4.24);
	bbpwFull = bpwFull/2.;
	#
	#Morel 1974 model for betaw taken from Twardowski et al (2007)
	#betawFull = (2.18e-4)*((longLam/450.)**(-4.32))* (1+ p90*((np.cos(theta)**2))) \
	#	*((1+0.3*S)/37);
	betawFull = ((1.38)*(longLam/500.0)**(-4.32)) * (1.0 + 0.3*S/37.0)*1e-4 *(1.0 + p90*((np.cos(theta)**2)));
	#
	betaWOut = 0.;
	bbWOut = 0.;
	#
	#
	for jj in range(0,len(longLam,)):
		for kk in range(0,len(sensorLam)):
			if longLam[jj] == sensorLam[kk]:
				#
				betaWOut = np.append(betaWOut,betawFull[jj]);
				bbWOut = np.append(bbWOut,bbpwFull[jj]);
				##
			##
		##
	#
	betaWOut = np.delete(betaWOut,0,0);
	bbWOut = np.delete(bbWOut,0,0);
	#
	return betaWOut,bbWOut
#
#-------------------------------------------------------------------------------
#Function: betaPLcorr
#
#Description: Function corrects for pathlength attenuation - bb9
#
def betaPLcorr(beta,a,sensorlam,sal,depth):
	#
	beta = beta.T
	length,width = beta.shape;
	betaCorr = np.zeros((length,width));
	bbp = np.zeros((length,width));
	bb = np.zeros((length,width));
	#
	#
	for j in range(0,len(a)):
		#
		#
		betaPLcorr = beta[j,:] * np.exp(0.0391*a[j,:]);
		#
		#Get seawater scattering values
		betaSW,bbSW = betaSwHs6(sensorlam,sal[j],117.);
		#
		#PERFORM SEAWATER CORRECTIONS....
		xibb9 = 1.1;
		#
		betaCorr[j,:] = betaPLcorr - betaSW;
		bbp[j,:] = 2.0*np.pi*xibb9*(betaCorr[j,:]);
		bb[j,:] = 2.0*np.pi*xibb9*(betaCorr[j,:]) + bbSW;
		
	#bbp = np.delete(bbp,0,0);
	#bb = np.delete(bb,0,0);
	#betaCorr = np.delete(betaCorr,0,0);
	#
	return bb,bbp, betaCorr
#
#-------------------------------------------------------------------------------
#Function: sigmCorr
#
#Description: Function corrects for pathlength attenuation - Hydroscat
#
def sigmaCorr(beta,a,c,sensorlam,sal,depth):
	#
	beta = beta.T;
	length,width = beta.shape;
	betaCorr = np.zeros((length,width));
	bbp = np.zeros((length,width));
	bb = np.zeros((length,width));
	#
	b = c - a;
	kexp = np.array((0.140,0.143,0.143,0.139,0.143,0.141));
	k1 = np.array((1,1,1,1,1,1));
	#
	for j in range(0,len(a)):
		#
		kbb = a[j,:] + 0.4* b[j,:];
		k = kexp*kbb;
		sigma = k1*np.exp(k);
		#
		#
		#Get seawater scattering values
		betaSW,bbSW = betaSwHs6(sensorlam,sal[j],140.);
		#
		betaCorr[j,:] = beta[j,:] * sigma - betaSW;
		#PERFORM SEAWATER CORRECTIONS....
		xiHs6 = 1.08;
		#
		bb[j,:] = 2.0*np.pi*xiHs6*(betaCorr[j,:]) + bbSW;
		bbp[j,:] = 2.0*np.pi*xiHs6*(betaCorr[j,:]);
	#
	#bbp = np.delete(bbp,0,0);
	#bb = np.delete(bb,0,0)
	#betaCorr = np.delete(betaCorr,0,0);
	#
	#
	return bb,bbp,betaCorr
#
#-------------------------------------------------------------------------------
#Function: ctdInterp
#
#Description: Interpolate the bb9 data to have the same depth dimensions as the ac9 
#				by interpolation and smoothing.
#
def ctdInterp(t,s,depth1,depth2):
	#
	temp = datainterp(depth1,depth2,t);
	sal = datainterp(depth1,depth2,s);
	#
	temp = signal.medfilt(temp,15);
	sal = signal.medfilt(sal,15);
	#
	return temp, sal
#
#-------------------------------------------------------------------------------
#Function: scatCor
#
#Description: Applies a scattering correction to the a-tube data to correct for 
# 				the small proportion of light that is scattered out of the tube and not 
# 				toward the detector (for ac-meter)
#
def scatCor(aCorr,cCorr,method):
	#
	bCorr = cCorr - aCorr;
	#
	#Create null array
	aScatCorr = np.zeros((9,));
	#
	#
	#----------
	if method == 1:
		#Apply the NIR offset correction
		#
		for i in range(0,len(aCorr)):
			#	
			#line= np.array((aCorr[i]),dtype=np.float64)
			#offset = line[8]
			#
			aScatterCorr = aCorr[i,:] - aCorr[i,8];
			#
			aScatCorr = np.vstack((aScatCorr,aScatterCorr));
			##
		##
	##
	#
	#---------
	if method == 2:
		#Apply the fixed proportion of b scattering correction
		#Set proportion to 0.3
		#
		for i in range(0,len(aCorr)):
			#
			#Case 1 - predominantly biological partilcs
			#aScatterCorr = aCorr[i] - 0.14*bCorr[i]
			#
			#Case 2 - predominantly suspended particles
			#!! Probably better for inshore Great Barrier Reef !!
			aScatterCorr = aCorr[i] - 0.14*bCorr[i];
			#
			aScatCorr = np.vstack((aScatCorr,aScatterCorr));
			#
	#---------
	if method == 3:
		#Apply a combination of methods 1 and 2
		#
		for i in range(0,len(aCorr)):
			#
			#Note the offset shold be between 0.07 and 0.35	
			aOffset = (aCorr[i,8]/bCorr[i,8]) * (bCorr[i]);
			#
			aScatterCorr =  aCorr[i] - aOffset;
			#
			aScatCorr = np.vstack((aScatCorr,aScatterCorr));
	#
	aScatCorr = np.delete(aScatCorr,0,0);
	#
	return aScatCorr
#	
#-------------------------------------------------------------------------------
#Function: bbpfit
#
#Description: fits Bbp using a power law using Levenberg-Marquardt least squares
#fitting routine. Fit bbp sensor wavelengths and ac-meter wavelengths
#
def bbpfit(p,sensorLam,dataRef,lamRef):
	fittedBbp = dataRef * (lamRef/sensorLam)**p[0];
	return fittedBbp
#
def residuals(p, sensorLam, bbpdata,dataRef,lamRef):
	bbpFitted = bbpfit(p,sensorLam,dataRef,lamRef);
	err = bbpdata - bbpFitted;
	return err
#
def bbpfitting(sensorLam,data,sensor,ac9lam):
	#
	if sensor == 'hs6':
		refLam = 590.0;
	elif sensor == 'bb9':
		refLam = 532.0;
	elif sensor == 'ac9':
		refLam = 555.0;
	#
	#print sensorLam
	lamDiff = sensorLam - refLam;
	idx = np.where(np.fabs(lamDiff) == np.amin(np.fabs(lamDiff)))[0][0];
	#
	#lamnew = np.delete(bb9lam,0,0)
	#datanew = np.delete(data,0,0)
	#
	guess = np.array((0.5));
	#
	pB, succes = leastsq(residuals, guess, args=(sensorLam,data,data[idx],refLam),\
		maxfev=1000, factor= 0.1,xtol=1.49012e-09);
		#
	#
	#Fit bb9 data to the spectral resolution of ac9
	#ac9lam = np.array((412.,440.,488.,510.,532.,555.,650.,676.,712.),dtype=np.float64)
	#
	#Fitted at bbp sensor wavelengths
	fittedBB = data[idx] * (refLam/sensorLam)**pB[0];
	#
	#fitted at ac9 sensor wavelengths
	fittedBBac9 =  data[idx] * (refLam/ac9lam)**pB[0];
	#
	#print pB
	return fittedBB,pB,fittedBBac9
#
#-------------------------------------------------------------------------------
#Function: depthBinFunc
#
#Description: Perform depth binning for a, c, beta9, beta6, t,s and fl channels.
#				Binsize is typically  set to 1 (for 1 metre depth bins).
#
def depthBinFunc(dataIn,depthIn,binsize):
	#
	maxDepth = np.amax(depthIn);
	#After the zeroth-depth, the increments in binArray represent 
	#the bin centres.
	binArray = np.arange(2.0,maxDepth,binsize);
	#print binArray
	#
	#Surface bin (zeroth bin) is from 0 - 0.3 m
	surfaceBinIndex = 0;
	#
	for k in range(0,len(dataIn.T)):
		#
	#	print k
		if depthIn[k] >= 0.0 and depthIn[k] <= 2.0:
			surfaceBinIndex = np.append(surfaceBinIndex,k);
			##
		##
	##
	surfaceBinIndex = np.delete(surfaceBinIndex,0,0);
	#print surfaceBinIndex
	#print depthIn
	data = dataIn[:, surfaceBinIndex];
	#
	#Create null bin arrays to be filled
	dataBin = np.zeros((len(binArray),));
	depthBin = np.zeros((len(binArray),));
	wid,length = dataIn.shape;
	#
	#
	dataBin = np.zeros(((len(binArray)),wid));
	#
	#
	#Values for surface bin determined by crude median average
	#Test if there is only a single set of surface data
	#How many data points in surface bin?
	dims =  dataIn[:, surfaceBinIndex].ndim;
	#If the number of samples in the surface bin is only 1
	if dims == 1:
		dataBin[0,:] = dataIn[:, surfaceBinIndex];
	#If the number of samples in the surface bin is more than 1
	else:
		dataBin[0,:] = np.mean(dataIn[:, surfaceBinIndex],axis=1);
	#
	#
	#
	for i in range(1,len(binArray)):
		#
		binCenter = binArray[i];
		binHalfUpper = binCenter - binsize/2.;
		binHalfLower = binCenter + binsize/2.;
		#
		binIndex = 0
		#
		for j in range(0,len(depthIn)):
		#
			if depthIn[j] >= binHalfUpper and depthIn[j] <= binHalfLower:
				binIndex = np.append(binIndex,j);
		#
		binIndex = np.delete(binIndex,0,0);
		#
		#DepthBinValues
		#print dataIn[:, binIndex]
		dataDep = np.mean(dataIn[:, binIndex],axis=1);
		depthDep = np.mean(depthIn[binIndex]);
		#
		dataBinned = (dataDep - dataBin[i-1,:])*(float(i) - depthBin[i-1]) / (depthDep - depthBin[i-1]) + dataBin[i-1,:]
		dataBin[i,:] = dataBinned;
		#
		#print dataBin
		#print depthDep
		#
	dataBin = dataBin.T;
	depthBin = binArray;
	#
	return depthBin, dataBin;
	#
##
#-------------------------------------------------------------------------------
#				HYDROLIGHT FILE WRITERS
#
#Below files are written out in Hydrolight input formats. This may help in optical 
#closure studies.
#
#Write Chlorophyll data to Hydrolight format file
def chlHLout(outChl,commonDepth,station):
	#
	print outChl
	print commonDepth
	#
	file_out = station + '_chl.dat';
	outChl = np.vstack((commonDepth,outChl));
	
	arrayin = outChl.T;
	#
	ofile = open(file_out,'wb');
	filewriter = csv.writer(ofile, delimiter='\t');
	#
	ofile.write(file_out + " Standard Hydrolight Input\n");
	ofile.write("Chloropyll data collected from Heron Island " + station+"\n");
	ofile.write("Input format for Hydrolight/Ecolight - 10 header rows \n");
	ofile.write("a \n");
	ofile.write("b \n");
	ofile.write("c \n");
	ofile.write("d \n");
	ofile.write("d \n");
	ofile.write("3 \n");
	ofile.write("depth  (m)   chl (ug/L) \n");
	filewriter.writerows(arrayin);
	ofile.write("-1.0  -1.0");
	#
	ofile.close();
	#
	return
#
#Write ac9 data to Hydrolight format file
def ac9HLout(ac9data,commonDepth,station):
	#
	file_out = station + '_ac9.dat';
	outAC9 = np.vstack((commonDepth,ac9data));
	
	arrayin = outAC9.T;
	#
	ofile = open(file_out,'wb')
	filewriter = csv.writer(ofile, delimiter='\t')
	#
	ofile.write(file_out + " Standard Hydrolight Input\n");
	ofile.write("ac9 data collected from Heron Island " + station+"\n");
	ofile.write("Input format for Hydrolight/Ecolight - 10 header rows \n");
	ofile.write("a \n");
	ofile.write("b \n");
	ofile.write("c \n");
	ofile.write("d \n");
	ofile.write("the depths shown are the average depth of the data in each bin.  Column headers (depth in m; and and c in 1/m):    \n")
	ofile.write("   depth    a412        a440        a488        a510        a532        a555        a650        a676        a715        a412        c440        c488        c510        c532        c555        c650        c676        c715 \n")
	ofile.write("Record 11 gives the number of wavelengths and the wavelengths \n")
	ofile.write("   9  412.0  440.0  488.0  510.0  532.0  555.0  650.0  676.0  715.0 \n")
	filewriter.writerows(arrayin)
	ofile.write("  -1.000   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00");
	#
	ofile.close()
	#
	return
	#
#
#Write HS6 data to Hydrolight format file
def hs6HLout(hs6data,commonDepth,station):
	#
	file_out = station + '_hs6.dat';
	outHS6 = np.vstack((commonDepth,hs6data));
	
	arrayin = outHS6.T;
	#
	ofile = open(file_out,'wb');
	filewriter = csv.writer(ofile, delimiter='\t');
	#
	ofile.write(file_out + " Standard Hydrolight Input\n");
	ofile.write("Hydroscat-6 data collected from Heron Island " + station+"\n");
	ofile.write("Input format for Hydrolight/Ecolight - 10 header rows \n");
	ofile.write("a \n");
	ofile.write("b \n");
	ofile.write("c \n");
	ofile.write("d \n");
	ofile.write("Column headers (depth in m and bb in 1/m): \n");
	ofile.write(" depth      bb420       bb442       bb470       bb510       bb590       bb620       bb700 \n");
	ofile.write("Record 11 gives the number of wavelengths and the wavelengths: \n");
	ofile.write("6 420.0  442.0  470  510.0  620.0  700.0 \n");
	filewriter.writerows(arrayin);
	ofile.write("  -1.000  0.0         0.0         0.0         0.0         0.0         0.0");
	#
	ofile.close();
	#
	return
#
#-------------------------------------------------------------------------------
#Function: ac9bb9hs6Cor
#
#Description: this function does the "heavy lifting" and calls most of the subroutines
3
def ac9bb9hs6Cor(ac9filename,bb9filename,hs6filename,ctdfilename,bb9Devfilename,binsize,station):
	#
	print hs6filename
	#
	#read in ac9 data
	time1, depth1, ac9lam, a, c  = ac9read(ac9filename)
	#
	#Read in bb9 data
	#bb9Devfilename = 'BB9-278.dev';
	time4, bb9lam, betaCal = bb9Read(bb9filename,bb9Devfilename);
	#Interpolate to get depth of bb9
	dummyVector1 = np.linspace(0,time1[len(time1)-1], num=len(time1))
	dummyVector2 = np.linspace(0,time1[len(time1)-1], num=len(time4))
	depth4 = datainterp(dummyVector2,dummyVector1,depth1);
	
	#
	#read in ctd data
	#time2, depth2, temp, sal, fl, bb, cstar = ctdRead(ctdfilename)
	time2, depth2, temp, sal,fl = ctdRead(ctdfilename)
	#
	#read hs6 data (calibrated fixed-angle volume scattering functions)
	hs6time, depth3, hs6lam, betabbUncorr = hs6read(hs6filename)
	#
	#Interpolate and smooth the ctd data to have the same dimensions as the ac9 data
	#temp,sal = ctdInterp(temp,sal,time1,time2)
	#
	#
	#Crop the ac9 data
	depth1 = signal.medfilt(depth1,201)
	depthMax1 = np.where(depth1 == np.amax(depth1));
	maxPoint1 = depthMax1[0][0];
	diffDepth1 = np.diff(depth1,2);
	dropPoint1 = np.where(diffDepth1 > 0.05)[0][0]  + 20;
	depthCrop1 = depth1[dropPoint1:maxPoint1];
	aCrop = a[:,dropPoint1:maxPoint1];
	cCrop = c[:,dropPoint1:maxPoint1];
	time1Crop = time1[dropPoint1:maxPoint1];
	#
	#----
	#Crop data so we only consider the down-cast
	#
	#Crop the ctd data
	#immediate Crop - allow for data logged when deploying at the surface
	depth2 = signal.medfilt(depth2,201)
	depthMax2 = np.where(depth2 == np.amax(depth2));
	maxPoint2 = depthMax2[0][0];
	diffDepth2 = np.diff(depth2,2);
	dropPoint2 = np.where(diffDepth2 > 0.02)[0][0] + 150;
	depthCrop2 = depth2[dropPoint2:maxPoint2];
	tempCrop = temp[dropPoint2:maxPoint2];
	salCrop = sal[dropPoint2:maxPoint2];
	flCrop = fl[dropPoint2:maxPoint2];
	#
	#Crop the hs6 data
	depth3 = signal.medfilt(depth3,3)
	depthMax3 = np.where(depth3 == np.amax(depth3));
	maxPoint3 = depthMax3[0][0];
	diffDepth3 = np.diff(depth3,2);
	dropPoint3 = np.where(diffDepth3 > 0.05)[0][0];
	depthCrop3 = depth3[dropPoint3:maxPoint3];
	hs6betaCrop = betabbUncorr[:,dropPoint3:maxPoint3];
	#
	#
	#crop the bb9 dataBin
	depth4 = signal.medfilt(depth4,7);
	depthMax4 = np.where(depth4 == np.amax(depth4));
	maxPoint4 = depthMax4[0][0];
	diffDepth4 = np.diff(depth4,2);
	dropPoint4 = np.where(diffDepth4 > 0.05)[0][0] ;
	depthCrop4 = depth4[dropPoint4:maxPoint4];
	bb9betaCrop = betaCal[:,dropPoint4:maxPoint4];
	time4Crop = time4[dropPoint4:maxPoint4];
	#
	#
	#add an offset from the pressure sensor to the sensors
	ac9Depth = depthCrop1 + 0.8;
	ctdDepth = depthCrop2;
	hs6Depth = depthCrop3 + 0.3;
	bb9Depth = depthCrop4 + 0.8;
	#
	#
	del depthCrop1,depthCrop2,depthCrop3, depthCrop4;
	#
	ctfCrop= np.vstack((salCrop,tempCrop))
	ctfCrop= np.vstack((ctfCrop,flCrop))
	#
	#--------------
	#Perform depth binning
	depthAc9Bin, aBin  = depthBinFunc(aCrop,ac9Depth,binsize);
	depthAc9Bin, cBin  = depthBinFunc(cCrop,ac9Depth,binsize);
	depthCtfBin, ctfBin  = depthBinFunc(ctfCrop,ctdDepth,binsize);
	depthBb9Bin, betaBb9Bin  = depthBinFunc(bb9betaCrop,bb9Depth,binsize);
	depthHs6Bin, betaHs6Bin  = depthBinFunc(hs6betaCrop,hs6Depth,binsize);
	#
	#---------------
	#Do all instruments reach the same depth?
	#If not, then go to the shallowest one.
	#
	maxDepthAc9 = np.amax(depthAc9Bin);
	maxDepthCtd = np.amax(depthCtfBin);
	maxDepthHs6 = np.amax(depthHs6Bin);
	maxDepthBb9 = np.amax(depthBb9Bin);
	#
	depthArrayAll = np.array((maxDepthAc9,maxDepthCtd,maxDepthBb9,maxDepthHs6));
	#depthArrayAll = np.array((maxDepthAc9,maxDepthCtd,maxDepthBb9));
	depthLimit = np.where(depthArrayAll == np.min(depthArrayAll))[0][0];
	depthLimit = depthArrayAll[depthLimit];
	#depthLimit = np.min(depthArrayAll);
	#print depthLimit
	#
	#Crop all datafile accordinly
	commonDepth = np.arange(2,depthLimit+1,binsize);
	#
	aBin = aBin[:,0:(depthLimit-1)];
	cBin = cBin[:,0:(depthLimit-1)];
	betaBb9Bin = betaBb9Bin[:,0:(depthLimit-1)];
	betaHs6Bin = betaHs6Bin[:,0:(depthLimit-1)];
	#
	salBin = ctfBin[0,0:(depthLimit-1)];
	tempBin = ctfBin[1,0:(depthLimit-1)];
	flBin = ctfBin[2,0:(depthLimit-1)];
	#
	#---------------
	#Call function to apply the temperature salinity corrections
	aCorr, cCorr = tsCorMerged(aBin,cBin,tempBin,salBin);
	#
	#---------------
	#Call function to apply scattering correction to a-tube data
	aScatCorr = scatCor(aCorr,cCorr,1);
	cBin = cBin.T;
	#
	#-------------
	#Interpolate ac9 absorption data to bb9 wavelength
	aBB9interp = np.zeros((len(aScatCorr),9));
	for i in range(0,len(aScatCorr)):
		aBB9interp[i,:] = datainterp(bb9lam,ac9lam,aScatCorr[i,:]);
		##
	#-------------
	#Interpolate ac9 absorption and attenuation to hs6 wavelengths
	cHS6interp = np.zeros((len(aScatCorr),6));
	aHS6interp = np.zeros((len(aScatCorr),6));
	for i in range(0,len(aScatCorr)):
		aHS6interp[i,:] = datainterp(hs6lam,ac9lam,aScatCorr[i,:]);
		cHS6interp[i,:] = datainterp(hs6lam,ac9lam,cCorr[i,:]);
		##
	#
	#------------
	#Apply pathlength attenuation correction to bb9data
	bbBB9, bbpBB9, betaBB9cor = betaPLcorr(betaBb9Bin,aBB9interp,bb9lam,salBin,commonDepth);
	#
	#------------
	#Apply pathlength attenuation correction for hs6 dat
	bbHS6,bbpHS6,betaHS6cor = sigmaCorr(betaHs6Bin,aHS6interp,cHS6interp,hs6lam,salBin,commonDepth);
	#
	#-------------
	#Fit the particulate HS6 bbp data to get spectral slopes
	fittedBBPhs6 = np.zeros((6,));
	fittedBBPac9 = np.zeros((9,));
	gammaBBPhs6 = 0.;
	#
	
	for i in range(0,len(bbHS6)):
		bbpHS6Fit,gammaHS6,bbpAC9Fit= bbpfitting(hs6lam,bbpHS6[i,:],'hs6',ac9lam);
		#
		fittedBBPhs6= np.vstack((fittedBBPhs6,bbpHS6Fit));
		gammaBBPhs6 = np.append(gammaBBPhs6,gammaHS6);
		fittedBBPac9= np.vstack((fittedBBPac9,bbpAC9Fit));
		#
	#
	gammaBBPhs6 = np.delete(gammaBBPhs6,0,0);
	fittedBBPhs6= np.delete(fittedBBPhs6,0,0);
	fittedBBPac9= np.delete(fittedBBPac9,0,0);
	#print gammaBBPhs6;
	#
	#-------------
	#Fit the total HS6 bb data to get spectral slopes
	fittedBBhs6 = np.zeros((6,));
	fittedBBac9 = np.zeros((9,));
	gammaBBhs6 = 0.;
	#
	#
	for i in range(0,len(bbHS6)):
		bbHS6Fit,gammaHS6,bbAC9Fit= bbpfitting(hs6lam,bbHS6[i,:],'hs6',ac9lam);
		#
		fittedBBhs6= np.vstack((fittedBBhs6,bbHS6Fit));
		gammaBBhs6 = np.append(gammaBBhs6,gammaHS6);
		fittedBBac9= np.vstack((fittedBBac9,bbAC9Fit));
		#
	#
	gammaBBhs6 = np.delete(gammaBBhs6,0,0);
	fittedBBhs6= np.delete(fittedBBhs6,0,0);
	fittedBBac9= np.delete(fittedBBac9,0,0);
	#print gammaBBhs6;
	
	##----------------
	#fit the attenuation (c) curve of the ac9
	slopeC = 0.;
	#
	for i in range(0,len(cBin)):
		cFit1,slopeCfit, cFit2= bbpfitting(ac9lam,cBin[i,:],'ac9',ac9lam);
		slopeC = np.append(slopeC,slopeCfit);
		##
	#
	slopeC = np.delete(slopeC,0,0);
	#-----------------
	#Fit the bb9 bbp data to get spectral slopes
	fittedBBPbb9 = np.zeros((9,));
	fittedBBPac92 = np.zeros((9,));
	gammaBBPbb9 = 0.;
	#
	#
	for i in range(0,len(bbBB9)):
		bbpBB9Fit,gammaBB9,bbpAC9Fit2= bbpfitting(bb9lam,bbpBB9[i,:],'bb9',ac9lam);
		#
		fittedBBPbb9= np.vstack((fittedBBPbb9,bbpBB9Fit));
		gammaBBPbb9 = np.append(gammaBBPbb9,gammaBB9);
		fittedBBPac92= np.vstack((fittedBBPac92,bbpAC9Fit2));
		#
	#
	gammaBBPbb9 = np.delete(gammaBBPbb9,0,0);
	fittedBBPbb9= np.delete(fittedBBPbb9,0,0);
	fittedBBPac92= np.delete(fittedBBPac92,0,0);
	#
	#-------------
	#Fit the HS6 data to the same wavelengths as the ac9
	bBin = cBin - aBin.T
	bbRatio =  fittedBBPac9 / bBin
	bbRatio2 = fittedBBPac92 / bBin
	#
	#-------------------
	#fit the attenuation (c) curve of the ac9
	slopeB = 0.;
	#
	for i in range(0,len(bBin)):
		bFit1,slopeBfit, bFit2= bbpfitting(ac9lam,bBin[i,:],'ac9',ac9lam);
		slopeB = np.append(slopeB,slopeBfit);
		#
	
	slopeB = np.delete(slopeB,0,0);
	#
	#-------------------
	#calculate IOP stats
	#-------------------
	#CHL
	avChl = np.median(flBin);
	stdChl = np.std(flBin);
	#
	#BBPGAMMA
	avBBPY = np.median(gammaBBPhs6)
	stdBBPY = np.std(gammaBBPhs6)
	#
	#BBGAMM
	avBBY = np.median(gammaBBhs6)
	stdBBY = np.std(gammaBBhs6)
	#
	#BBP550 
	avBBP550 = np.median(fittedBBPac9[:,5]);
	stdBBP550 = np.std(fittedBBPac9[:,5]);
	#
	#B550
	avB550 = np.median(bBin[:,5]);
	stdB550 = np.std(bBin[:,5]);
	#
	#BSLOPE
	avBSLOPE = np.median(slopeB);
	stdBSLOPE = np.std(slopeB);
	#
	#B550RATIO
	avBBRATIO555 = np.median(bbRatio[:,5]);
	stdBBRATIO555 = np.std(bbRatio[:,5]);
	#
	#A440
	avA440 = np.median(aBin[:,1]);
	stdA440 = np.std(aBin[:,1]);
	#
	#C550
	avC550 = np.median(cBin[:,5]);
	stdC550 = np.std(cBin[:,5]);
	#
	#Slope C
	avCSLOPE = np.median(slopeC);
	stdCSLOPE = np.std(slopeC);
	#
	print ' '
	print '	.......IOP SUMMARY: ' + station + '.......'
	print 'a440 = ' + str(avA440) + ' +/- ' + str(stdA440) + ' m^-1';
	print ' '
	print 'b555 = '+  str(avB550) + ' +/- ' + str(stdB550) + ' m^-1';
	print 'b slope = ' + str(avBSLOPE) + ' +/- ' + str(stdBSLOPE) + ' nm^-1';
	print ' '
	print 'c555 = ' + str(avC550) + ' +/- ' + str(stdC550) + ' m^-1' ;
	print 'c slope = ' + str(avCSLOPE) + ' +/- ' + str(stdCSLOPE) + ' nm^-1';
	print ' ';
	print 'bbp555 = ' + str(avBBP550) + ' +/- ' + str(stdBBP550) + ' m^-1';
	print  'bbp slope = ' + str(avBBPY) + ' +/- ' + str(stdBBPY) + ' nm^-1';
	print 'bb/b 555 = '+ str(avBBRATIO555) + ' +/- ' + str(stdBBRATIO555);
	print ' ';
	print 'chl = ' + str(avChl) + ' +/- ' + str(stdChl) + ' ug/L';
	print '	...........................'
	##
	file_out = station + '_summary.dat';
	ofile = open(file_out,'wb');
	filewriter = csv.writer(ofile, delimiter='\t');
	#
	ofile.write(file_out + " Summary of IOPs " + station + "\n");
	ofile.write("parameter    "  + 'median    ' + 'stdev    '+ "\n");
	ofile.write("chl    " + str(avChl) + '    ' + str(stdChl) + "\n");
	ofile.write("a440    "  + str(avA440) + '    ' + str(stdA440) + "\n");
	ofile.write("c550    " + str(avC550) + '    ' + str(stdC550) + "\n");
	ofile.write("b550    " + str(avB550) + '    ' + str(stdB550) + "\n");
	ofile.write("bbp550    " + str(avBBP550) + '    ' + str(stdBBP550) + "\n");
	ofile.write("slope_b    " + str(avBSLOPE) + '    ' + str(stdBSLOPE)  + "\n");
	ofile.write("slope_c    " + str(avB550) + '    ' + str(stdB550) + "\n");
	ofile.write("slope_bbp    " + str(avBBPY) + '    ' + str(stdBBPY)  + "\n");
	ofile.write("bb555/b555    " + str(avBBRATIO555) + '    ' + str(stdBBRATIO555)  + "\n");
	#
	ofile.close();
	##---------------------
	#save outputs to Hydrolight compatible formats
	bbHS6 = bbHS6.T;	
	ac9data = np.vstack((aBin,cBin.T));
	
	chlHLout(flBin,commonDepth,station);
	ac9HLout(ac9data,commonDepth,station);
	hs6HLout(bbHS6,commonDepth,station);
	
	#---------------------
	#Plot data
	
	plt.figure(1)
	plt.plot(bb9lam,bbpBB9.T,'k:')
	plt.plot(hs6lam,bbpHS6.T,'b:')
	plt.plot(hs6lam,fittedBBPhs6.T,'-bo')
	plt.plot(ac9lam,fittedBBPbb9.T,'-ko')
	
	#plt.plot(bb9lam,bbpBB9.T,'-ko')
	#plt.plot(ac9lam,fittedBBPac92.T,'-ko')
	plt.xlabel('wavelength (nm)')
	plt.ylabel('particulate backscatter (m$^{-1}$)')
	plt.title('b$_{bp}$ ' + station)
	
	plt.figure(2)
	plt.plot(hs6lam,bbHS6)
	plt.plot(hs6lam,fittedBBhs6.T,'bo')
	plt.plot(ac9lam,fittedBBac9.T,'-ro')
	
	plt.plot(bb9lam,bbpBB9.T)
	plt.plot(ac9lam,fittedBBPac92.T,'ko')
	plt.xlabel('wavelength (nm)')
	plt.ylabel('total backscatter (m$^{-1}$)')
	plt.title('b$_{b}$ ' + station)
	
	plt.figure(3)
	plt.plot(ac9lam,bbRatio.T,'r');
	plt.plot(ac9lam,bbRatio2.T,'k');
	plt.xlabel('wavelength (nm)')
	plt.ylabel('particulate backscatter ratio')
	plt.title('b$_{bp}$/b$_{p} $' + station + ' (red is hs6/bb9; black is bb9/ac9) ' )
	
	plt.figure(4)
	plt.plot(gammaBBhs6,-commonDepth,'-ok',label='bb9')
	plt.plot(gammaBBPhs6,-commonDepth,'-or',label='hs6')
	plt.legend();
	plt.ylabel('depth (m)')
	plt.xlabel('particulate backscatter coefficient slope $\gamma_{bbp}$')
	plt.title('$\gamma_{b_{bp}}$ ' + station)
	
	plt.figure(5)
	plt.plot(ac9lam,aCorr.T,'-o')
	plt.xlabel('wavelength (nm)')
	plt.ylabel('absorption (m$^{-1}$)')
	plt.title('a ' + station)
	
	plt.figure(6)
	plt.plot(ac9lam,cBin.T,'-o')
	plt.xlabel('wavelength (nm)')
	plt.ylabel('attenuation (m$^{-1}$)')
	plt.title('c ' + station)
	
	plt.figure(7)
	plt.plot(ac9lam,bBin.T ,'-o')
	plt.xlabel('wavelength (nm)')
	plt.ylabel('scattering (m$^{-1}$)')
	plt.title('b ' + station)
	
	plt.figure(8)
	plt.plot(flBin,-commonDepth,'g-o')
	plt.xlabel('chlorophyll-a concentration ($\mu gL^{-1}$)')
	plt.ylabel('depth (m)')
	plt.title('chl ' + station)
	plt.show()
	
	output = 1.0
	return output
#
#-------------------------------------------------------------------------------
# Function: iopCor
#
# Descritption: function depth bins bb9, ac9 and ctd data.
# Applies temperature-salinity corrections to ac9 data, and pathlength
# attenuation correction for bb9 instrument
#
def iopCor(ac9filename,ctdfilename,bb9filename,hs6filename,bb9DevFilename, binsize,station):
	#
	#
	output = ac9bb9hs6Cor(ac9filename,bb9filename,hs6filename,ctdfilename,bb9DevFilename,binsize,station);
	#
	return output
	#
#-------------------------------------------------------------------------------
'''
station = 'stn5';
ac9file = station + '/' + station + '.ac9';
ctdfile = station + '/' + station + '.ctd';
bb9file = station + '/' + station + '.bb9';
hs6file = station + '/' + station + '.hs6';
bb9Dev = 'BB9-278.dev';

output = iopCor(ac9file,ctdfile,bb9file,hs6file,bb9Dev ,1,station);
#
'''
