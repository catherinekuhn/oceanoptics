Python script: 
IOP.py

Description: 
This code contains a series of functions for reading and processing oceanographic 
instrument data. Specifially, the code was developed to read and depth merge 
IOP instrument data (ac-9, bb-9 and HS-6) with SBE 19 CTD data). Corrections for
effects including: temperature-salinity effects, pathlength attenuation, and 
ac-meter tube scattering. After corrections are applied, the data are depth binned.
Finally, functions exist to write the outputs out in Hydrolight/Ecolight file format
for the purposes of radiative transfer modelling (closure studies).


Background:
The code was written for data collected from a small boat (4.5 m length) where tethering 
instruments into a single profile cage was impracticable. Thus, merging of data using 
an instrument logger (e.g. WetLabs dh4) was only performed for ac-9 and bb9 instruments.
Coincicdent ata collected using a HS6 and SBE 19plus were logged separated and are thus 
recorded in separated files which need merging/binning.


Pre-requisites:
This code was developed using Python 2.7

The following python modules are required by IOP.py:
numpy, datetime, os, sys, csv, shutil, scipy, matplotlib;


Usage:
It is unlikely that IOP.py will be completely compatible with end-user needs. 
Nontheless, basic read/write functions are included and functions used to peform data 
corrections and depth binning are included. Functions can be called from your own 
python script by by importing IOP.py as a module within your code.


Example/test:
I have included a set of test files for processing in the directory "stn 5".
To run the test you can call "main.py" from your terminal using the following command:

"python main.py"

This should process a test station (stn5), produce a series of plots and output
Hydrolight-compatible files. Note, the quality of the bb9 data collected is highly 
questionable.


Disclaimer:
This code is distributed freely. It may be used and modified for non-commercial 
purposes (i.e. education and research). This code is provided "as is" WITHOUT 
ANY EXPRESS OR IMPLIED WARRANTY whatsover. The code is not guaranteed to be fault 
tolerant. The author accepts no responsibility for potential errors, faults or bugs 
within this code. The user must assume the entire risk of using this code.

Lachlan McKinna (Curtin University).
IOP.py was developed as part of an Australian Research Coucil-funded project.





