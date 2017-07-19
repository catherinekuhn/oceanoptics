%read in data from ACS easily for all files in a given directory
clc
clear all
close all
%put directory name here
directoryName = '/Users/Akroma/Dropbox (MIT)/oceanoptics2017_students/Lab03/raw_data/ACS_022/';

%Read all file names in the directory
fileNames = dir(directoryName);

%Assign length of file names to variable
lengthFileName = length(fileNames);

% COUNTS loop iterations
COUNT = 1 ;

%Where the first data file occurs 
%get min and max sizes for data storage

for I=4:lengthFileName %Start reading at data files
    
    
   FileName=fileNames(I).name; %switch to different files
   [t,dat,wla,wlc] = rd_wetview_acs_022([directoryName FileName]); % read in data starting from 11t2th row and first column
   
   %create structure
   dat_structure{COUNT} = dat;
   t_structure{COUNT} = t;
   wla_structure{COUNT} = wla;
   wlc_structure{COUNT} = wlc;
   combinedWL_Structure{COUNT} = [wla wlc];
   filename_structure{COUNT} = fileNames(I).name;

  COUNT=COUNT+1;
end 

%Data from 83:164 is a side, data from 1:82 is c side
%scatter(combinedWL_Structure{1},dat_structure{1}(1,:))
plot(wla_structure{3},dat_structure{3}(1,83:164))
title(fileNames(6).name)
