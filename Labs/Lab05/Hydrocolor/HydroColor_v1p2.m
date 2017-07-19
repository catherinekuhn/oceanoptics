function varargout = HydroColor_v1p2(filename)
%
%Imports data from the HydroColor v1.2 data file (images optional)
%
%   The HydroColor data file is a text file that contains information about
%   each measurement. The data file, along with corresponding images, can
%   be downloaded from your apple device via iTunes.
%
%   [HCData,HCData_Headers] = HydroColor_v1p2('HydroColor_Datafile.txt') 
%   will import the core HydroColor data (Date, Lat, Lon, RGB Rrs,   
%   turbidity, SPM, bb_red) as a single matrix (HCData). The headers for   
%   the data will be contained in the cell array HCData_Headers. The user   
%   entered name for each measurement will appear in the first column of 
%   HCData_Headers.
%
%   [~,~,HCData_un] = HydroColor_v1p2('HydroColor_Datafile.txt') will 
%   provide a matrix the same size as HCData containing the uncertainty 
%   values for each parameter. The uncertainty units are the same as the 
%   variable except for latitude and longitude, where the units are meters.
%
%   [~,~,~,QC,QC_Headers] = HydroColor_v1p2('HydroColor_Datafile.txt') will
%   output a matrix of data that can be used for quality control. Headers
%   for the data are contained in QC_Headers.
%
%   [~,~,~,~,~,Images] = HydroColor_v1p2('HydroColor_Datafile.txt') will 
%   import the images saved with each measurement. The images must be 
%   located in the same directory as the text file. The three images from  
%   each measurement are stored in a structure where each row corresponds   
%   to the same row in the matrices HCData and QC.
%
%   FUNCTION WITH ALL OUPUTS:
%
%   [HCData,HCData_Headers,HCData_un,QC,QC_Headers,Images] = ...
%   HydroColor_v1p2('HydroColor_Datafile.txt');
%
%
%   Created by THOMAS LEEUW
%   thomas.leeuw@umit.maine.edu
%   May 23, 2014
%
%
%¸.·´¯`·.¸><(((º>.·´¯`·.´¯`·.¸¸.·´¯`·.¸><(((º>.·´¯`·.´¯`·.¸¸.·´¯`·.¸><(((º>


% Import text file headers
fid = fopen(filename);
C = textscan(fid,'%s');
Headers = C{1}(1:29)';

% Import text file data
frewind(fid)
C = textscan(fid,['%s%s%s' repmat('%f',1,20) '%s%s%s%s%s%s'], ...
    'headerlines',1);
% Convert date to matlab datenum
dt = datenum(strcat(C{1},C{2}),'mmddyyyyHHMMSS');
numrows = length(dt);

% Split up rrs and rrs uncertainty, convert string to number
rrs = [];
rrs_un = [];
for ii = 24:26
    tmp = char(C{ii});
    rrs = [rrs str2num(tmp(:,1:5))];
    rrs_un = [rrs_un str2num(tmp(:,8:12))];
end

% Do the same for the rest of the data products
turbidity = [];
turbidity_un = [];
SPM = [];
SPM_un = [];
bb_red = [];
bb_red_un = [];
for ii = 1:numrows
    
    % Turbidity
    tmp = C{27}{ii};
    idx1 = strfind(tmp,'>');
    idx2 = strfind(tmp,'±');
    if isempty(idx1)
        turbidity = [turbidity; str2num(tmp(1:idx2-2))];
    else
        turbidity = [turbidity; str2num(tmp(2:idx2-2))];
    end
    turbidity_un = [turbidity_un; str2num(tmp(idx2+1:end))];
    
    % SPM
    tmp = C{28}{ii};
    idx1 = strfind(tmp,'>');
    idx2 = strfind(tmp,'±');
    if isempty(idx1)
        SPM = [SPM; str2num(tmp(1:idx2-2))];
    else
        SPM = [SPM; str2num(tmp(2:idx2-2))];
    end
    SPM_un = [SPM_un; str2num(tmp(idx2+1:end))];
    
    % bb red
    tmp = C{29}{ii};
    idx1 = strfind(tmp,'>');
    idx2 = strfind(tmp,'±');
    if isempty(idx1)
        bb_red = [bb_red; str2num(tmp(1:idx2-2))];
    else
        bb_red = [bb_red; str2num(tmp(2:idx2-2))];
    end
    bb_red_un = [bb_red_un; str2num(tmp(idx2+1:end))];
    
end

% Create the HCData_Header cell array
HChead = cell(numrows+1,7);
HChead(1,:) = {'Datenum' Headers{24:29}};
HChead(2:numrows+1,1) = C{3};

% Fill in output variables
varargout{1} = [dt C{4} C{5} rrs turbidity SPM bb_red];
varargout{2} = HChead;
varargout{3} = [dt C{6} C{6} rrs_un turbidity_un SPM_un bb_red_un];
varargout{4} = cell2mat(C(6:23));
varargout{5} = Headers(6:23);

% If there are six ouputs variables, fill the last one with images
if nargout == 6
    idx = find(filename == '/' | filename == '\');
    if isempty(idx)
        path = '';
    else
        path = filename(1:idx(end));
    end
    for ii = 1:numrows
        images(ii,1).card = imread(strcat(path,C{1}{ii},'_',C{2}{ii},'_Ed.png'));
        images(ii,1).sky = imread(strcat(path,C{1}{ii},'_',C{2}{ii},'_Ls.png'));
        images(ii,1).water = imread(strcat(path,C{1}{ii},'_',C{2}{ii},'_Lw.png'));
    end
    varargout{6} = images;
end

return






