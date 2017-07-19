function [ wl, var, delta_var] = hp_import( filename_val, filename_dark, var_name )
%HP_IMPORT import HyperPro file
%   Return Lu or Ed

% Check input
if nargin > 3
   error('Too many input arguments')
elseif nargin < 2
   error('Not enough input arguments')
end

% Set param
start_line = 4;

% Read val file
[wl_g, var_g] = import_file(filename_val, var_name, start_line);

% Read dark file
[wl_k, var_k] = import_file(filename_dark, var_name, start_line);

% Compute median, delta
wl = wl_g';
var_g = median(var_g,2);
var_k = median(var_k,2);
delta_var = std(var_g,0,2);

% Interpolate dark
var_k = interp1(wl_k, var_k, wl);
var = var_g - var_k;

end

function [wl, var] = import_file(filename, var_name, start_line)
% Read file for a
[fid, errmsg] = fopen(filename);
if fid < 0; error(errmsg); end;

% Skip first part
for foo=1:start_line;
  tline = fgetl(fid);
end;

% Get wavelength
data = strsplit(tline, ',');
wl=[]; i_wl = [];
for i=5:size(data, 2);
  foo = sscanf(data{i}, sprintf('%s(%%f)', var_name));
  if ~isempty(foo)
    wl(end+1) = foo;
    i_wl(end+1) = i;
  else
    %warning('Wavelength not recognise %s', data{i});
  end;
end;

% Skip unit line
tline = fgetl(fid);

% Get data
var=[];
tline = fgetl(fid);
while ischar(tline)
  data = strsplit(tline, ',');
  var(:, end+1) = str2double(data(i_wl));
  tline = fgetl(fid);
end;
fclose(fid);

end