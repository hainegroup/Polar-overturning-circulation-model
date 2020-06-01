function [Station, DateTime, Latitude, Longitude, Press, Temp, Salinity] = read_Fram_Strait_section(filename)

%% Initialize variables.

delimiter = '\t';

%% Format string for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: text (%s)
%	column4: text (%s)
%   column5: text (%s)
%	column6: text (%s)
%   column7: text (%s)
%	column8: text (%s)
%   column9: text (%s)
%	column10: text (%s)
%   column11: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
data_inds = 130:length(dataArray{1}) ;   % NOTE!  This number depends on the header length for the cruise/datafile in question!
Station   =            dataArray{ 1}(data_inds) ;
DateTime  =            dataArray{ 2}(data_inds) ;
Latitude  = str2double(dataArray{ 3}(data_inds)) ;
Longitude = str2double(dataArray{ 4}(data_inds)) ;
Elevation = str2double(dataArray{ 5}(data_inds)) ; 
Depth     = str2double(dataArray{ 6}(data_inds)) ;
Press     = str2double(dataArray{ 7}(data_inds)) ;
Temp      = str2double(dataArray{ 8}(data_inds)) ;
Potemp    = str2double(dataArray{ 9}(data_inds)) ;
Salinity  = str2double(dataArray{10}(data_inds)) ;
Sigma_t   = str2double(dataArray{11}(data_inds)) ;


%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;