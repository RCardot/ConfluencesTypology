% Script to format computed Qr TimeTables by:
% (1) Rounding the timestamps in XX:59:59 to (XX+1):00:00
% (2) Adding a 'Day of year' field for further representations/plots
% (3) removing the first months of the dataset
%
% Tested on 18 Jan 2018 : Works perfectly on the already computed Qr TimeTables
%                         It is VERY slow however
%
%******************************************** R. CARDOT - 17 Jan 2018 *****


% get directories with input and output files 
% FolderWithQrTimeTables = uigetdir(pwd,'Folder with IN files');
FolderWithQrTimeTables = 'C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\ComputedQr_Mar2018_v3';

% FolderForFormattedFiles = uigetdir(pwd,'Folder with Output files');
FolderForFormattedFiles = [FolderWithQrTimeTables,'\','BetterFormattedTimeTables'];

if ~exist(FolderForFormattedFiles,'dir')
   mkdir(FolderForFormattedFiles) 
end

tic
clc
% Display summary of input and output directories
disp(['Input directory: ',FolderWithQrTimeTables,char(10),char(10),...
    'Output directory: ',FolderForFormattedFiles,char(10)])

% get list of 'Qr*.mat' files in the input directory
ListFiles = dir([FolderWithQrTimeTables,'\Qr*.mat']);

% One iteration per Qr file
for i = 1:length(ListFiles)
    
    % display the current processing file
    disp([num2str([i,length(ListFiles)],'Processing file %d on %d'),' (',ListFiles(i).name,')'])
    
    % Load the Qr file
    load([FolderWithQrTimeTables,'\',ListFiles(i).name]); % Qr TimeTables are saved under a variable named 'Qr'
   
    % Remove dates before the 1st Jan 1991
    Qr = Qr(Qr.Time >= datetime(1991,01,01),:);

    % Retime (interpolate) for 'wrong time stamps correction
    Qr = retime(Qr,'hourly','linear'); % resampling at the exact time (i.e. XX:00:00)
    
    % Add the day of year field
    Qr.DayOfYear = DayOfYear(Qr.Time); % Computing the Day of year
    
    % Save the new Qr file
    save([FolderForFormattedFiles,'\',ListFiles(i).name],'Qr')
      
end

toc