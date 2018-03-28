% Script to plot Stacked timeseries by cluster

% User selection dialog  box for the file resulting from clustering
[FileWithClusters,PathOfFilesWithClusters] = uigetfile('*.mat',...
    'Select the file with the clusters');

load([PathOfFilesWithClusters,FileWithClusters])


% get the folder where the Qr TimeTables are stored as MAT-files
WorkingFolder = uigetdir(pwd,'Folder with Qr files');

%
% OutPutFIGFolder = uigetdir(pwd,'Folder for axes export');

tic % Start timer

% Get list of Qr files
ListFiles = dir([WorkingFolder,'\Qr*.mat']);

% Load the PixelsID XLSX-File to display the confluence name on the plots
PixID = readtable('C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\PIxels_ID_CorrectedMar2018_corrected_v2.xlsx');

% Store and change directory
StorePWD = pwd;
cd(WorkingFolder);

% Creation of the figure
fig = figure;
MaximizeFigureWindow;
% set(fig,'Visible','off')

% Creation of a tab group
TabGroup = uitabgroup;

ListCluster = unique(ClusteringResultAndScores.Cluster);

% Loop parameters
StartLoop = 1;
EndLoop = length(ListCluster); % 2;
NbFiles = height(ClusteringResultAndScores);
FileNumber = 1;
hWaitBar = waitbar(0,'');


% One iteration per Cluster
for i = StartLoop : EndLoop
    
    timerCluster = tic;
    
    ListJunctionNumber = sortrows(ClusteringResultAndScores.JunctionNumber(ClusteringResultAndScores.Cluster == ListCluster(i)));
    
    % Creates a new tab with empty axes
    TabName = num2str(ListCluster(i),'Cluster n°%i');
    %     FigOrFiles
    disp('******************************************')
    disp(['********* Processing ',TabName,' *********'])
    Tab = uitab(TabGroup,'title',TabName);
    Tab.UserData = ListCluster(i);
    ax = axes('Parent',Tab);
    
    for j = 1:length(ListJunctionNumber)
        % Get the strings of the file name
        FileName = num2str(ListJunctionNumber(j),'Qr_%d.mat');
        idxFile = strcmp({ListFiles.name},FileName);
        
        %     JunctionNumber = sscanf(TabName,'Qr_%d');
        PlotName = TabName;
        
        waitbar(FileNumber/NbFiles,hWaitBar,num2str([FileNumber,NbFiles],'Processing File %d on %d'))
        
        
        
        
        disp([num2str([FileNumber,NbFiles],'Processing File %d on %d'),' (',ListFiles(idxFile).name,')'])
        
        
        
        
        % Load the Qr Time Table
        load(ListFiles(idxFile).name)
        Qr(:,3) = array2timetable(repmat(FileNumber,height(Qr),1),'RowTimes',Qr.Properties.RowTimes);
        Qr.Properties.VariableNames(3) = {'FileNumber'};
        
        if j == 1
            TempQr = Qr;
        else
            TempQr = [TempQr;Qr];
        end
        
        FileNumber= FileNumber + 1;
    end
    
    
    disp([char(10),'*** COMPUTING QUANTILES AND PLOTTING... ***'])
    
    % Plot the values as a tome serie within the predefined axes
    PlotStackedYearlyTimeseriesForStackedClusters(ax,TempQr,4,PlotName);
    
    %     clear TempQr Qr
    
    % %     savefig(ax,[OutPutFIGFolder,'\',PlotName]);
    
    % Iteration of the processed file
    
    disp('**************** Elapsed Time for this cluster :')
    toc(timerCluster)
    
end
fig.Visible ='on';
cd(StorePWD)

disp('*************** Total elapsed time:')
toc