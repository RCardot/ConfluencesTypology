% Script with clean PCA Workflow appliable on the different datasets

% Work In Progress

% To Do: Code the choice to enter the values interactively or directly in
% the code

%% Clear WorkSpace and command window
close all

% % % If at least one of the parameters already exists, keep them all for
% % % pre-fill
% % if exist('CutoffMethodForClustering','var')
% %     clearvars -except CutoffMethodForClustering SaveResultsAndFigures ...
% %         WorkWithAllVariables standardization NbComponentsToKeep
% % end
clearvars
clc


%% Select interactive or 'hard-coding' mode for parameters selection

% Interactive='no'; % Not implemented yet

%% Parameters for computation

SaveResultsAndFigures='yes';
if strcmpi(SaveResultsAndFigures,'yes')
    QstSuffix=questdlg('Use a suffix for exported files ?','Suffix');
    if strcmpi(QstSuffix,'yes')
        SuffixForExport=cell2mat(inputdlg(...
            'Which suffix do you want to use?',...
            'Suffix string'));
    else
        SuffixForExport='';
    end
else
    SuffixForExport='';
end


% WorkWithAllVariables='yes';
% StringPatternsOfVarsToIgnore={'Slope','Elev_1','90'}; %Fill this cell array with patterns of Variables names, which have to be ignored in th analysis
standardization='yes';
% NbComponentsToKeep=3; % Number of principal componenets for plotting and clustering

% Setup for clustering
CutoffMethodForClustering='distance'; % put either 'maxclust' or 'distance'

% NbOfClusters=7; % Should be set regarding the dendrogram

ScreenSize=get(groot,'ScreenSize'); % get screen size to fit figures to screen

RoundingLevel=4; % Number of digits to keep in exported results (e.g. PC loadings)

clear QstSuffix
%% Path of files for import and export
PathOfFilesAndPrefix=readtable('C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\ComputedQr_Mar2018_v3\PathOfFilesForPCA_v3_WithWeeks.xls');

PathForExport='C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\ComputedQr_Mar2018_v3\StatsOnQr\Figures\PCA\';

if not(exist(PathForExport,'dir'))
    mkdir(PathForExport)
end

 
% PathOfMATFileVarsToRemove='C:\Users\rcardot\Documents\SIG\TOPKAPI_ETH\StatsAtJunctions\VariablesToRemoveGeometry.mat';

ObsNamesFieldName='JunctionNumber'; % Name of the field which contains the name of observations in the original table dataset

% Shape Files
PathOfSHPFileOfJunctionsPoints='C:\Users\rcardot\Documents\SIG\RHONE_CATCHMENT\Hydrographic Network\gwn_n_selected_Mar2018_1903Plus.shp';
PathOfSHPFileOfHydroNetwork='C:\Users\rcardot\Documents\SIG\RHONE_CATCHMENT\Hydrographic Network\StreamsKanalLakesWithOrder.shp';
PathOfSHPFileBasinPolyg='C:\Users\rcardot\Documents\SIG\ShapeFiles\RhoneCatchmentLimits1903Plus.shp';
PathOfSHPFileBasinPolyLine='C:\Users\rcardot\Documents\SIG\ShapeFiles\RhoneCatchmentLimitsPolyLine1903Plus.shp';


%% Process whole analysis on each dataset
%  (one loop iteration('idxAnalysis' incrementation) per dataset)


nTotalFiles=height(PathOfFilesAndPrefix);

for idxAnalysis=1:nTotalFiles
    
    % Display the processing progress (i.e. 'File n/nTotalFiles
    disp(num2str([idxAnalysis,nTotalFiles],'Processing file %d/%d'))
    
    
    % Check if there is a Cell array defined with Variables to ignore
    WorkWithAllVariables='yes'; %initialisation of the variable
    try 
        isnan(PathOfFilesAndPrefix.PathOfVariablesToIgnoreArray(idxAnalysis));
    catch % if there is an error using 'isnan', it means that there is a string in the excel table
        WorkWithAllVariables='no';
    end
    
    
    %% Load dataset and other parameters from the 'PathOfFilesAndPrefix' table
    
    % load the observation/values dataset
    load(PathOfFilesAndPrefix.PathOfDataSet{idxAnalysis})
    
    % Load the patterns of variable strings to ignore (stored as a .mat
    % file) and set 'WorkWithAllVariables' state to 'yes' or 'no'
    if strcmp(WorkWithAllVariables,'no') 
        load(PathOfFilesAndPrefix.PathOfVariablesToIgnoreArray{idxAnalysis}); % load VarsToIgnore cell array        
    end
    
    % Removal of rows with NaN and/or Inf values
    if any(strcmp('QrMean',SummTable.Properties.VariableNames))
        SummTable((isnan([SummTable.QrMean])|isinf([SummTable.QrMean])),:)=[];
    end
    
    % Get the prefix to use for export and add 'standardized' if
    % standardization is 'yes'
    PrefixForExport=PathOfFilesAndPrefix.PrefixForExport{idxAnalysis};
    if strcmpi(standardization,'yes')==1
        PrefixForExport=strcat(PrefixForExport,'_standardized');
    end
    
    clear DataSetFilePath
    
    %% Variables Selection
    
    % Retrieve index of variables to process (+ exclude Variable with
    % observations labels)
    if strcmpi(WorkWithAllVariables,'no')==1
        StringPatternsOfVarsToIgnore{1,end+1}=ObsNamesFieldName;
    else
        StringPatternsOfVarsToIgnore={ObsNamesFieldName};
    end
    
    % Compute the index of variables to keep in the processed dataset
    idxVarsToKeep=not(contains(SummTable.Properties.VariableNames,...
        StringPatternsOfVarsToIgnore));
    
    % Prepare the dataset to be processed and give RowNames corresponding to
    % observation labels
    DataSet=SummTable(:,idxVarsToKeep);
    DataSet.Properties.RowNames=cellstr(num2str(SummTable{:,ObsNamesFieldName}));
    
    clear SummTable idxVarsToKeep WorkWithAllVariables ...
        StringPatternsOfVarsToIgnore
    %% Prepare cell arrays for labelling observations and variables
    ObsNames=DataSet.Properties.RowNames;
    VarNames=DataSet.Properties.VariableNames;
    
    %% Plotting the boxplots of all variables
    figure('Visible','off')
    
    if strcmpi(standardization,'yes') % If standardization is asked
        % boxplot is plotted with zscore of the dataset (i.e.: standardized
        % values)
        boxplot(zscore(DataSet{:,:}),'orientation','horizontal',...
            'labels',VarNames); %,'Notch','on','PlotStyle','compact')
    else
        % else with the original dataset
        boxplot(DataSet{:,:},'orientation','horizontal','labels',VarNames,...
            'Notch','on'); %,'PlotStyle','compact')
    end
    
    % Formatting the axes
    title([PrefixForExport,SuffixForExport],'Interpreter','none')
    ax=gca;
    Yaxis=ax.YAxis;
    Yaxis.FontSize=8;
    
    % Export the Figure if asked
    if strcmpi(SaveResultsAndFigures,'yes')==1
        export_fig([PathForExport,PrefixForExport,...
            'BoxPlotVariables_',datestr(now,'YYYYmmDD'),...
            SuffixForExport],'-pdf','-p0.1');
    end
    
    clear YAxis
    %% Compute PCA
    % C = corr(DataSet,DataSet);
    
    if strcmpi(standardization,'yes') % If standardization is asked
        
        % PCA is computed with zscore of the dataset (i.e.: standardized
        % values)
        [wcoeff,score,latent,tsquared,explained,mu] = ...
            pca(zscore(DataSet{:,:}),'Centered','on');
    else
        % else with the original dataset
        [wcoeff,score,latent,tsquared,explained,mu] = ...
            pca(DataSet{:,:},'Centered','on');
    end
    
    %% Create orthonormal matrix of coefficient (IS THIS STEP NEEDED?!!!)
    % (/!\ THIS SHOULD BE VERIFIED: IS IT NEEDED TO ORTHONORM COEFF + DOES IT CHANGE THE INTERPRETATION)
    if strcmpi(standardization,'yes')==1
        coefforth = inv(diag(std(zscore(DataSet{:,:}))))*wcoeff;
    else
        coefforth = inv(diag(std(DataSet{:,:})))*wcoeff;
    end
    
    %% Plot a pareto graph to show explained variance by Principal components
    figure(2)
    
    [~,ParAxes]=pareto(explained);
    
    % Formating the figure
    ParAxes(2).YTickLabel={}; % Remove the default right y-axis label
    xlabel('Principal Components')
    ylabel('Variance Explained (%)')
    title([PrefixForExport,SuffixForExport],'Interpreter','none')
    
    % Ask user for the number of components to keep in further clustering
    % regarding the pareto graph
    NbComponentsToKeep=str2double(inputdlg(...
        'Number of components to keep ?',...
        'Number of components'));
    
    % Export figure
    if strcmpi(SaveResultsAndFigures,'yes')
        export_fig([PathForExport,PrefixForExport,'Pareto_',...
            datestr(now,'YYYYmmDD'),SuffixForExport],...
            '-pdf','-p0.1');
        close
    end
    
    
    clear ParAxes
    %% Summarize loadings of variables on chosen components
    
    % Create a table
    PCLoadings=table();
    
    % Fill the first columns with the loadings of variables (rows) on each
    % Principal Component (columns)
    for i=1:NbComponentsToKeep
        PCLoadings{:,i}=round(coefforth(:,i),RoundingLevel);
        PCLoadings.Properties.VariableNames{i}=num2str(i,'PC%i');
    end
    
    % Fill the next columns with the squared loadings of variables (rows) on each
    % Principal Component (columns)
    for i=NbComponentsToKeep+1:2*NbComponentsToKeep
        PCLoadings{:,i}=round(coefforth(:,i-NbComponentsToKeep).^2,...
            RoundingLevel);
        PCLoadings.Properties.VariableNames{i}=num2str(i-NbComponentsToKeep,...
            'squaredloads_PC%i');
    end
    
    % Fill the next columns with the computed t-value (p=0.95) of variables
    % loadings (rows) on each Principal Component (columns)
    for i=2*NbComponentsToKeep+1:3*NbComponentsToKeep
        
        % Define a function which computes the t-student value for a given
        % correlation coefficient 'r' and a number of observations 'n'
        tStudentOnLoading=@(r,n) r.*sqrt((n-2)./(1-r.^2));
        
        % Fill columns of the table with a rounded value of the computed
        % t-student
        PCLoadings{:,i}=round(tStudentOnLoading(coefforth(:,...
            i-2*NbComponentsToKeep),...
            height(DataSet)),RoundingLevel);
        
        % Give a name to the column
        PCLoadings.Properties.VariableNames{i}=...
            num2str(i-2*NbComponentsToKeep,'tStud_PC%i');
    end
    
    % Computes the t-student limit for our number of observations
    tLimit=tinv(0.95,height(DataSet)-2);
    rLimit=sqrt(tLimit.^2/(tLimit.^2 + height(DataSet) - 2));
    
    % Fill the last column with the t-student limit computed for the number
    % of observations           (/!\ See for vectorization)
    for i=1:height(PCLoadings)
        PCLoadings{i,3*NbComponentsToKeep+1}=round(rLimit,RoundingLevel);
    end
    
    PCLoadings.Properties.VariableNames{3*NbComponentsToKeep+1}='rLimit';
    
    % Fill the row names of the table with the variables names
    PCLoadings.Properties.RowNames=VarNames;
    
    % Export in XLS- and MAT-files
    if strcmpi(SaveResultsAndFigures,'yes')==1
        writetable(PCLoadings,[PathForExport,PrefixForExport,...
            'PCLoadings_',datestr(now,'YYYYmmDD'),...
            SuffixForExport,'.xls'],'WriteRowNames',1);
        save([PathForExport,PrefixForExport,'PCLoadings_',...
            datestr(now,'YYYYmmDD'),SuffixForExport,'.mat'],'PCLoadings')
    end
    
    
    
    %% Plot bargraphes of loadings and squared loadings
    
    % SubPlot with Loadings
    figure('Position',[0 35 ScreenSize(3) ScreenSize(4)-110],'Visible','off')
    subplot(2,1,1)
    bar(PCLoadings{:,1:NbComponentsToKeep},1,'hist');
    BarGraph1=gca;
    title([PrefixForExport,'_loadings',SuffixForExport],'Interpreter','none')
    
    %     hold on
    %     line(
    
    % Formating the axes
    BarGraph1.TickLabelInterpreter = 'none';
    BarGraph1.TitleFontSizeMultiplier=1;
    BarGraph1.XLim=[0 height(PCLoadings)+1];
    BarGraph1.XTick=(1:height(PCLoadings));
    BarGraph1.XTickLabelMode='manual';
    BarGraph1.XTickLabelRotation=90;
    BarGraph1.XLabel.Interpreter='none';
    BarGraph1.XAxis.FontSize=12;
    BarGraph1.XTickLabel=VarNames;
    BarGraph1.XGrid='on';
    BarGraph1.XAxis.TickLength=[0 0];
    
    % Create the legend labels strings
    for i=1:NbComponentsToKeep
        legendText{i}=num2str(i,'PC%d');
    end
    
    % Display the legend
    legend(legendText,'Box','off','Position',[0.021 0.496 0.057 0.122]);
    
    % SubPlot with Squared loadings
    subplot(2,1,2)
    bar(PCLoadings{:,NbComponentsToKeep+1:2*NbComponentsToKeep},1,'hist');
    BarGraph2=gca;
    title([PrefixForExport,'_squared_loadings',SuffixForExport],...
        'Interpreter','none')
    
    % Formating the axes
    BarGraph2.TickLabelInterpreter = 'none';
    BarGraph2.TitleFontSizeMultiplier=1;
    BarGraph2.XLim=[0 height(PCLoadings)+1];
    BarGraph2.XTick=(1:height(PCLoadings));
    BarGraph2.XTickLabelMode='manual';
    BarGraph2.XTickLabelRotation=90;
    BarGraph2.XLabel.Interpreter='none';
    BarGraph2.XAxis.FontSize=12;
    BarGraph2.XTickLabel=VarNames;
    BarGraph2.XGrid='on';
    BarGraph2.XAxis.TickLength=[0 0];
    
    % Export the Figure if asked
    if strcmpi(SaveResultsAndFigures,'yes')==1
        export_fig([PathForExport,PrefixForExport,...
            'BarGraphLoadings_',datestr(now,'YYYYmmDD'),...
            SuffixForExport],'-pdf','-p0.1');
    end
    
    clear BarGraph1 BarGraph2
    
    
    %% Plot the dendrogram
    figure
    
    % Compute the distance between observations on the selected
    % Principal Components
    Z=linkage(score(:,1:NbComponentsToKeep),'ward','euclidean');
    dendrogram(Z,0)
    
    set(gca,'YGrid','on','YMinorGrid','on')
    title([PrefixForExport,SuffixForExport],'Interpreter','none')
    
    set(gcf,'Visible','on')
    % set(gca,'ColorOrder',distinguishable_colors(30))
    
    % Ask for the cutoff distance or Nb of clusters (in function of the
    % chosen cut-off method (see header lines of the script)
    
    % if the chosen method is distance between clusters
    if strcmpi(CutoffMethodForClustering,'distance')
        
        % interactively ask for a distance
        CutOffDistance=str2double(inputdlg(['Which distance'...
            '(y-axis of the dendrogram) for clustering cut-off ?'],...
            'Choose clustering cut-off'));
        
        % draw the colored dendrogram
        dendrogram(Z,0,'colorthreshold',CutOffDistance)
        
        % draw a horizontal line corresponding to the chosen cut-off distance
        line([0,length(ObsNames)+1],[CutOffDistance CutOffDistance],...
            'LineStyle','--','Color','k');
        
        % If the cut-off method is the number of clusters
    elseif strcmpi(CutoffMethodForClustering,'maxclust')==1
        
        % Interactively ask for the number of clusters
        NbOfClusters=str2double(inputdlg('How many clusters?',...
            'Choose clustering cut-off'));
        
        % get the value of the distance between the 'NbOfClusters'-th
        % cluster and the previous one, in the output of the 'linkage'
        % function (i.e. the third column of Z, Aat row
        % 'end-NbOfClusters+1') and add a small amount (0.0001) to cut
        % above this value.
        CutOffDistance=Z(end-NbOfClusters+1,3)+0.0001;
        
        % draw the colored dendrogram with selected clusters
        dendrogram(Z,0,'colorthreshold',CutOffDistance)
        
        % draw a horizontal line at the selected cut-off distance
        line([0,length(ObsNames)+1],[CutOffDistance CutOffDistance],...
            'LineStyle','--','Color','k');
    end
    
    % Formatting the axes
    set(gca,'YGrid','on','YMinorGrid','on')
    set(gca,'XTickLabel',{})
    title([PrefixForExport,SuffixForExport],'Interpreter','none')
    
    % Export the figure if asked
    if strcmpi(SaveResultsAndFigures,'yes')
        export_fig([PathForExport,PrefixForExport,'Dendrogram_',...
            datestr(now,'YYYYmmDD'),SuffixForExport],'-pdf','-p0.1');
    end
    
    %% Clustering the data with the threshold chosen above
    
    % if the cut-off method is 'distance'
    if strcmpi(CutoffMethodForClustering,'distance')==1
        T=clusterdata(score(:,1:NbComponentsToKeep),'cutoff', CutOffDistance,...
            'criterion','distance','linkage','ward');   %'maxclust',30);
        
        % if the cut-off method is 'maxclust'
    elseif strcmpi(CutoffMethodForClustering,'maxclust')==1
        T=clusterdata(score(:,1:NbComponentsToKeep),'maxclust',NbOfClusters,...
            'linkage','ward');
    end
    
    % %      % Next block of the code is used for display in the command window
    % %      % (maybe useless)
    % %     clusters=unique(T);
    % %
    % %     for i=1:length(clusters)
    % %
    % %         CountByCluster(2*i-1)=i;
    % %         CountByCluster(2*i)=sum(T==clusters(i));
    % %
    % %     end
    % %
    % %     disp(num2str(CountByCluster,'Cluster n°%d : %d confluences \n'))
    
    % Create a summary table of clusters allocated to each observation
    ClusteringResultAndScores=table();
    ClusteringResultAndScores{:,1}=str2double(ObsNames);
    ClusteringResultAndScores{:,2}=T;
    
    % Filling the table with the observation scores on desired Principal
    % Components
    for i=1:NbComponentsToKeep
        ClusteringResultAndScores{:,i+2}=score(:,i);
        AdditionalVarNames{i}=num2str(i,'PC%d');
    end
    
    % set the name of columns
    ClusteringResultAndScores.Properties.VariableNames =...
        [{'JunctionNumber','Cluster'},AdditionalVarNames];
    
    % Export the summary table if asked
    if strcmpi(SaveResultsAndFigures,'yes')
        writetable(ClusteringResultAndScores,...
            [PathForExport,PrefixForExport,'_Clusters_',...
            datestr(now,'YYYYmmDD'),SuffixForExport,'.xls'])
        save([PathForExport,PrefixForExport,...
            '_ClustersAndScores_',datestr(now,'YYYYmmDD'),SuffixForExport,...
            '.mat'],'ClusteringResultAndScores')
    end
    
    close(gcf)
    
    
    %% Get unique clusters and associate a color
    
    ListOfClusters=unique(T);
    
    % Store the number of clusters
    NbClusters=length(unique(T));
    
    if NbClusters <= 8
        ListOfColors = brewermap(NbClusters,'Set1');
    elseif NbClusters <= 12
        ListOfColors = brewermap(NbClusters,'Paired');
    else
        ListOfColors=distinguishable_colors(length(ListOfClusters));
    end
    
    
    
    %% Computation of the gravity center of each cluster
    
    % Create a table to summarize the gravity centers coordinates of each
    % cluster
    GravityCenters=table();
    GravityCenters.Cluster=ListOfClusters;
    
    for i = 1:NbClusters
        idx=[ClusteringResultAndScores.Cluster]==ListOfClusters(i);
        
        for j=1:NbComponentsToKeep
            GravityCenters{i,j+1}=mean(ClusteringResultAndScores{...
                idx,AdditionalVarNames{j}});
        end
    end
    
    % Set the name of the columns of the table
    GravityCenters.Properties.VariableNames=[{'Cluster'},strcat('MeanOn'...
        ,AdditionalVarNames(:)')];
    
    clear AdditionalVarNames idx
    
    
    %% Find the representative observation for each Cluster
    
    % Pre-allocate a table
    ClosestObsName=table(nan(NbClusters,1),nan(NbClusters,1));
    ClosestObsName.Properties.VariableNames={'Cluster','ObsName'};
    
    % Look for the closest observation to the gravity center of the i-th cluster
    for i=1:NbClusters
        
        idxObservations = [ClusteringResultAndScores.Cluster]==ListOfClusters(i);
        
        % Coordinates of the observations for the i-th cluster
        PointCloudCoord=ClusteringResultAndScores{idxObservations,...
            3:3+NbComponentsToKeep-1};
        
        % Coordinates of the gravity center for the i-th cluster
        GravityCenterCoord=GravityCenters{...
            [GravityCenters.Cluster]==ListOfClusters(i),2:end};
        
        % Find the idx of the closest observation to the gravity center
        idxClosestObs = ClosestPoint(PointCloudCoord,GravityCenterCoord);
        
        if sum(idxClosestObs) ~= 1
            disp(num2str([sum(idxClosestObs),num2str(ListOfClusters(i),'%d')],...
                'There are %d representative elements for Cluster n°%d, The first is kept.'))
            idxOBS = find(idxClosestObs == 1,1,'first');
            idxClosestObs(:) = deal(0);
            idxClosestObs(idxOBS) = 1;
        end
        
        % Fill the table with the name of the found closest observation
        ClusterObsNames=ClusteringResultAndScores{idxObservations,'JunctionNumber'};
        ClosestObsName.ObsName(i)=ClusterObsNames(idxClosestObs);
        ClosestObsName.Cluster(i)=ListOfClusters(i);
        
        clear PointCloudCoord GravityCenterCoordinates idxObservations ...
            idxClosestObs ClusterObsNames
    end
    
     % Export the summary table if asked
    if strcmpi(SaveResultsAndFigures,'yes')
        writetable(ClosestObsName,...
            [PathForExport,PrefixForExport,'_RepresentativObs_',...
            datestr(now,'YYYYmmDD_HHMMSS'),SuffixForExport,'.xls'])
        save([PathForExport,PrefixForExport,'_RepresentativObs_',...
            datestr(now,'YYYYmmDD_HHMMSS'),SuffixForExport,...
            '.mat'],'ClosestObsName')
    end
    
    
    % Export %
    
    
    %% Mapping clusters on a study site map
    
    %     if ~exist('JunctionsPoint','var')
    JunctionsPoint=shaperead(PathOfSHPFileOfJunctionsPoints);
    %     end
    % JunctionsPoint.Cluster=zeros(size(JunctionsPoint,1),1);
    
    %     j=1;
    for i=1:length(ObsNames)
        try
            JunctionsPoint([JunctionsPoint.NUMERO]==str2double(ObsNames(i,:))).Cluster=T(i);
            %             emptyPoints=false;
        catch
            disp(['No point for confluence n°',ObsNames(i,:)])
            %             RowsToRemove(j)=find([JunctionsPoint.NUMERO]==str2double(ObsNames(i,:)));
            %             j=j+1;
            %             emptyPoints=true;
        end
    end
    
    % Remove points which are not included in the analysis
    j=1;
    for i=1:length(JunctionsPoint)
        if ~ismember([JunctionsPoint(i).NUMERO],str2double(DataSet.Properties.RowNames))
            RowsToRemove(j)=i;j=j+1;
        end
    end
    
    
    if exist('RowsToRemove','var')
        JunctionsPoint(RowsToRemove)=[];
    end
    
    clear RowsToRemove
    %% Mapping the results
    
    % Initialize the figure
    figure('Visible','on')
    MaximizeFigureWindow
%     set(gcf,'Position',[0 35 ScreenSize(3) ScreenSize(4)-110]);
    hold on
    
    % Load the Hydrographic Network shape file if it is not already loaded
    if ~exist('HydroNetwork','var')
        HydroNetwork=shaperead(PathOfSHPFileOfHydroNetwork);
        
        % Define the line style and color for Lakes and the rest of the
        % network. The field 'FLOZ' corresponds to the Strahler Order of
        % the reach
        HydroNetworkStyle=makesymbolspec('Line',...
            {'OBJECTVAL','See','Color',[0.53 0.81 0.98]},...
            {'FLOZ',1,'LineWidth',0.25},...
            {'FLOZ',2,'LineWidth',0.5},...
            {'FLOZ',3,'LineWidth',0.75},...
            {'FLOZ',4,'LineWidth',1},...
            {'FLOZ',5,'LineWidth',1.25},...
            {'FLOZ',6,'LineWidth',1.5},...
            {'FLOZ',7,'LineWidth',1.75},...
            {'FLOZ',8,'LineWidth',2.5},...
            {'Default','Color',[0.25 0.25 0.25]});
    end
    
    % Load the catchment shape file if it is not already loaded
    if not( exist('CatchmentPolyg','var') || exist('CatchmentPolyLine','var') )
        
        
        CatchmentShapeDisplay = 1; % 1 : Outline, 2 : Only Face, 3 : Outline + Face
        
        
        switch CatchmentShapeDisplay
            case 1
                % Load the catchment outline SHP-File
                CatchmentPolyLine=shaperead(PathOfSHPFileBasinPolyLine); 
                
                % Define the line style and color for The catchment outline
                CatchmentPolyLineStyle=makesymbolspec('Line',{'Default','Color',[0.7 0.2 0.2]});
                
                % Draw the catchment polygon
                mapshow(CatchmentPolyLine,'SymbolSpec',CatchmentPolyLineStyle);
                
            case 2
                % Load the catchment polygon SHP-File
                CatchmentPolyg=shaperead(PathOfSHPFileBasinPolyg);
                
                % Define the line style and Facecolor for The catchment polygon
                CatchmentPolygStyle=makesymbolspec('Patch',...
                    {'Default','FaceColor',[0.98 0.98 0.85],...
                    'EdgeColor','none'});
                
                 % Draw the catchment polygon
                 mapshow(CatchmentPolyg,'SymbolSpec',CatchmentPolygStyle);
            
            case 3
                % Load the catchment polygon SHP-File
                CatchmentPolyg=shaperead(PathOfSHPFileBasinPolyg);
                
                % Define the line style and Facecolor for The catchment polygon
                CatchmentPolygStyle=makesymbolspec('Patch',...
                    {'Default','FaceColor',[0.98 0.98 0.85],...
                    'EdgeColor',[0 0 0]});
                
                 % Draw the catchment polygon
                 mapshow(CatchmentPolyg,'SymbolSpec',CatchmentPolygStyle);
        end
        
    end
    
   
    
    % Draw the Hydrographic network
    mapshow(HydroNetwork,'SymbolSpec',HydroNetworkStyle);
    
    
    
    
    % Add colored points corresponding to junctions (One color per cluster)
    % One loop iteration per cluster
    
    idxMappedClusters=1;
    for i=1:NbClusters
        
        % Set the points color to the color affected to the selected
        % cluster
        StyleOfPoints=makesymbolspec('Point',...
            {'Default','Marker','o',...
            'MarkerFaceColor',ListOfColors(i,:),...
            'MarkerEdgeColor',ListOfColors(i,:)});
        
        % Draw junctions points alocated to the selected cluster
        
        %         if emptyPoints(i)
        %             continue
        %         else
        
        MappedPoints(idxMappedClusters)=mapshow(JunctionsPoint(...
            [JunctionsPoint(:).Cluster]==ListOfClusters(i)),...
            'SymbolSpec',StyleOfPoints);
        
        
        %         end
        
        % Add a 'DisplayName' to the drawn points for the map legend
        
        NbOfObservationsInCluster=sum([JunctionsPoint(:).Cluster]==ListOfClusters(i));
        
        MappedPoints(idxMappedClusters).DisplayName=num2str([ListOfClusters(i),...
            NbOfObservationsInCluster],...
            'Cluster %d (%d)');
        if not(isempty([JunctionsPoint(:).Cluster]==ListOfClusters(i)))
            idxMappedClusters=idxMappedClusters + 1;
        end
        
    end
    
    % Add textlegend with junction number
    
    for idxTxt = 1: length(JunctionsPoint)
        text(JunctionsPoint(idxTxt).X + 40,JunctionsPoint(idxTxt).Y + 30, ...
            num2str(JunctionsPoint(idxTxt).NUMERO),'Color',[0.8 0.7 0.7],'FontWeight','bold',...
            'FontSize',10)
    end
    
    
    % Formatting the axes
    mapClusters=gca;
    mapClusters.XLabel.String='X - meters (CH1903+)';
    mapClusters.YLabel.String='Y - meters (CH1903+)';
    title(['Classification of confluences based on ',PrefixForExport,...
        SuffixForExport,char(10),...
        num2str(...
        NbComponentsToKeep,'Number of Principal Components: %i')],...
        'Interpreter','none')
    
    % Show the legend for mapped points
    try % Added on 2017.09.12 to avoid an new error occuring with a empty cluster (i.e. no junction associated)
        hleg=legend(MappedPoints);
        hleg.Location='southeast';
    catch
        disp('Error while creating the map legend')
    end
    %     hleg.
    
    % Set AspectRatio to 1 (i.e. same scale for x and y coordinates)
    set(gca,'DataAspectRatio',[1 1 1])
    
    % Export the map if asked
    if strcmpi(SaveResultsAndFigures,'yes')
        export_fig([PathForExport,PrefixForExport,'ClustersMap_',...
            datestr(now,'YYYYmmDD'),SuffixForExport],...
            '-pdf','-p0.03');
    end
    
    clear hleg mapClusters MappedPoints
    
    
    %% Plotting Biplots 
    % This version uses the internal matlab 'biplot' function, which
    % plots together variables and observations and (/but) scales data
    
    figure
    
    coefforth = PCLoadings{:,1:2};
    score = ClusteringResultAndScores{:,3:4};
    VarNames = PCLoadings.Properties.RowNames;
    T = ClusteringResultAndScores.Cluster;
    ListOfClusters = unique(T);
    
    BiPlot=biplot(coefforth,'scores',score,'MarkerSize',6,'varlabels',VarNames);%,'obslabels',str2mat(strcat('Junction n°',ObsNames,'-',num2str(T,'Clust.n°%i'))));%,'UserData',T);
    
    ApplyColorToBiplot(BiPlot,T,ListOfClusters,ListOfColors,12) % Apply a color grouping on observation points
    
    % Changing the colors of variable lines and markers
    for i=1:2*length(VarNames)
    
        BiPlot(i).Color=[0.85 0.85 0.85];
        BiPlot(i).LineStyle='--';
    
    
    end
    
    % Changing the color and size of Variable labels
    
    for i=2*length(VarNames)+1:3*length(VarNames)
    
        BiPlot(i).Color=[0.4 0.4 0.4];
        BiPlot(i).FontSize=6;
    
    end
    
    axis square
    
    
    title([PrefixForExport,SuffixForExport],'Interpreter','none')
    
    if strcmpi(SaveResultsAndFigures,'yes')
        export_fig([PathForExport,PrefixForExport,'_BiplotsWithVars_',datestr(now,'YYYYmmDD'),SuffixForExport],'-pdf','-p0.05');
    end
    
    if NbComponentsToKeep>2
        figure()
        coefforth = PCLoadings{:,1:3};
    score = ClusteringResultAndScores{:,3:5};
    
        BiPlot2=biplot(coefforth(:,[1,3]),'scores',score(:,[1,3]),'MarkerSize',6,'varlabels',VarNames,'obslabels',str2mat(strcat('Junction n°',ObsNames,'-',num2str(T,'Clust.n°%i'))));
        ylabel('Component 3')
    
        ApplyColorToBiplot(BiPlot2,T,ListOfClusters,ListOfColors,12) % Apply a color grouping on observation points
    
        % Changing the colors of variable lines and markers
        for i=1:2*length(VarNames)
    
            BiPlot2(i).Color=[0.85 0.85 0.85];
            BiPlot2(i).LineStyle='--';
    
        end
    
        % Changing the color and size of Variable labels
    
        for i=2*length(VarNames)+1:3*length(VarNames)
    
            BiPlot2(i).Color=[0.4 0.4 0.4];
            BiPlot2(i).FontSize=6;
    
        end
    
        axis square
        title(PrefixForExport)
        if strcmpi(SaveResultsAndFigures,'yes')
            export_fig([PathForExport,PrefixForExport,'_BiplotsWithVars_',...
                datestr(now,'YYYYmmDD'),SuffixForExport],'-pdf','-p0.05','-append');
        end
    end
    
    clear coefforth score VarNames T
   
    %% Plotting biplots (2x2 axes on one figure)
    % This version only plots observations but not variables, and uses
    % a call to the 'scatter' function instead of the 'biplot' function
    
    % Computes the number of iterations for plotting biplots
    % 'endLoop' corresponds to the number of principal components minus one
    %
    % For each principal component, one biplot is drawn with this component
    % as x-axis and each following one (i.e. with higher number) as y-axis
    
    endLoop=sum(1:NbComponentsToKeep-1);
    
    
    figure
    
    % First indexing of subplot
    idxSubPlot=1;
    
    % First indexing for export (this index is used to append all the drawn
    % biplots in the same pdf file)
    idxExport=1;
    
    % Loop to plot observations on the different Factorial plane/design
    
    % index 'i' is for the principal component on x-axis
    for i=1:NbComponentsToKeep-1
        
        % index 'j' is for the following principal components (y-axis)
        for j=(i+1):NbComponentsToKeep
            
            % reset the subplot index if it is over 4
            if idxSubPlot==5 || idxSubPlot==1
                idxSubPlot=1;
                %                 figure('Visible','off')
            end
            
            % Select the subplot, where to draw the plot
            subplot(2,2,idxSubPlot)
            hold on
            
            % Plot the observations points on the (i-th PC, j-th PC) plane
            % there is one 'k' iteration per cluster
            for k=1:length(ListOfClusters)
                
                % Get the index of the corresponding gravity center
                idxGravCentersThisCluster=...
                    [GravityCenters.Cluster]==ListOfClusters(k);
                
                % Draw the gravity center of the k-th cluster
                GravCenterMark=plot(GravityCenters{...
                    idxGravCentersThisCluster,i+1},...
                    GravityCenters{...
                    idxGravCentersThisCluster,j+1});
                % Formatting the markers style
                GravCenterMark.Marker='s';
                GravCenterMark.LineStyle='none';
                GravCenterMark.LineWidth=1;
                GravCenterMark.MarkerEdgeColor=[0.5 0.5 0.5];
                GravCenterMark.MarkerSize=8;
                GravCenterMark.MarkerFaceColor=ListOfColors(k,:);
                
                % Get the index of the corresponding observations
                idxObsThisCluster = ...
                    [ClusteringResultAndScores.Cluster]==ListOfClusters(k);
                
                % Draw the observations included in the k-th cluster
                BiPlot=plot(ClusteringResultAndScores{...
                    idxObsThisCluster,num2str(i,'PC%d')},...
                    ClusteringResultAndScores{...
                    idxObsThisCluster,num2str(j,'PC%d')},...
                    'DisplayName',num2str(ListOfClusters(k),...
                    'Cluster n°%d'));
                % Formatting the markers style
                BiPlot.Marker='o';
                BiPlot.LineStyle='none';
                BiPlot.MarkerSize=5;
                BiPlot.MarkerFaceColor=ListOfColors(k,:);
                BiPlot.MarkerEdgeColor='w';
            end
            
            % Formatting the axes
            ax=gca;
            ax.XLimMode='auto';
            ax.YLimMode='auto';
            axis square
            set(ax,'YAxisLocation','origin','XAxisLocation','origin')
            
            ax.XLabel.String=num2str([i,explained(i)],'PC%d(%.1f%%)');
            ax.YLabel.String=num2str([j,explained(j)],'PC%d(%.1f%%)');
            ax.YLabel.FontSize=11;
            ax.XLabel.FontSize=11;
            ax.FontSize=ax.XLabel.FontSize;
            %             ax.Title.FontSize=13;
            
            
            FigureIsFull = ...
                mod(idxSubPlot,4)==0 || (i*j)/2 == endLoop ;
            
            % Add a title to the figure each 4 iterations
            if FigureIsFull
                suptitle([PrefixForExport,SuffixForExport]) %,'k',13)%,'none')
                %                 hold off
            end
            
            % Export each fourth iterations if asked
            if strcmpi(SaveResultsAndFigures,'yes') && FigureIsFull
                
                % Create a new file if its the first export
                if idxExport==1
                    export_fig([PathForExport,PrefixForExport,...
                        'BiplotsGravCenters_',datestr(now,'YYYYmmDD'),...
                        SuffixForExport],'-pdf','-p0.05');
                    idxExport=idxExport+1;
                    close
                    % Else, append to the already existing file
                else
                    export_fig([PathForExport,PrefixForExport,'BiplotsGravCenters_',datestr(now,'YYYYmmDD'),SuffixForExport],'-pdf','-p0.05','-append');
                    idxExport = idxExport + 1;
                    close
                end
            end
            % Incrementation of the subplot index
            idxSubPlot=idxSubPlot+1;
        end
    end
    
    
    %% Create a log file with the main parameters chosen in processing
    
    if idxAnalysis==1
        SummaryOfProcessing=table();
    end
    
    SummaryOfProcessing{idxAnalysis,1}={PrefixForExport};
    SummaryOfProcessing{idxAnalysis,2}=NbComponentsToKeep;
    SummaryOfProcessing{idxAnalysis,3}={CutoffMethodForClustering};
    SummaryOfProcessing{idxAnalysis,4}=NbClusters;
    SummaryOfProcessing{idxAnalysis,5}=CutOffDistance;
    
    SummaryOfProcessing.Properties.VariableNames={'DataSet',...
        'NbOfPrincComponents','CutOffMethod','NbOfClusters','CutOffDistance'};
    % Close figures and clear memory before processing a new file
%     close all
    clear ListOfClusters ListOfColors T   BiPlot CountByCluster ...
        CutOffDistance DataSet GravCenterMark legendText mapClusters ...
        score tsquared VarNames wcoeff Z clusters coefforth ax ...
        endLoop explained figTitle i k latent mu ...
        NbComponentsToKeep ObsNames PrefixForExport ...
        StringPatternsOfVarsToIgnore tLimit Yaxis idxClosestObs ...
        idxExport idxSubPlot j StyleOfPoints % ClusterObsNames
end

writetable(SummaryOfProcessing,[PathForExport,datestr(now,'YYYYmmDD'),...
    SuffixForExport,'.xls'])