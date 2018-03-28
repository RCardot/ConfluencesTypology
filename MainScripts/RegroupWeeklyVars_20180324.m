% Script to regroup variables within the SummTable
%

PathOfFilesAndPrefix=readtable('C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\ComputedQr_Mar2018_v3\PathOfFilesForPCA_v3_WithWeeks.xls');
load(PathOfFilesAndPrefix.PathOfDataSet{1})


if exist('RegroupSummTable','var')
    clear RegroupSummTable
end


AvgingIndexes = {2,[59,7:16],17,[18:20],[21:23],24,[25:26],[27:28],[29:30],[31:32],[33:35],[36:40],[41:44],[45:46],[47,48],[49:50],[51:52],[53:54],[55:58]};

% correct an error while typing the above cell array
% for i = 2:length(AvgingIndexes)
%     AvgingIndexes{i} = AvgingIndexes{i} - 1 ;
% end

RegroupSummTable(:,1) = SummTable.JunctionNumber;

% Computing Variable Names

RegroupVarNames = cell(1,length(AvgingIndexes) + 1);
RegroupVarNames{1} = 'JunctionNumber';

for i =1 : length(AvgingIndexes)
    
    if numel(AvgingIndexes{i}) == 1
        RegroupVarNames(i+1) = SummTable.Properties.VariableNames(AvgingIndexes{i});
    else
        VarsNumbers = AvgingIndexes{i};
        StartVar= VarsNumbers(1);
        EndVar= VarsNumbers(end);
        
        idxStartNumString = regexp(SummTable.Properties.VariableNames{StartVar},'\d');
        strStartWeek = SummTable.Properties.VariableNames{StartVar}(idxStartNumString);
        
        idxEndNumString = regexp(SummTable.Properties.VariableNames{EndVar},'\d');
        strEndWeek = SummTable.Properties.VariableNames{EndVar}(idxEndNumString);
        
        idxBaseString = 1:idxStartNumString-1;
        BaseVarName = SummTable.Properties.VariableNames{StartVar}(idxBaseString);
        
        NewVarName = strcat(BaseVarName,strStartWeek,'_',strEndWeek);
        
        RegroupVarNames{i+1} = NewVarName;
        
        clear BaseVarName strEndWeek strStartWeek NewVarName
    end
    
    
end

% Compute means (NOTE: It is not correct to compute mean of STD. A step should be added in the pre processing,
% where n, sum of x, and sum of x^2, are computed. it will allow a better computation of aggregated STD and
% mean)
for i=1: length(AvgingIndexes)
    
    RegroupSummTable(1:height(SummTable),i+1) = mean(SummTable{:,AvgingIndexes{i}},2);
    
end

RegroupSummTable = array2table(RegroupSummTable);
RegroupSummTable.Properties.VariableNames = RegroupVarNames;

RegroupSummTable = [RegroupSummTable,SummTable(:,end-11:end)];
% Add Standard deviations

