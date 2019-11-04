% Script To parse Spatial statistics from a SHP-File containing them in 


SHPFilePath='C:\Users\rcardot\Documents\SIG\TOPKAPI_ETH\StatsAtJunctions\';
SHPFileName='AllCatchmentsWithElevAndLandCoverStats_Mar2018.shp';

ShapeToProcess=shaperead(strcat(SHPFilePath,SHPFileName));
ShapeToProcess = struct2table(ShapeToProcess);
ShapeToProcess (:,1:4) = [];
ShapeToProcess.Properties.VariableNames{1} = 'PourBasinID';

VarNames = ShapeToProcess.Properties.VariableNames;
MainVarNames = strcat('Main_',VarNames([1,4:end]));
TribVarNames = strcat('Trib_',VarNames([1,4:end]));
NewVarNames = {'JunctionNumber',MainVarNames{:},TribVarNames{:},'AreaRatio'};

nbJunctions = length(unique(ShapeToProcess.JunctionNu));
nbVars = length(NewVarNames);


% Pre-allocate a table
% % NaNArray = nan(nbJunctions,1);
SummTable = table();
% % 
% % for i = 1 : nbVars
% %    SummTable(:,i) = NaNArray;  
% % end
% % 
% % clear NaNArray
% % SummTable.Properties.VariableNames = NewVarNames;


%% Parse data (trib/main)
idxFill = 1;

for i = 1: max(ShapeToProcess.JunctionNu)
    SummTable(idxFill,1) = {i};
    
    idxJunction = ShapeToProcess.JunctionNu == i;
    idxMain = idxJunction & ShapeToProcess.IsTrib == 0;
    idxTrib = idxJunction & ShapeToProcess.IsTrib == 1;
    
    if sum(idxMain) == 1 && sum(idxTrib) == 1
        SummTable(idxFill,2) = ShapeToProcess(idxMain,'PourBasinID');
        SummTable(idxFill,3:20) = ShapeToProcess(idxMain,4:21);
        
        SummTable(idxFill,21) = ShapeToProcess(idxTrib,'PourBasinID');
        SummTable(idxFill,22:39) = ShapeToProcess(idxTrib,4:21);
        
        SummTable(idxFill,40) = {ShapeToProcess{idxMain,'Area'} ./ ShapeToProcess{idxTrib,'Area'}};
        
    idxFill = idxFill + 1;
    end
    
    
end

SummTable.Properties.VariableNames = NewVarNames;