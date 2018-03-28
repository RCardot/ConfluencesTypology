clearvars

close all



% New version: WorkInProgress



tic



PathofFiles='C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\ComputedQr_Mar2018_v3\BetterFormattedTimeTables\';

FileNamePattern='Qr_*.mat';



Files=dir([PathofFiles,FileNamePattern]);

nbFiles = length(Files);
% EmptyArrayForFilling = nan(nbFiles,1);


SummTable=array2table(nan(nbFiles,111));
VarNames = cell(1,111);
VarNames(1:7) ={'JunctionNumber','QrMean','QrMedian','Qr5ile','Qr25ile','Qr75ile','Qr95ile'};



fillCell = 8;
for i = 1:52
    VarNames(fillCell) = {num2str(i,'QrWeek%d')};
    fillCell = fillCell +1;
    
    VarNames(fillCell) = {num2str(i,'STDWeek%d')};
    fillCell = fillCell +1;
end

SummTable.Properties.VariableNames = VarNames;

clear VarNames fillCell
%%

YearStart=1991;

YearEnd=2008;



warning('off')



for i=1:nbFiles
    
    
    
    if mod(i,10)==0
        
        disp([num2str(i),'/',num2str(length(Files))])
        
    end
    
    
    
    % Loads a hourly timetable of discharge ratio into a variable called Qr
    
    load([PathofFiles,Files(i).name])
    
    idxYears = year(Qr.Time) >= YearStart & year(Qr.Time) <= YearEnd;
    
    Qr = Qr(idxYears,1);
    
    % retrieve the junction number from the filename
    
    JunctionNumber=regexp(Files(i).name,'[0-9]','match');
    
    JunctionNumber={string(str2double([JunctionNumber{1:end}]))}; % Converts into string
    
%     SummTable.Properties.RowNames(i) = JunctionNumber;
    
    SummTable.JunctionNumber(i) = str2double(JunctionNumber{:});
    
    clear JunctionNumber
    
    
    %% Total
    
    
    
    SummTable.QrMean(i)=nanmean(table2array(Qr));
    
    SummTable.QrMedian(i)=median(table2array(Qr));
    
    QrPrctiles=prctile(table2array(Qr),[5,25,75,95]);
    
    SummTable.Qr5ile(i)=QrPrctiles(1);
    
    SummTable.Qr25ile(i)=QrPrctiles(2);
    
    SummTable.Qr75ile(i)=QrPrctiles(3);
    
    SummTable.Qr95ile(i)=QrPrctiles(4);
    
    clear QrPrctiles
    
    
    %% Weekly
    
    QrWeekly = retime(Qr(:,1),'weekly','mean');
    
    Qr.WeekNumber = week(Qr.Time);
    
    
    %     QrSTD = retime(Qr(:,1),'weekly',@std);
    %     QrWeekNumber = retime(Qr(:,2),'weekly',@mode);
    
    %     QrWeekly.STD = QrSTD(:,1).Variables;
    %     QrWeekly.WeekNumber = QrWeekNumber(:,1).Variables;
    
    clear QrSTD QrWeekNumber %Qr
    
    %% Fill Summary Table
    
    for weekNb = 1 : 52
        
        SummTable.(num2str(weekNb,'QrWeek%d'))(i) = nanmean(Qr.Qr(Qr.WeekNumber == weekNb))./SummTable.QrMean(i);
        
        SummTable.(num2str(weekNb,'STDWeek%d'))(i) = nanstd(Qr.Qr(Qr.WeekNumber == weekNb))./SummTable.QrMean(i);
    end
    
    
    for idxMonth = 1:12
        
        SummTable.(num2str(idxMonth,'STDMonth%d'))(i) = nanstd(Qr.Qr(month(Qr.Time) == idxMonth))./SummTable.QrMean(i);
        
    end
    
    
    
end

% Reorder variables
SummTable = SummTable(:,[1:7 8:2:end-12 9:2:end-12 end-11:end]);

clear YearStart YearEnd FileNamePattern i nbFiles PathofFiles Files

toc