function PlotStackedYearlyTimeseriesForStackedClusters(ax,QrTimeTable,MovMeanSize,PlotTitle)


disp('**************** averaging values (movmean)')
QrTimeTable.DateTime = datetime(2000,01,01) + days(QrTimeTable.DayOfYear) + hours(hour(QrTimeTable.Time));
QrTimeTable.AveragedQr = movmean(QrTimeTable.Qr,MovMeanSize);

uniqDateTime = unique(QrTimeTable.DateTime);
HeightSynthesis = length(uniqDateTime);


Synthesis = table(uniqDateTime,nan(HeightSynthesis,1),nan(HeightSynthesis,1),...
    nan(HeightSynthesis,1),nan(HeightSynthesis,1),nan(HeightSynthesis,1),nan(HeightSynthesis,1));
Synthesis.Properties.VariableNames = {'DateTime','Mean','Median','Perc95ile','Perc05ile','Perc25ile','Perc75ile'};


ax.NextPlot = 'add';
hLegendEntries = gobjects(1,7);

uniqFileNb = unique(QrTimeTable.FileNumber);
disp(num2str(length(uniqFileNb),'**************** length of uniqFileNb : %d'))

% FileWithClusters
minYear = min(year(QrTimeTable.Time));
maxYear = max(year(QrTimeTable.Time));

disp(num2str([minYear,maxYear],'**************** minYear: %d - maxYear: %d'))

for File = 1:length(uniqFileNb)
    
    idxFile = QrTimeTable.FileNumber == uniqFileNb(File);
    
    disp(num2str([uniqFileNb(File),sum(idxFile)],'**************** Plotting FileNumber: %d, (nb. of entries: %d)'))
    
    for y =  minYear: maxYear
        
        idxYear = idxFile & (year(QrTimeTable.Time) == y);
        
        
        % Creates a legend entry if this is the first plotted year
        if File == 1 && y == minYear
            hLegendEntries(1) = plot(ax,QrTimeTable.DateTime(idxYear),...
                QrTimeTable.AveragedQr(idxYear),'Color',[0.7 0.7 0.8],...
                'LineWidth',0.6);%,'DisplayName',num2str(y,'Actual Year :%d'),...
                %'UserData',y);
        else
            plot(ax,QrTimeTable.DateTime(idxYear),QrTimeTable.AveragedQr(idxYear),...
                'Color',[0.7 0.7 0.8],'LineWidth',0.6);%,...
%                 'DisplayName',num2str(y,'Actual Year :%d'),...
%                 'UserData',y);
        end
        
    end
    
end

% myPercentiles = @x prctile(x,[50 95 5 25 75]); (This is an idea for improvments, but each year must be a
% variable of the TimeTable and each row should correspond  to DateTime) # 26 mar 2018


% % % % QrTimeTable.Time = QrTimeTable.DateTime;


disp('**************** COMPUTING PERCENTILES *****************')
for i = 1:height(Synthesis)
    
    if mod(i-1,24*26) == 0
       disp([num2str(round((i+1)/24),'**************** Percentile for day %d/366:  '),datestr(Synthesis.DateTime(i),'dd mmm')]) 
        
    end
    idxTime = QrTimeTable.DateTime == Synthesis.DateTime(i);
    
    Synthesis.Mean(i) = nanmean(QrTimeTable.AveragedQr(idxTime));
    Percentiles = prctile(QrTimeTable.AveragedQr(idxTime),[50 95 5 25 75]);
    
    Synthesis.Median(i) = Percentiles(1);
    Synthesis.Perc95ile(i) = Percentiles(2);
    Synthesis.Perc05ile(i) = Percentiles(3);
    Synthesis.Perc25ile(i) = Percentiles(4);
    Synthesis.Perc75ile(i) = Percentiles(5);
end

hLegendEntries(3) = plot(ax,Synthesis.DateTime,Synthesis.Median,...
    'Color',[1 0.1 0.1],'LineWidth',0.7);

hLegendEntries(7) = plot(ax,Synthesis.DateTime,Synthesis.Perc95ile,...
    'Color',[1 0.7 0.1],'LineWidth',0.7);

hLegendEntries(4) = plot(ax,Synthesis.DateTime,Synthesis.Perc05ile,...
    'Color',[1 0.7 0.1],'LineWidth',0.7);

hLegendEntries(2) = plot(ax,Synthesis.DateTime,Synthesis.Mean,...
    'Color',[1 0.5 0.2],'LineWidth',0.5);

hLegendEntries(5) = plot(ax,Synthesis.DateTime,Synthesis.Perc25ile,...
    'Color',[0.9 0.9 0.1],'LineWidth',0.7);

hLegendEntries(6) = plot(ax,Synthesis.DateTime,Synthesis.Perc75ile,...
    'Color',[0.9 0.9 0.1],'LineWidth',0.7);



legend(hLegendEntries(1:5),{'Qr time series','Mean','Median','5/95-ile','25/75-ile'});%,'75','95'});

% xlim([datetime(2000,01,01) datetime(2000,12,31)])
ylim([0 1.05*nanmax(Synthesis.Perc95ile)])
ax.XMinorTick = 'on';
ax.YGrid='on';
% dynamicDateTicks(ax);

% datetick('x','dd mmm') % Removes the display of the year (e.g. 2000), but disable the automatic zoom on x-axids

title(PlotTitle,'Interpreter','none')