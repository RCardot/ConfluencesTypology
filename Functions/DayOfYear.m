function outDayOfYear = DayOfYear(DateTime)

% This functions return the day of year (e.g.: 2 for the 2nd January, 33 for the 2nd of February)
%
% It takes leap year into account, which leads to a 366th day every 4 years
%
% NOTE: This function is very slow as I was not able to vectorize the processing
%
%**************************************************** R. CARDOT - 17 Jan 2018


DaysInMonthsNormalYear = [0,31,28,31,30,31,30,31,31,30,31,30]; % no need for december
DaysInMonthsLeapYear = [0,31,29,31,30,31,30,31,31,30,31,30]; % no need for december
outDayOfYear = nan(size(DateTime));
LeapYear = isleap(year(DateTime));


for i=1:length(DateTime)
 
    if LeapYear(i)
        outDayOfYear(i) = sum(DaysInMonthsLeapYear(1:month(DateTime(i)))) + day(DateTime(i));
    else
        outDayOfYear(i) = sum(DaysInMonthsNormalYear(1:month(DateTime(i)))) + day(DateTime(i));
    end   
    
end

end