function [Qr]=DischargeRatioTimeTable(Qmain,Qtrib,MainStemPosition,VariableName)
% Computes the discharge ratio 'Qtributary/Qmainstem'
% The two Inputs are timestable variables with the same time step
% The function finds the common timerange and computes the Qr
% 'MainStemPosition' takes two values 'U' (Upstream - Default) or 'D' (Downstream)
% VariableName is optionnal. It is a string for the desired VariableName.
%
% **************************************** R. CARDOT - 01.02.2017


% If there is no third argument it assumed that the mainstem discharge
% value is given upstream
if nargin<4
    VariableName='Qr';
    
    if nargin<3
        MainStemPosition='U';
    end
end

% Get the Start and End time of each time serie
StartQ1=Qmain.Time(1);
EndQ1=Qmain.Time(end);
StartQ2=Qtrib.Time(1);
EndQ2=Qtrib.Time(end);


% Determination of the common date interval
if StartQ1<=StartQ2
    if EndQ1>=EndQ2
        TimeRange=timerange(StartQ2,EndQ2);
    else
        TimeRange=timerange(StartQ2,EndQ1);
    end
else
    if EndQ1>=EndQ2
        TimeRange=timerange(StartQ1,EndQ2);
    else
        TimeRange=timerange(StartQ1,EndQ1);
    end
end


%% Computation of the discharge ratio

% By substraction, then division, ifthe Main Stem discharge is given
% downstream
if strcmpi(MainStemPosition,'D')
    Qr=Qtrib{TimeRange,1}./(Qmain{TimeRange,1}-Qtrib{TimeRange,1});
    
    % By only dividing, else
else
    Qr=Qtrib{TimeRange,1}./Qmain{TimeRange,1};
end

% Output as timetable
Qr=timetable(Qmain.Time(TimeRange),Qr); %,'VariableNames',{VariableName});


end