clearvars

tic % run a timer

% Load the matrix with all hourly run off data
load('C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\For_Romain\Q_nat.mat')

% converts to timetable
Q_dateTime=datetime(Q_nat.time(:),'ConvertFrom','datenum');
Q_natTable=timetable(Q_dateTime(:));

% Fill the timetable
for i=1:length(Q_nat.FID_sections)
    Q_natTable.(i)=Q_nat.Hourly_Discharge(:,i);
end

% Get the pixel numbers from simulation results
VarNames=Q_nat.FID_sections;
clear Q_nat

% Loading pixels ID from a file containing most of GIS extracted informations about studied confluences
PixID=importPixelsIDOnly('C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\PIxels_ID_CorrectedMar2018_corrected_v2.xlsx');

% % % % Line added for tests
% % % PixID = PixID(PixID.NUMERO == 86,:);
% % % 

% Computes the number of confluences
NbConfluences = height(PixID);

% One loop iteration per confluence (Computation of the discharge ration at each timestep
for i=1:NbConfluences
    
    disp(num2str([i NbConfluences],'Processing confluence %d on %d'))
    
    % Load cute runoff time serie at the pixel corresponding to the mainstem
    QMain = Q_natTable(:,VarNames==PixID.MainStemUp(i)); 
    
    % Loads cute runoff time serie at the pixel corresponding to the mainstem
    QTrib = Q_natTable(:,VarNames==PixID.TribUp(i)); 
    
    % Computes the discharge ratio
    Qr = DischargeRatioTimeTable(QMain,QTrib,'U');
    
    
    % Save computed data into MAT-File
    
%     save(['C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\ComputedQ_Nat_MAINTRIB_Jan2018\QNat_MAIN_',num2str(PixID.NUMERO(i)),'.mat'],'QMain') ;
%     save(['C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\ComputedQ_Nat_MAINTRIB_Jan2018\QNat_TRIB_',num2str(PixID.NUMERO(i)),'.mat'],'QTrib') ;
    save(['C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\ComputedQr_Mar2018_v3\Qr_',num2str(PixID.NUMERO(i)),'.mat'],'Qr') ;
    
    
    % Save computed data into CSV-File
%     writetable(timetable2table(QMain),['C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\ComputedQ_Nat_MAINTRIB_Jan2018\QNat_MAIN_',num2str(PixID.NUMERO(i)),'.csv'],'DateLocale','fr_FR')
%     writetable(timetable2table(QTrib),['C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\ComputedQ_Nat_MAINTRIB_Jan2018\QNat_TRIB_',num2str(PixID.NUMERO(i)),'.csv'],'DateLocale','fr_FR')
    writetable(timetable2table(Qr),['C:\Users\rcardot\Documents\Hydrologie\TOPKAPI_ETHZ\Simulation Results\ComputedQr_Mar2018_v3\Qr_',num2str(PixID.NUMERO(i)),'.csv'],'DateLocale','fr_FR')
    
    
    clear QMain QTrib Qr
end


toc