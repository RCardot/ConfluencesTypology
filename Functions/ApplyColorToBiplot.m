function ApplyColorToBiplot(hbi,Grouping,ListOfGroups,ColorSet,SizeForMarkers)

% Function to apply colors on data points plotted on a biplot
% Manipulate plot colors

% Code inspired by the thread:
% URL:[https://stackoverflow.com/questions/19052317/configuring-biplot-in-matlab-to-distinguish-in-scatter]
% and the code : https://pastebin.com/KHUj3DnA
%
% INPUTS:
% hbi: handle of an already plotted boxplot of the data.
% while plotting 'hbi' the parameter 'obslabel' should be filled with the
% same argument 'groups'

% Grouping : is a matrix of n eelements corresponding to each point in the
% bi plot. Each element of the array is the group associated to the
% observation point, which is linked to the color which has to be
% represented

% ListOfGroups: is the list of m groups given as 'obslabel' when plotting the
% biplot. It is a an array containing groups names as numbers

% ColorSet is the set of colors to be associated to the list of groups. It
% is a mx3 matrix of colors defined as a RGB triplet between 0 and 1.
% Each row corresponds to desired color for each group in 'groups'
%
% SizeForMarkers: Is a scaler defining the size for observation points
% markers. The Unit used is the unit defined withe the parameter 'Units'
% (default: 'points')


% NOTE: This script works if there is no text attached to observations.
% Else it will change the color of the text instead of the color of the
% observation markers.



idxPoint=1;
for ii = length(hbi)-length(Grouping):length(hbi)-1
    groupOfThePoint = Grouping(idxPoint);
    if ~isempty(groupOfThePoint)
        hbi(ii).Color=ColorSet(ListOfGroups==groupOfThePoint,:);
        hbi(ii).MarkerSize=SizeForMarkers;
    else
        hbi(ii).Color=[0.3 0.3 0.3];
        hbi(ii).MarkerSize=SizeForMarkers;
    end
    idxPoint=idxPoint+1;
end

% disp(num2str(idxPoint-1))