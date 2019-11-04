% Plot crossBiplots (from older script)

PathForFiguresExport = 'C:\Users\rcardot\Documents\MATLAB\Figures\ConfluencesClassification\CrossComparisonBetweenDataSet\';

SaveResultsAndFigures='yes';

SuffixForExport='_NewQrAndQMetricsMAINTRIB';

NbComponentsToKeep = size(HydrologyClusters(:,3:end),2);

for ProcessingOrder=1:2
    
    if ProcessingOrder==1
        ListOfClusters = unique(Clusters.MORPHO_Cluster);
        ListOfClusters(isnan(ListOfClusters)) = [];
        PrefixForExport = 'Geomorpho_On_QrNat';
        DataSet = GeomorphClusters(:,{'MORPHO_JunctionNumber','MORPHO_Cluster'});
        DataSet.Properties.VariableNames = {'JunctionNumber','Cluster'};
    elseif ProcessingOrder==2
        ListOfClusters = unique(Clusters.GEOMET_Cluster);
        ListOfClusters(isnan(ListOfClusters)) = [];
        PrefixForExport = 'Geometric_On_QrNat';
        DataSet = GeometryClusters(:,{'GEOMET_JunctionNumber','GEOMET_Cluster'});
        DataSet.Properties.VariableNames = {'JunctionNumber','Cluster'};
    end
    
    %% Get unique clusters and associate a color
    
    % Store the number of clusters
    NbClusters=length(ListOfClusters);
    
    if NbClusters <= 8
        ListOfColors = brewermap(NbClusters,'Set1');
    elseif NbClusters <= 12
        ListOfColors = brewermap(NbClusters,'Paired');
    else
        ListOfColors = distinguishable_colors(length(ListOfClusters));
    end
    
    
    
    
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
        Xminimum=min(HydrologyClusters{...
            :,num2str(i,'HYDROL_PC%d')});
        
        Xmaximum=max(HydrologyClusters{...
            :,num2str(i,'HYDROL_PC%d')});
        % index 'j' is for the following principal components (y-axis)
        for j=(i+1):NbComponentsToKeep
            Yminimum=min(HydrologyClusters{...
                :,num2str(j,'HYDROL_PC%d')});
            Ymaximum=max(HydrologyClusters{...
                :,num2str(j,'HYDROL_PC%d')});
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
                
                %                 % Get the index of the corresponding gravity center
                %                 idxGravCentersThisCluster=...
                %                     [GravityCenters.Cluster]==ListOfClusters(k);
                %
                %                 % Draw the gravity center of the k-th cluster
                %                 GravCenterMark=plot(GravityCenters{...
                %                     idxGravCentersThisCluster,i+1},...
                %                     GravityCenters{...
                %                     idxGravCentersThisCluster,j+1});
                %                 % Formatting the markers style
                %                 GravCenterMark.Marker='s';
                %                 GravCenterMark.LineStyle='none';
                %                 GravCenterMark.LineWidth=1;
                %                 GravCenterMark.MarkerEdgeColor=[0.5 0.5 0.5];
                %                 GravCenterMark.MarkerSize=8;
                %                 GravCenterMark.MarkerFaceColor=ListOfColors(k,:);
                %
                % Get the index of the corresponding observations
                idxObsThisCluster = ...
                    [DataSet.Cluster] == ListOfClusters(k);
                
                % Get the Name of observations in this clusters
                
                Names = DataSet{idxObsThisCluster,'JunctionNumber'};
                
                % Find the corresponding idx in the Qr_nat dataset
                
                idxObsToPlot=ismember([HydrologyClusters.HYDROL_JunctionNumber],Names);
                
                if sum(idxObsToPlot) == 0
                    continue
                end
                
                % Draw the observations included in the k-th cluster
                BiPlot=plot(HydrologyClusters{...
                    idxObsToPlot,num2str(i,'HYDROL_PC%d')},...
                    HydrologyClusters{...
                    idxObsToPlot,num2str(j,'HYDROL_PC%d')},...
                    'DisplayName',num2str(ListOfClusters(k),...
                    'Cluster n°%d'));
                % Formatting the markers style
                BiPlot.Marker='o';
                BiPlot.LineStyle='none';
                BiPlot.MarkerSize=5;
                BiPlot.MarkerFaceColor=ListOfColors(k,:);
                BiPlot.MarkerEdgeColor = 'w';
                %                 BiPlot.MarkerEdgeColor='w';
            end
            
            % Formatting the axes
            ax=gca;
            xlim([Xminimum-abs(0.1*Xminimum) Xmaximum+abs(0.1*Xmaximum)])
            ylim([Yminimum-abs(0.1*Yminimum) Ymaximum+abs(0.1*Ymaximum)])
            %             ax.XLimMode='auto';
            %             ax.YLimMode='auto';
            axis square
            set(ax,'YAxisLocation','origin','XAxisLocation','origin')
            
            ax.XLabel.String=num2str(i,'PC%d');
            ax.YLabel.String=num2str(j,'PC%d');
            ax.YLabel.FontSize=11;
            ax.XLabel.FontSize=11;
            ax.FontSize=ax.XLabel.FontSize;
            %             ax.Title.FontSize=13;
            
            
            FigureIsFull = 1 ;%...
            %                 mod(idxSubPlot,4)==0 || (i*j)/2 == endLoop ;
            
            % Add a title to the figure each 4 iterations
            if FigureIsFull
                % title([PrefixForExport,SuffixForExport],'Color','k','FontSize',13,'Interpreter','none')
                
            end
            
            % Export each fourth iterations if asked
            if strcmpi(SaveResultsAndFigures,'yes') && FigureIsFull
                
                % Create a new file if its the first export
                if idxExport==1
                    export_fig([PathForFiguresExport,PrefixForExport,...
                        'Biplots_',datestr(now,'YYYYmmDD'),...
                        SuffixForExport],'-pdf','-p0.1');
                    idxExport=idxExport+1;
                    %                     close
                    % Else, append to the already existing file
                else
                    export_fig([PathForFiguresExport,PrefixForExport,'BiplotsGravCenters_',datestr(now,'YYYYmmDD'),SuffixForExport],'-pdf','-p0.1','-append');
                    idxExport = idxExport + 1;
                    %                     close
                end
                
            end
            % Incrementation of the subplot index
            idxSubPlot=idxSubPlot+1;
        end
    end
    
end

