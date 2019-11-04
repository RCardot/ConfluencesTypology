% Aggregate clusters results from the three dataset

GeomorphClusters = load('C:\Users\rcardot\Documents\MATLAB\ConfluencesTypology\SampleFigures\PCAResults\Geomorph\Geomorph_standardized_ClustersAndScores_20180404.mat');
GeomorphClusters = GeomorphClusters.ClusteringResultAndScores;
GeomorphClusters.Properties.VariableNames = strcat('MORPHO_',GeomorphClusters.Properties.VariableNames);

GeometryClusters = load('C:\Users\rcardot\Documents\MATLAB\ConfluencesTypology\SampleFigures\PCAResults\Geometry\Geometry_standardized_ClustersAndScores_20180404.mat');
GeometryClusters = GeometryClusters.ClusteringResultAndScores;
GeometryClusters.Properties.VariableNames = strcat('GEOMET_',GeometryClusters.Properties.VariableNames);

HydrologyClusters = load('C:\Users\rcardot\Documents\MATLAB\ConfluencesTypology\SampleFigures\PCAResults\Qr_natural\Qr_natural_standardized_ClustersAndScores_20180404.mat');
HydrologyClusters = HydrologyClusters.ClusteringResultAndScores;
HydrologyClusters.Properties.VariableNames = strcat('HYDROL_',HydrologyClusters.Properties.VariableNames);

% Loop initiaélisation
maxJunctionNumber = max(unique([GeomorphClusters.MORPHO_JunctionNumber;GeometryClusters.GEOMET_JunctionNumber;HydrologyClusters.HYDROL_JunctionNumber]));
Clusters = table;
idxFill = 1;

% Loop : One iteration per junction
for i = 1:maxJunctionNumber
    
    % get index of the corresponding row in each table
    idxMorpho = GeomorphClusters.MORPHO_JunctionNumber == i;
    idxGeomet = GeometryClusters.GEOMET_JunctionNumber == i;
    idxHydrol = HydrologyClusters.HYDROL_JunctionNumber ==i;
    
    % There is at least one cluster value amongst the three dataset
    if sum([idxMorpho;idxGeomet;idxHydrol]) ~= 0
        Clusters(idxFill,1) = {i};
        
        % Fill with Morpho cluster if not empty, else fill with NaN
        if sum(idxMorpho) ~= 0
            Clusters(idxFill,2) = {GeomorphClusters.MORPHO_Cluster(idxMorpho)};
        else
            Clusters(idxFill,2) = {nan};
        end
        
        % Fill with Geomet cluster if not empty, else fill with NaN
        if sum(idxGeomet) ~= 0
            Clusters(idxFill,3) = {GeometryClusters.GEOMET_Cluster(idxGeomet)};
        else
            Clusters(idxFill,3) = {nan};
        end
        
        % Fill with Hydro cluster if not empty, else fill with NaN
        if sum(idxHydrol) ~= 0
            Clusters(idxFill,4) = {HydrologyClusters.HYDROL_Cluster(idxHydrol)};
        else
            Clusters(idxFill,4) = {nan};
        end
        
        idxFill = idxFill + 1;
    else
        continue
    end
end

Clusters.Properties.VariableNames = {'JunctionNumber','MORPHO_Cluster','GEOMET_Cluster','HYDROL_Cluster'};

