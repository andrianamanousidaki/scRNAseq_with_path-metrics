function [Geometric_Perturbation, EmbeddingCentroids] = GeometricPerturbation(X,U,Labels)

% input: data set X
% input: data set U, which is a low dimensional embedding of X
% input: Labels, predetermined groups in X
% output: measurement of global geometric perturbation of the groups when
% data forced into lower dimensions

k=length(unique(Labels));
uniqueGroups = unique(Labels);
EmbeddingCentroids = zeros(k, size(U,2) );
ClusterCentroidDistances = zeros(k,k);
ClusterCentroidDistancesFull = zeros(size(U,1),size(U,1));
PMClusterCentroidDistances = zeros(k,k);
PMClusterCentroidDistancesFull = zeros(size(U,1),size(U,1));
for i=1:length(uniqueGroups)
    ind_i = find(Labels==uniqueGroups(i));
    centroid_i = mean(X(ind_i,:));
    EmbeddingCentroids(i,:) = mean(U(ind_i,:));
    for j=i+1:length(uniqueGroups)
        ind_j = find(Labels==uniqueGroups(j));
        centroid_j = mean(X(ind_j,:));
        pm_centroid_j = mean(U(ind_j,:));
        ClusterCentroidDistances(i,j) = norm(centroid_i - centroid_j);
        ClusterCentroidDistances(j,i) = ClusterCentroidDistances(i,j);
        ClusterCentroidDistancesFull(ind_i,ind_j)=ClusterCentroidDistances(i,j);
        ClusterCentroidDistancesFull(ind_j,ind_i)=ClusterCentroidDistances(i,j);
        PMClusterCentroidDistances(i,j) = norm(EmbeddingCentroids(i,:) - pm_centroid_j);
        PMClusterCentroidDistances(j,i) = PMClusterCentroidDistances(i,j);
        PMClusterCentroidDistancesFull(ind_i,ind_j)=PMClusterCentroidDistances(i,j);
        PMClusterCentroidDistancesFull(ind_j,ind_i)=PMClusterCentroidDistances(i,j);
    end
end 
c = sum(sum( ClusterCentroidDistancesFull.*PMClusterCentroidDistancesFull ))/ sum(sum( PMClusterCentroidDistancesFull.^2 ));
Geometric_Perturbation = sum(sum( (ClusterCentroidDistancesFull - c*PMClusterCentroidDistancesFull).^2 ))/sum(sum( ClusterCentroidDistancesFull.^2));

end

