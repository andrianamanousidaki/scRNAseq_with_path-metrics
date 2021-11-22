function [ExtendedDbscanLabels,NumberDbscanClusters] = RunDBSCAN(X,epsilon,minpts)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

    DbscanLabels = dbscan(X,epsilon,minpts);

    labeled_idx = find(DbscanLabels~=-1);
    outlier_idx = find(DbscanLabels==-1);

    %% Use NN classifier to assign outliers

    ExtendedDbscanLabels = zeros(length(DbscanLabels),1);
    ExtendedDbscanLabels(labeled_idx) = DbscanLabels(labeled_idx);
    CorePoints = X(labeled_idx,:);
    CorePointLabels = DbscanLabels(labeled_idx);
    [NN1, D_NN1] = knnsearch(CorePoints,X(outlier_idx,:),'k',1);

    for i=1:length(outlier_idx)
        ExtendedDbscanLabels(outlier_idx(i)) = CorePointLabels(NN1(i));
    end
    
    NumberDbscanClusters = length(unique(ExtendedDbscanLabels));
    
end

