% This script shows how the data sets were preprocessed; output stored in scRNAseq/DataSets/ProcessedDataSets

clear all
addpath(genpath('../DataSets'))

% Save results?
SaveResults = 'yes';

% Choose Data Set (for running stand-alone, not needed for RunPathMetrics)
%DataSet = 'RNAMix1';
%DataSet = 'RNAMix2';
%DataSet = 'CellMix';
%DataSet ='CellMixSng';
%DataSet = 'Beta';
%DataSet = 'Beta2';
%DataSet = 'BaronPanc';
%DataSet = 'TMPanc';
%DataSet = 'TMLung';
%DataSet='PBMC3k';
%DataSet='PBMC3k_SingleR'
DataSet='PBMC4k_main'
%DataSet='PBMC4k_extented'
%DataSet = 'Manifold';

% Choose Normilaztion
Normalization = 'Basic';
%Normalization = 'Linnorm';
%Normalization = 'Seurat';
%Normalization = 'SCT';

% Plot figs?
PlotFigs = 0;

% Processing & Denoising Options (Defaults used in paper)
DataOpts.EliminateGenes = 'yes'; %Options: 'yes' or 'no'
DataOpts.NumberOfGenes = 2000; %If DataOpts.EliminateGenes = 'yes', number of genes to keep
if strcmp(Normalization,'Basic')
    DataOpts.RescaleHighVarGenes = 'yes'; %Options: 'yes' or 'no'
else
    DataOpts.RescaleHighVarGenes = 'no';
end
DataOpts.PCA = 'yes'; %Options: 'yes' or 'no'; reduce data dimension via PCA
DataOpts.DM = 'no'; %Options: 'yes' or 'no'; reduce data dimension via Diffusion Maps
DataOpts.NumberOfPCs = 40; %If DataOpts.PCA='yes', number of principal components to use
DenoisingOpts.RemoveOutliers = 'no'; %Options: 'yes' or 'no' -> should outliers be removed?
DenoisingOpts.OutlierRemovalKNN = 10;
DenoisingOpts.LocalAvgAllPts = 'no'; %Options: 'yes' or 'no'
DenoisingOpts.LocalAvgNbhdSize = 12;

LoadRNAData

%% Remove tiny clusters (optional, used for Pancreatic data)
if strcmp(DataSet,'PBMC4k_main')
    Labels(Labels==3)=2;  
end
%%
H=hist(Labels,unique(Labels));
if exist('MinClusterSampleSize','var')
    u=unique(Labels);
    clusters_to_keep = u(find(H>=MinClusterSampleSize));
    points_to_keeps_idx = find(Labels==clusters_to_keep(1));
    for i=2:length(clusters_to_keep)
        points_to_keeps_idx = union(points_to_keeps_idx, find(Labels==clusters_to_keep(i)));
    end
    Data = Data(points_to_keeps_idx,:);
    Labels = Labels(points_to_keeps_idx,:);
else
    clusters_to_keep = find(H>=1); %Keep all clusters
end
%%
  
% Rename the labels to be consecutive integers (needed for  correct calc of
% LabelsAligned in GetAccuracies, which is needed for geometric fidelity)

for i=1:length(clusters_to_keep)
    Labels( Labels==clusters_to_keep(i) ) = i;
end

%% Optional: restrict to high variance genes; rescale super high variance genes

if strcmp(DataOpts.EliminateGenes,'yes')
    sorted_var = sort(var(Data),'descend');
    cutoff = sorted_var(DataOpts.NumberOfGenes+1);
    high_var_genes = find((var(Data)>cutoff));
    X = Data(:,high_var_genes);
else
    X = Data;
end

if PlotFigs==1
  
    figure
    plot(sort(var(X)))
    title('Original variances')
    
end

%% Optional: rescale super high variance genes

if strcmp(DataOpts.RescaleHighVarGenes,'yes')
    % See which genes are outliers in terms of very high variance:
    if isempty(DataOpts.VarCutoff)
        var_cutoff = prctile(var(X),75)+1.5*(prctile(var(X),75)-prctile(var(X),25));
    else
        var_cutoff = DataOpts.VarCutoff;
    end
    super_high_var_genes = find(var(X)>var_cutoff);
    X(:,super_high_var_genes) = sqrt(var_cutoff)*X(:,super_high_var_genes)*diag(1./sqrt(var(X(:,super_high_var_genes)))); 
    
    if PlotFigs==1
        figure
        plot(sort(var(X)))
        title('Rescale Very High Var Genes')
    end
end


%% Optional: use PCA to reduce dimension

[coeff2,score2,latent2] = pca(X,'NumComponents',min(DataOpts.NumberOfPCs,size(X,2)));

if PlotFigs==1
    figure
    view(3)
    grid on
    hold on
    axis equal
    rotate3d on
    % Plot each group individually: 
    for q = 1:length(uniqueGroups)
          % Get indices of this particular unique group:
          ind = Labels==uniqueGroups(q); 
          % Plot only this group: 
          plot3(score2(ind,1),score2(ind,2),score2(ind,3),'.','markersize',10,'Color',Colors(q,:)); 
    end
    title('PCA embedding of processed data','fontsize',14)
    xlabel('PC_1','fontsize',14)
    ylabel('PC_2','fontsize',14)
    zlabel('PC_3','fontsize',14)
    %add 1 to the labels to match numbering of https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-019-0425-8/MediaObjects/41592_2019_425_MOESM1_ESM.pdf
    %legend('group 1-2','group 3','group 4','group 5','group 6','group 7','group 8')
end

if strcmp(DataOpts.PCA,'yes')
    X = score2(:,1:DataOpts.NumberOfPCs);
elseif strcmp(DataOpts.DM,'yes')
    X = DiffusionMapsEmbedding(X,DataOpts.NumberOfPCs,15);
end


%% Denoise Data (optional)

if strcmp(DenoisingOpts.RemoveOutliers,'yes')
    
    [IDX, D_KNN] = knnsearch(X,X,'k',DenoisingOpts.OutlierRemovalKNN);
    sorted_KNN = sort(D_KNN(:,end));
    %cutoff =prctile(sorted_KNN,75)+1.5*(prctile(sorted_KNN,75)-prctile(sorted_KNN,25));
    cutoff =prctile(sorted_KNN,97.5);
    %cutoff =0.8*cutoff;
    if PlotFigs == 1
        figure
        plot(sorted_KNN)
        hold on
        plot(ones(1,length(sorted_KNN)).*cutoff)
        title('Sorted Euclidean KNNs and Outlier Removal Cutoff','Fontsize',14)
    end
    
    points_to_keep_idx = find(D_KNN(:,end)<cutoff);
    outlier_idx = find(D_KNN(:,end)>=cutoff);
    X_before_denoising = X; %store in order to classify outliers later
    Labels_before_denoising = Labels;
    X = X(points_to_keep_idx,:);
    Labels = Labels(points_to_keep_idx,:);
    K=size(X,1);
    
end

if strcmp(DenoisingOpts.LocalAvgAllPts, 'yes')
    
    [IDX, D_KNN] = knnsearch(X,X,'k',DenoisingOpts.LocalAvgNbhdSize);
    LocalAverages = zeros(size(X));
    for i=1:size(X,1)
        LocalAverages(i,:) = mean(X(IDX(i,:),:));
    end
    X_before_local_avg = X;
    X = LocalAverages;
    
end

%% Save results

if strcmp(SaveResults,'yes')
    
    SaveProcessedResults
    
end
