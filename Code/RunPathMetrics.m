%% Choose Data Set:

clear all

%%%%%%%%% RNA DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DataSet = 'RNAMix1';
%DataSet = 'RNAMix2';
%DataSet = 'TMLung';
%DataSet = 'Beta2';
%DataSet = 'TMPanc';
%DataSet = 'BaronPanc';
%DataSet = 'PBMC4k_main'
%DataSet = 'CellMixSng';


%%%%%%%%% MANIFOLD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DataSet = 'Balls';
%DataSet = 'ElongatedWithBridge';
%DataSet = 'SwissRoll';
%DataSet = 'GL_Manifold';

if ((strcmp(DataSet,'Balls') || strcmp(DataSet,'ElongatedWithBridge')) || strcmp(DataSet,'SwissRoll')) || strcmp(DataSet,'GL_Manifold') == 1
    DataSetType = 'Manifold';
else
    DataSetType = 'RNA';
end

if strcmp(DataSet, 'CellMixSng')
    LabelNames = {'HCC827', 'H1975' , 'H838',  'H2228' , 'A549'};
elseif strcmp(DataSet, 'TMLung')
    LabelNames = {'lung', 'alveolar' , 'mast & im.',  'circ. mono.' , 'inv. mono.','multicil.','dendritic'};
elseif strcmp(DataSet, 'TMPanc')
    LabelNames={'acinar','beta','stellate','pancreatic A','ductal','pancreatic D','pancreatic PP'};
end

if strcmp(DataSetType,'RNA')
    
    addpath(genpath('ProcessedDataSets'))
    LoadProcessedData; DenoisingOpts.RemoveOutliers = 'no';
    
elseif strcmp(DataSetType,'Manifold')
    
    addpath(genpath('ManifoldDataSets'))
    LoadManifoldData;
    
end
    

%% Path metric options

AdjOpts.p =1.5;
AdjOpts.NumberOfNN = min(size(X,1),500); %Parameter defining the kNN graph in which path metrics are computed; When NumberOfNN = size(Data,1), the computed path metrics are exact. Otherwise they maybe approximations. NumberOfNN should scale at least like log(sample size). 
Silhouette.Criterion='no';% If 'yes' then prediction of optimal number of clusters with silhouette(silhouette([],clustering_labels,squareform(PD_Dis)))
run.adjusted.kmeans='yes';% If 'yes' tiny clusters are merged to the closest cluster
klist=20; % Find silhouette scores for number of clusters in 1:klist
SpectralOpts.EmbeddingDimension = 40;%RNAMIX spec.emb.dim=2
SpectralOpts.LearnMDSEmbeddingDimension = 'yes'; %options 'yes' or 'no'
RunComparisons = 'yes'; %Options: 'yes' or 'no'; compare clustering and geometric fidelity with other methods
Comparisons.HigherDimension = 'no'; %Options: 'yes' or 'no'; Default: 'no'; if yes, use PM dimension for umap and tsne; if no, use d=2 for umap and tsne. Note: parameters for clustering comparisons are chosen to yield right number of clusters for 'no' 
PlotFigs = 0;
MakeIndividualPlots = 0;
ShowTitles = 1;

%% Add to Matlab path:

addpath(genpath('cmccomb-rand_index-10a43e9'))
addpath(genpath('MunkresAlgorithm'))
addpath(genpath('umapFileExchange (1.5.2)'))
addpath(genpath('CPM-master'))

%% Compute shortest paths:

[IDX, D_KNN] = knnsearch(X,X,'k',min(AdjOpts.NumberOfNN,size(X,1)));
n = size(X,1);
Base=ones(n,min(AdjOpts.NumberOfNN,size(X,1)));

for i=1:n
    Base(i,:)=i*Base(i,:);
end

D_KNN=D_KNN';
IDX=IDX';
Base=Base';

A=sparse(Base(:),IDX(:),D_KNN(:),n,n);
A = max(A, A'); 

p=AdjOpts.p;
Ap = A.^p; %weight edges by ||x_i-x_j||^p

% If kNN graph is disconnected, add small edges to obtain a connected graph
Ap = ConnectGraph(Ap, X, p);

% Compute l_p distances from kNN graph
[PD_Dis] = graphallshortestpaths(Ap);
PD_Dis = PD_Dis.^(1/p);    

%% Compute path metric MDS embedding:

[U1, Eigvals]=cmdscale(PD_Dis,SpectralOpts.EmbeddingDimension);

% Select MDS embedding dimension: 

if strcmp(SpectralOpts.LearnMDSEmbeddingDimension, 'yes')
    
    ratios=Eigvals(1:SpectralOpts.EmbeddingDimension-1)./Eigvals(2:end);
    max_dim = max(find(Eigvals./max(Eigvals) > .01))-1; %want BOTH eigvals to not be too small
    ratio_dim_est = find(max(ratios(3:max_dim))==ratios(3:max_dim))+2; %at least 3 dimensions
    %ratio_dim_est = find(max(ratios(2:max_dim))==ratios(2:max_dim))+1; %at least 2 dimensions
    U = U1(:,1:ratio_dim_est);
    
end

%% Run clustering:

uniqueGroups = unique(Labels);
Colors=linspecer(length(uniqueGroups));
                
k =length(unique(Labels)); %number of clusters

%% Path metric clustering (basic/SCT scaling)

rng(6); %for reproducibility, since output of k-means is random
if strcmp(run.adjusted.kmeans,'yes')
    pm_labels = adjusted_kmeans(U,k);
else
    pm_labels = kmeans(U,k,'replicates',20);
end    
%pm_labels = kmeans(U1(:,1:5),k,'replicates',20);
[pm_OA, ~, ~, ~, pm_labels] = GetAccuracies(pm_labels,Labels,k);
PM_evaluation = clustering_evaluation(pm_labels, Labels);

if strcmp(RunComparisons, 'yes')


    %% k-mean clustering (basic scaling)

    rng(5); %for reproducibility, since output of k-means is random
    if strcmp(run.adjusted.kmeans,'yes')
        kmeans_labels = adjusted_kmeans(X,k);
    else
        kmeans_labels = kmeans(X,k,'replicates',20);
    end
    [kmeans_OA, ~, ~, ~, kmeans_labels] = GetAccuracies(kmeans_labels,Labels,k); 
    kmeans_evaluation = clustering_evaluation(kmeans_labels,Labels);

    %% DBSCAN clustering (basic/SCT scaling)
    minpts = 10;
    [dbscan_labels,~] = RunDBSCAN(X,epsilon,minpts);
    [dbscan_OA, ~, ~, ~, dbscan_labels] = GetAccuracies(dbscan_labels,Labels,k); 
    dbscan_evaluation = clustering_evaluation(dbscan_labels,Labels);

    %% UMAP + clustering:
    rng(7); %for reproducibility, since output is random
    [umap_coordinates, ~, ~, ~]=run_umap(XLinND);%change to X, XLin, XND
    if strcmp(Comparisons.HigherDimension,'yes')
        [umap_coordinates_rd, ~, ~, ~]=run_umap(XLinND,'n_components',ratio_dim_est); %keep same number of dimensions
    end
    umap_labels = RunDBSCAN(umap_coordinates,epsilon_umap,minpts);
    [umap_OA, ~, ~, ~, umap_labels] = GetAccuracies(umap_labels,Labels,k); 
    umap_evaluation = clustering_evaluation(umap_labels,Labels)

    %% t-sne + clustering:
    rng(8); %for reproducibility, since output is random
    tsne_coordinates = tsne(XLinND); %change to X, XLin, XND
    if strcmp(Comparisons.HigherDimension,'yes') && (size(XLinND,2) >= ratio_dim_est)
        tsne_coordinates_rd = tsne(XLinND,'NumDimensions',ratio_dim_est);% keep same number of dimensions
    end
    if strcmp(run.adjusted.kmeans,'yes')
        tsne_labels = adjusted_kmeans(tsne_coordinates,k);
    else
        tsne_labels = kmeans(tsne_coordinates,k,'replicates',20);
    end
    [tsne_OA, ~, ~, ~, tsne_labels] = GetAccuracies(tsne_labels,Labels,k); 
    tsne_evaluation = clustering_evaluation(tsne_labels,Labels);
    
    %% Make clustering comparison table

    Method = {'kmeans','dbscan','umap+dbscan','tsne+kmeans',strcat('pm', num2str(p))}; Method=Method';
    NumClusters = [length(unique(kmeans_labels)); length(unique(dbscan_labels)); length(unique(umap_labels)); length(unique(tsne_labels)); length(unique(pm_labels))];
    ARI = [kmeans_evaluation.ARI; dbscan_evaluation.ARI; umap_evaluation.ARI; tsne_evaluation.ARI; PM_evaluation.ARI]; ARI = round(ARI,3);
    ECP = [kmeans_evaluation.ECP; dbscan_evaluation.ECP; umap_evaluation.ECP; tsne_evaluation.ECP; PM_evaluation.ECP]; ECP = round(ECP,3);
    ECA = [kmeans_evaluation.ECA; dbscan_evaluation.ECA; umap_evaluation.ECA; tsne_evaluation.ECA; PM_evaluation.ECA]; ECA = round(ECA,3);
    clustering_results = table(Method,NumClusters,ARI,ECP,ECA)


elseif strcmp(RunComparisons, 'no')

    Method = {strcat('pm', num2str(p))};
    ARI = round(PM_evaluation.ARI,3);
    ECA = round(PM_evaluation.ECA,3);
    ECP = round(PM_evaluation.ECP,3);
    NumClusters = length(unique(pm_labels));
    clustering_results = table(Method,NumClusters,ARI,ECP,ECA)

end

%% Compute geometric perturbation

if strcmp(RunComparisons,'yes')
    
    [pca_geo_pert, pca_EmbeddingCentroids] = GeometricPerturbation(X,X,Labels);
    [pm_geo_pert, pm_EmbeddingCentroids] = GeometricPerturbation(X,U,Labels);
    [pm_geo_pert_2d, pm_EmbeddingCentroids_2d] = GeometricPerturbation(X,U(:,1:2),Labels);
    [umap_geo_pert, umap_EmbeddingCentroids] = GeometricPerturbation(XLinND,umap_coordinates,Labels);
    [tsne_geo_pert,tsne_EmbeddingCentroids] = GeometricPerturbation(XLinND,tsne_coordinates,Labels);
    if strcmp(Comparisons.HigherDimension,'yes')
        [umap_rd_geo_pert, umap_rd_EmbeddingCentroids] = GeometricPerturbation(XLinND,umap_coordinates_rd,Labels);
        if size(XLinND,2) >= ratio_dim_est
            [tsne_rd_geo_pert, tsne_rd_EmbeddingCentroids] = GeometricPerturbation(XLinND,tsne_coordinates_rd,Labels);
        else
            tsne_rd_geo_pert = 1000; %Note: geometric perturbation of 1000 means 'NA'
        end
    end
        
    
    if PlotFigs==1

        PlotEmbeddings

        PlotClusterings
        
    end
%%  
    if strcmp(Comparisons.HigherDimension,'no')
        Method = {'2d umap','2d tsne',strcat('2d pm', num2str(p)),strcat('rd pm', num2str(p))};  Method=Method';
        GeoPert = [umap_geo_pert; tsne_geo_pert; pm_geo_pert_2d; pm_geo_pert]; GeoPert=round(GeoPert,3);
        geo_pert_results = table(Method,GeoPert)
    elseif strcmp(Comparisons.HigherDimension,'yes')
         Method = {'2d umap','rd umap','2d tsne','rd tsne',strcat('2d pm', num2str(p)),strcat('rd pm', num2str(p))};  Method=Method';
         GeoPert = [umap_geo_pert; umap_rd_geo_pert; tsne_geo_pert; tsne_rd_geo_pert; pm_geo_pert_2d; pm_geo_pert]; GeoPert=round(GeoPert,3);
         geo_pert_results = table(Method,GeoPert)
    end

    if PlotFigs==1
    
        PlotTrees 
        
    end
    
    
elseif strcmp(RunComparisons,'no') 
    
    [pca_geo_pert, pca_EmbeddingCentroids] = GeometricPerturbation(X,X,Labels);
    [pm_geo_pert, pm_EmbeddingCentroids] = GeometricPerturbation(X,U,Labels);
    [pm_geo_pert_2d, pm_EmbeddingCentroids_2d] = GeometricPerturbation(X,U(:,1:2),Labels);
     Method = {strcat('2d pm', num2str(p)),strcat('rd pm', num2str(p))};  Method=Method';
     GeoPert = [pm_geo_pert_2d; pm_geo_pert]; GeoPert=round(GeoPert,3);
     geo_pert_results = table(Method,GeoPert)
    
    if PlotFigs==1
    
        PlotTrees 
        
    end 
    
end

%%
if (size(U,2)>2) && (PlotFigs==1)
    
    Make3dPlots 
    
end  
