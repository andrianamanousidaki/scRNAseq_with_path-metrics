%% Choose Data Set (or run ProcessData to create in workspace):

clear all

%%%%%%%%% RNA DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DataSet = 'RNAMix1';
%DataSet = 'RNAMix2';
%DataSet = 'TMLung';
%DataSet = 'Beta2';
%DataSet = 'TMPanc';
%DataSet = 'BaronPanc';
%DataSet = 'PBMC3k';
%DataSet = 'PBMC3k_SingleR';
%DataSet='PBMC4k_main'
%DataSet='PBMC4k_extented'
DataSet ='CellMixSng';

% Defunct:
%DataSet = 'CellMix';
%DataSet = 'Beta';

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
    

%% Option 1 for RNA (loads all normalizations, basic/SCT and Linnorm, with and without denoising...): 

addpath(genpath('../DataSets/ProcessedDataSets'))
LoadProcessedData; DenoisingOpts.RemoveOutliers = 'no';

%% Option 2 for RNA (create 1 normalization of choice; set parameters in ProcessData as desired):

%addpath(genpath('../DataSets'))
%ProcessData

%% Option 3 for Manifold(just load anything):

%addpath(genpath('../DataSets'))
%LoadManifoldData

%% Path metric options

AdjOpts.p = 4;%AdjOpts.p = 4;%for PBMC4K
AdjOpts.NumberOfNN = min(size(X,1),500); %Parameter defining the kNN graph in which path metrics are computed; When NumberOfNN = size(Data,1), the computed path metrics are exact. Otherwise they maybe approximations. NumberOfNN should scale at least like log(sample size). 
Silhouette.Criterion='yes';% If 'yes' then prediction of optimal number of clusters with silhouette(silhouette([],clustering_labels,squareform(PD_Dis)))
run.adjusted.kmeans='yes';% If 'yes' tiny clusters are merged to the closest cluster
klist=20; % Find silhouette scores for number of clusters in 1:klist
SpectralOpts.EmbeddingDimension = 40;%RNAMIX spec.emb.dim=2
SpectralOpts.LearnMDSEmbeddingDimension = 'yes'; %options 'yes' or 'no'
RunComparisons = 'no'; %Options: 'yes' or 'no'; compare clustering and geometric fidelity with other methods
RunCPM = 'no';
PlotFigs = 1;
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
                
if strcmp(Silhouette.Criterion,'yes')
    
    sil_scores =zeros(1,klist);
    matrix_pm_labels= zeros(n,klist);

    for k= 1:klist;
        if strcmp(run.adjusted.kmeans,'yes')
             matrix_pm_labels(:,k) = adjusted_kmeans(U,k);
        else
            matrix_pm_labels(:,k) = kmeans(U,k,'replicates',20);
        end
        sil_scores(1,k)= mean(silhouette([],matrix_pm_labels(:,k),pdist(U)));
    end   
    figure
    plot(1:klist,sil_scores)
    optimal_num_clusters= find(sil_scores == max(sil_scores))
    pm_labels=matrix_pm_labels(:,optimal_num_clusters);
    PM_evaluation=clustering_evaluation(pm_labels, Labels)

elseif strcmp(Silhouette.Criterion,'no')
    
    k =length(unique(Labels)); %number of clusters
    
    if strcmp(RunComparisons, 'yes')
        
        
        %% k-mean clustering (basic scaling)
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
        [length(unique(Labels)) length(unique(dbscan_labels))]
        [dbscan_OA, ~, ~, ~, dbscan_labels] = GetAccuracies(dbscan_labels,Labels,k); 
        dbscan_evaluation = clustering_evaluation(dbscan_labels,Labels);
        
        %% Path metric clustering (basic/SCT scaling)
        if strcmp(run.adjusted.kmeans,'yes')
            pm_labels = adjusted_kmeans(U,k);
        else
            pm_labels = kmeans(U,k,'replicates',20);
        end    
        %pm_labels = kmeans(U1(:,1:5),k,'replicates',20);
        [pm_OA, ~, ~, ~, pm_labels] = GetAccuracies(pm_labels,Labels,k);
        PM_evaluation = clustering_evaluation(pm_labels, Labels);
        
        if strcmp(DataSetType,'RNA')

            %% k-mean clustering (linnorm scaling)
            if strcmp(run.adjusted.kmeans,'yes')
              kmeans_labels_lin = adjusted_kmeans(XLin,k);
            else
              kmeans_labels_lin = kmeans(XLin,k,'replicates',20);
            end    
            [kmeans_OA_lin, ~, ~, ~, kmeans_labels_lin] = GetAccuracies(kmeans_labels_lin,Labels,k); 
            kmeans_evaluation_lin = clustering_evaluation(kmeans_labels_lin,Labels);

            %% DBSCAN clustering (linnorm scaling)
            [dbscan_labels_lin,~] = RunDBSCAN(XLin,epsilon_lin,minpts);
            [dbscan_OA_lin, ~, ~, ~, dbscan_labels_lin] = GetAccuracies(dbscan_labels_lin,Labels,k); 
            dbscan_evaluation_lin = clustering_evaluation(dbscan_labels_lin,Labels);

            %% UMAP + clustering:
            rng(7); %for reproducibility, since output is random
            [umap_coordinates, ~, ~, ~]=run_umap(XLinND);%change to X, XLin, XND
            %[umap_coordinates, ~, ~, ~]=run_umap(XLinND,'n_components',ratio_dim_est); %keep same
            %number of dimensions
            umap_labels = RunDBSCAN(umap_coordinates,epsilon_umap,minpts);
            [umap_OA, ~, ~, ~, umap_labels] = GetAccuracies(umap_labels,Labels,k); 
            umap_evaluation = clustering_evaluation(umap_labels,Labels)

            %% t-sne + clustering:
            %tsne_coordinates = tsne(XLinND); %change to X, XLin, XND
            tsne_coordinates = tsne(XLinND,[],ratio_dim_est);% keep
            %same number of dimensions
            if strcmp(run.adjusted.kmeans,'yes')
                tsne_labels = adjusted_kmeans(tsne_coordinates,k);
            else
                tsne_labels = kmeans(tsne_coordinates,k,'replicates',20);
            end
            [tsne_OA, ~, ~, ~, tsne_labels] = GetAccuracies(tsne_labels,Labels,k); 
            tsne_evaluation = clustering_evaluation(tsne_labels,Labels);

            Method = {'kmeans','kmeans_lin','dbscan','dbscan_lin','umap+dbscan','tsne+kmeans','pm'}; Method=Method';
            NumClusters = [length(unique(kmeans_labels)); length(unique(kmeans_labels_lin)); length(unique(dbscan_labels)); length(unique(dbscan_labels_lin)); length(unique(umap_labels)); length(unique(tsne_labels)); length(unique(pm_labels))];
            ARI = [kmeans_evaluation.ARI; kmeans_evaluation_lin.ARI; dbscan_evaluation.ARI; dbscan_evaluation_lin.ARI; umap_evaluation.ARI; tsne_evaluation.ARI; PM_evaluation.ARI];
            ECP = [kmeans_evaluation.ECP; kmeans_evaluation_lin.ECP; dbscan_evaluation.ECP; dbscan_evaluation_lin.ECP; umap_evaluation.ECP; tsne_evaluation.ECP; PM_evaluation.ECP];
            ECA = [kmeans_evaluation.ECA; kmeans_evaluation_lin.ECA; dbscan_evaluation.ECA; dbscan_evaluation_lin.ECA; umap_evaluation.ECA; tsne_evaluation.ECA; PM_evaluation.ECA];
            clustering_results = table(Method,NumClusters,ARI,ECP,ECA)
            
        elseif strcmp(DataSetType,'Manifold')
            
            
            Method = {'kmeans','dbscan','pm'}; Method=Method';
            NumClusters = [length(unique(kmeans_labels)); length(unique(dbscan_labels)); length(unique(pm_labels))];
            ARI = [kmeans_evaluation.ARI; dbscan_evaluation.ARI; PM_evaluation.ARI];
            ECP = [kmeans_evaluation.ECP; dbscan_evaluation.ECP; PM_evaluation.ECP];
            ECA = [kmeans_evaluation.ECA; dbscan_evaluation.ECA; PM_evaluation.ECA];
            clustering_results = table(Method,NumClusters,ARI,ECP,ECA)
            
        end
            
        
    elseif strcmp(RunComparisons, 'no')
        
        if strcmp(DenoisingOpts.RemoveOutliers,'no')
        
            % Path metric clustering (basic scaling)
            if strcmp(run.adjusted.kmeans,'yes')
                pm_labels = adjusted_kmeans(U,k);
            else    
                pm_labels = kmeans(U,k,'replicates',20);
            end
            [pm_OA, ~, ~, ~, pm_labels] = GetAccuracies(pm_labels,Labels,k);
            PM_evaluation = clustering_evaluation(pm_labels, Labels)
        
        elseif strcmp(DenoisingOpts.RemoveOutliers,'yes')
            if strcmp(run.adjusted.kmeans,'yes')
                pm_labels = adjusted_kmeans(U,k);
            else
                pm_labels = kmeans(U,k,'replicates',20);
            end
            PM_evaluation_core = clustering_evaluation(pm_labels, Labels)
            
            % Extend Labels
            
            ExtendLabels
            
%             pm_labels_extended = zeros(size(X_before_denoising,1),1);
%             pm_labels_extended(points_to_keep_idx) = pm_labels;
%             [NN1, D_NN1] = knnsearch(X,X_before_denoising(outlier_idx,:),'k',1);
%             for i=1:length(outlier_idx)
%                 pm_labels_extended(outlier_idx(i)) = mode(pm_labels(NN1(i,:)));
%             end
            
            PM_evaluation=clustering_evaluation(categorical(pm_labels_extended), categorical(Labels_before_denoising))
            
        end
    
    end

end

if strcmp(RunComparisons,'yes')

    %% Compute geometric perturbation
    
    [pca_geo_pert, pca_EmbeddingCentroids] = GeometricPerturbation(X,X,Labels);
    [pm_geo_pert, pm_EmbeddingCentroids] = GeometricPerturbation(X,U,Labels);
    [pm_geo_pert_2d, pm_EmbeddingCentroids_2d] = GeometricPerturbation(X,U(:,1:2),Labels);
    
    if strcmp(DataSetType,'RNA')
    
        [umap_geo_pert, umap_EmbeddingCentroids] = GeometricPerturbation(XLinND,umap_coordinates,Labels);
        [tsne_geo_pert,tsne_EmbeddingCentroids] = GeometricPerturbation(XLinND,tsne_coordinates,Labels);
        
    elseif strcmp(DataSetType,'Manifold')
    
        %% UMAP:
        [umap_coordinates, ~, ~, ~]=run_umap(X);
        %[umap_coordinates, ~, ~, ~]=run_umap(X,'n_components',ratio_dim_est); %keep same
        %number of dimensions

        %% t-sne + clustering:
        tsne_coordinates = tsne(X); %change to X, XLin, XND
        %tsne_coordinates = tsne(X,[],ratio_dim_est);% keep
        %same number of dimensions
        
        [umap_geo_pert, umap_EmbeddingCentroids] = GeometricPerturbation(X,umap_coordinates,Labels);
        [tsne_geo_pert,tsne_EmbeddingCentroids] = GeometricPerturbation(X,tsne_coordinates,Labels);

    end

    PlotEmbeddings

    PlotClusterings

    if strcmp(RunCPM,'yes')
    
        ydata= cpm(U,2,0);

        figure
        grid on
        hold on
        axis equal
        for q = 1:length(uniqueGroups)
              % Get indices of this particular unique group:
              ind = Labels==uniqueGroups(q); 
              % Plot only this group: 
              plot(ydata(ind,1),ydata(ind,2),'.','markersize',10,'Color',Colors(q,:)); 
        end
        title('Path Metric + MDS + CPM','FontSize',14)

        [cpm_pm_geo_pert, cpm_pm_EmbeddingCentroids] = GeometricPerturbation(X,ydata,Labels);

        Method = {'pm','2d pm','umap','tsne','pm+cpm'};  Method=Method';
        GeoPert = [pm_geo_pert; pm_geo_pert_2d; umap_geo_pert; tsne_geo_pert; pm_cpm_geo_pert];
        geo_pert_results = table(Method,GeoPert)
        
    else
        
        Method = {'pm','2d pm','umap','tsne'};  Method=Method';
        GeoPert = [pm_geo_pert; pm_geo_pert_2d; umap_geo_pert; tsne_geo_pert];
        geo_pert_results = table(Method,GeoPert)
        
    end
        
    
    PlotTrees   
    
    
elseif strcmp(RunComparisons,'no') 
    
    [pca_geo_pert, pca_EmbeddingCentroids] = GeometricPerturbation(X,X,Labels);
    [pm_geo_pert, pm_EmbeddingCentroids] = GeometricPerturbation(X,U,Labels);
    [pm_geo_pert_2d, pm_EmbeddingCentroids_2d] = GeometricPerturbation(X,U(:,1:2),Labels);
    pm_geo_pert
    pm_geo_pert_2d
    
    PlotTrees 
    
end

%%
if (size(U,2)>2) && (PlotFigs==1)
    
    Make3dPlots 
    
end  
