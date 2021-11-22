DenoisingOpts.LocalAvgAllPts = 'yes'; %Options: 'yes' or 'no'
DenoisingOpts.LocalAvgNbhdSize = 12;
DenoisingOpts.RemoveOutliers='no';

%% Load Data (reduce with PCA if high-dimensional):

if strcmp(DataSet,'ElongatedWithBridge')
    
    load('ElongatedGaussiansWithBridge3.mat'); %loads X, Labels
    epsilon = 0.2; epsilon_umap = 0.8; %values for denoised data
    
elseif strcmp(DataSet, 'Balls')
    
    load('Balls_n400_r53.mat')  %loads X, Labels
    epsilon = 0.08687; epsilon_umap = 0.64;
    
elseif strcmp(DataSet,'SwissRoll')
    
    load('SwissRoll1.mat'); %loads Data, Labels
    X = Data;
    epsilon = 3; epsilon_umap = 5; %values for denoised data
    
elseif strcmp(DataSet,'GL_Manifold')
    
    load('GLmanifold_d9_N3000_k3_sig0075.mat'); %loads Data, Labels
    X = Data;
    [coeff2,score2,latent2] = pca(X,'NumComponents', 40 );
    X = score2;
    epsilon = .15; epsilon_umap = 3; %values for denoised data
    
end

%% Denoise Data:

if strcmp(DenoisingOpts.LocalAvgAllPts, 'yes')
    
    [IDX, D_KNN] = knnsearch(X,X,'k',DenoisingOpts.LocalAvgNbhdSize);
    LocalAverages = zeros(size(X));
    for i=1:size(X,1)
        LocalAverages(i,:) = mean(X(IDX(i,:),:));
    end
    XND = X;
    X = LocalAverages;
    
end

%% For compatibility with RunPathMetrics, define XLin, XLin_ND

XLin = X;
if strcmp(DenoisingOpts.LocalAvgAllPts, 'yes')
    XLinND = XND;
end
