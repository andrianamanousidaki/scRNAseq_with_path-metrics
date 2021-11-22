[IDX2, D_KNN2] = knnsearch(X_before_denoising,X_before_denoising,'k',min(AdjOpts.NumberOfNN,size(X_before_denoising,1)));
n2 = size(X_before_denoising,1);
Base=ones(n2,min(AdjOpts.NumberOfNN,size(X_before_denoising,1)));

for i=1:n2
    Base(i,:)=i*Base(i,:);
end

D_KNN2=D_KNN2';
IDX2=IDX2';
Base=Base';

A2=sparse(Base(:),IDX2(:),D_KNN2(:),n2,n2);
A2 = max(A2, A2'); 

p=AdjOpts.p;
A2p = A2.^p; %weight edges by ||x_i-x_j||^p

% If kNN graph is disconnected, add small edges to obtain a connected graph
A2p = ConnectGraph(A2p, X_before_denoising, p);

% Compute l_p distances from kNN graph
[PD_Dis2] = graphallshortestpaths(A2p);
PD_Dis2 = PD_Dis2.^(1/p); 

%%

% Extend Labels
pm_labels_extended = zeros(size(X_before_denoising,1),1);
pm_labels_extended(points_to_keep_idx) = pm_labels;
for i=1:length(outlier_idx)
     [D, I] = sort( PD_Dis2(outlier_idx(i),points_to_keep_idx) );
    pm_labels_extended(outlier_idx(i)) = pm_labels(I(1));
end