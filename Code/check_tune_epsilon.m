%% Check epsilon for DBSCAN clustering (basic/SCT scaling)

% minpts = 10;
% [dbscan_labels,~] = RunDBSCAN(X,epsilon,minpts);
% [length(unique(Labels)) length(unique(dbscan_labels))]

% %% Check epsilon_lin for DBSCAN clustering (Linnorm scaling)
% 
% minpts = 10;
% [dbscan_labels_lin,~] = RunDBSCAN(XLin,epsilon_lin,minpts);
% [length(unique(Labels)) length(unique(dbscan_labels_lin))]

%% Check epsilon_umap for DBSCAN clustering (Linnorm scaling)

minpts = 10;
%[umap_coordinates, ~, ~, ~]=run_umap(XLinND); %change to X, XLin, XND
%umap_labels = kmeans(umap_coordinates,k,'replicates',20);
umap_labels = RunDBSCAN(umap_coordinates,epsilon_umap,minpts);
[length(unique(Labels)) length(unique(umap_labels))]
    
%% For tuning the epsilon

data_input = umap_coordinates;

typ_dis = prctile(pdist(data_input),10)

%eps_values = linspace(typ_dis/4, typ_dis, 20);
eps_values = linspace(.2, 10, 20);
%eps_values = linspace(5.5, 6, 20);
%eps_values = linspace(.0865, .0875, 20);
NumClusters = zeros(size(eps_values));
for w = 1:length(eps_values)
    epsilon_temp = eps_values(w);
    [data_input_labels,~] = RunDBSCAN(data_input,epsilon_temp,minpts);
    [length(unique(Labels)) length(unique(data_input_labels))];
    NumClusters(w) = length(unique(data_input_labels));
end
figure
plot(eps_values, NumClusters)
hold on
plot(eps_values, length(unique(Labels))*ones(size(NumClusters)))
max(eps_values(find(NumClusters==length(unique(Labels)))))  