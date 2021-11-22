%% u-map coordinates colored by various clustering alg's

figure
subplot(2,2,1)
grid on
hold on
axis equal
for q = 1:length(uniqueGroups)
      % Get indices of this particular unique group:
      ind = Labels==uniqueGroups(q); 
      % Plot only this group: 
      plot(umap_coordinates(ind,1),umap_coordinates(ind,2),'.','markersize',10,'Color',Colors(q,:)); 
end
title('Ground Truth','FontSize',14)

subplot(2,2,2)
grid on
hold on
axis equal
for q = 1:length(uniqueGroups)
      % Get indices of this particular unique group:
      ind = pm_labels==uniqueGroups(q); 
      % Plot only this group: 
      plot(umap_coordinates(ind,1),umap_coordinates(ind,2),'.','markersize',10,'Color',Colors(q,:)); 
end
title('Path Metric','FontSize',14)

subplot(2,2,3)
grid on
hold on
axis equal
for q = 1:length(uniqueGroups)
      % Get indices of this particular unique group:
      ind = dbscan_labels==uniqueGroups(q); 
      % Plot only this group: 
      plot(umap_coordinates(ind,1),umap_coordinates(ind,2),'.','markersize',10,'Color',Colors(q,:)); 
end
title('DBSCAN','FontSize',14)

subplot(2,2,4)
grid on
hold on
axis equal
for q = 1:length(uniqueGroups)
      % Get indices of this particular unique group:
      ind = kmeans_labels==uniqueGroups(q); 
      % Plot only this group: 
      plot(umap_coordinates(ind,1),umap_coordinates(ind,2),'.','markersize',10,'Color',Colors(q,:)); 
end
title('k-means','FontSize',14)