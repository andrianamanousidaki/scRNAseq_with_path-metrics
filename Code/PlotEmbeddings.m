%% Plot PCA, path metric, umap, and t-sne coordinates
%% Make 2d plots

figure
if MakeIndividualPlots == 0
    subplot(2,2,1)
end
grid on
hold on
axis equal
for q = 1:length(uniqueGroups)
      % Get indices of this particular unique group:
      ind = Labels==uniqueGroups(q); 
      % Plot only this group: 
      plot(X(ind,1),X(ind,2),'.','markersize',10,'Color',Colors(q,:)); 
end
axis square
if strcmp(DataSet, 'CellMix') || strcmp(DataSet, 'TMLung') || strcmp(DataSet, 'TMPanc')
    legend(LabelNames,'FontSize',10)
end
if ShowTitles == 1
    title('PCA Coordinates','FontSize',14)
end
%%
if MakeIndividualPlots == 0
    subplot(2,2,2)
else
    figure
end
grid on
hold on
axis equal
for q = 1:length(uniqueGroups)
      % Get indices of this particular unique group:
      ind = Labels==uniqueGroups(q); 
       %ind = pm_labels==uniqueGroups(q); 
      % Plot only this group: 
      plot(U(ind,1),U(ind,2),'.','markersize',10,'Color',Colors(q,:)); 
      %xlabel('PC_1','fontsize',14)
      %ylabel('PC_2','fontsize',14)
end
axis square
if  strcmp(DataSet, 'CellMix') || strcmp(DataSet, 'TMLung') || strcmp(DataSet, 'TMPanc')
    legend(LabelNames,'FontSize',10)
end
if ShowTitles == 1
    title('Path Metric + MDS Coordinates','FontSize',14)
end
%%
if MakeIndividualPlots == 0
    subplot(2,2,3)
else
    figure
end
grid on
hold on
axis equal
for q = 1:length(uniqueGroups)
      % Get indices of this particular unique group:
      ind = Labels==uniqueGroups(q); 
      % Plot only this group: 
      plot(umap_coordinates(ind,1),umap_coordinates(ind,2),'.','markersize',10,'Color',Colors(q,:)); 
end
axis square
%gscatter(umap_coordinates(:,1),umap_coordinates(:,2),Labels,'bgmkyrc','',[20 20 20 20],'off')
if  strcmp(DataSet, 'CellMix') || strcmp(DataSet, 'TMLung') || strcmp(DataSet, 'TMPanc')
    legend(LabelNames,'FontSize',10)
end
if ShowTitles == 1
    title('UMAP Coordinates','FontSize',14)
end
%%
if MakeIndividualPlots == 0
    subplot(2,2,4)
else
    figure
end
grid on
hold on
axis equal
for q = 1:length(uniqueGroups)
      % Get indices of this particular unique group:
      ind = Labels==uniqueGroups(q); 
      % Plot only this group: 
      plot(tsne_coordinates(ind,1),tsne_coordinates(ind,2),'.','markersize',10,'Color',Colors(q,:)); 
end
axis square
if  strcmp(DataSet, 'CellMix') || strcmp(DataSet, 'TMLung') || strcmp(DataSet, 'TMPanc')
    legend(LabelNames,'FontSize',10)
end
if ShowTitles == 1
    title('t-sne Coordinates','FontSize',14)
end
    
    