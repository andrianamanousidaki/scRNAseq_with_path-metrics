%% Make 3d plots

figure
% PCA Plot
if MakeIndividualPlots == 0
    subplot(1,2,1)
end
if size(X,2)==2
    grid on
    hold on
    for q = 1:length(uniqueGroups)
          % Get indices of this particular unique group:
          ind = Labels==uniqueGroups(q); 
          % Plot only this group: 
          plot(X(ind,1),X(ind,2),'.','markersize',10,'Color',Colors(q,:)); 
    end
    axis equal
  %  axis square
    if strcmp(DataSet, 'CellMix')
        legend(LabelNames,'FontSize',12)
    end
    if ShowTitles == 1
        title('PCA Coordinates','FontSize',14)
    end
else
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
          plot3(X(ind,1),X(ind,2),X(ind,3),'.','markersize',10,'Color',Colors(q,:)); 
    end
    axis square
    if strcmp(DataSet, 'CellMix')
        legend(LabelNames,'FontSize',12)
    end
    if ShowTitles == 1
        title('PCA embedding of processed data','fontsize',14)
    end
    xlabel('PC_1','fontsize',14)
    ylabel('PC_2','fontsize',14)
    zlabel('PC_3','fontsize',14)
end

%% PM Plot
if MakeIndividualPlots == 0
    subplot(1,2,2)
else
    figure
end
view(3)
grid on
hold on
axis equal
rotate3d on
% Plot each group individually: 
for q = 1:length(uniqueGroups)
      % Get indices of this particular unique group:
      ind = Labels==uniqueGroups(q); 
      %ind = pm_labels==uniqueGroups(q); 
      % Plot only this group: 
      plot3(U(ind,1),U(ind,2),U(ind,3),'.','markersize',10,'Color',Colors(q,:)); 
end
axis square
if strcmp(DataSet, 'CellMix')
    legend(LabelNames,'FontSize',12)
end
if ShowTitles == 1
    title('Path Metric embedding of processed data','fontsize',14)
end
xlabel('PM_1','fontsize',14)
ylabel('PM_2','fontsize',14)
zlabel('PM_3','fontsize',14)