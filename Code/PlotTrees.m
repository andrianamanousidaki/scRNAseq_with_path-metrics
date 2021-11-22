Z_pca = linkage(pca_EmbeddingCentroids,'average');
Z_pm = linkage(pm_EmbeddingCentroids,'average');
if strcmp(RunComparisons,'yes') && strcmp(DataSetType,'RNA') == 1
    Z_umap = linkage(umap_EmbeddingCentroids,'average');
    Z_tsne = linkage(tsne_EmbeddingCentroids,'average');
end
%         Z_pca = linkage(pca_EmbeddingCentroids,'single');
%         Z_pm = linkage(pm_EmbeddingCentroids,'single');
%         Z_umap = linkage(umap_EmbeddingCentroids,'single');
%         Z_tsne = linkage(tsne_EmbeddingCentroids,'single');
%    figure
if strcmp(RunComparisons,'yes') && strcmp(DataSetType,'RNA') == 1
    
    if strcmp(DataSet, 'CellMixSng')
        figure
        if MakeIndividualPlots == 0
            subplot(2,2,1)
        end
        dendrogram(Z_pca,'Labels',LabelNames)
        axis square
        set(gca,'YTick',[],'box','on')
        if ShowTitles == 1
            title('PCA (40 dims)','Fontsize',14)
        end
        if MakeIndividualPlots == 0
            subplot(2,2,2)
        else
            figure
        end
        dendrogram(Z_pm,'Labels',LabelNames)
        axis square
        set(gca,'YTick',[],'box','on')
        if ShowTitles == 1
            title('Path Metrics','Fontsize',14)
        end
        if MakeIndividualPlots == 0
            subplot(2,2,3)
        else
            figure
        end
        dendrogram(Z_umap,'Labels',LabelNames)
        axis square
        set(gca,'YTick', [],'box','on')
        if ShowTitles == 1
            title('Umap','Fontsize',14)
        end
        if MakeIndividualPlots == 0
            subplot(2,2,4)
        else
            figure
        end
        dendrogram(Z_tsne,'Labels',LabelNames)
        axis square
        set(gca,'YTick', [],'box','on')
        if ShowTitles == 1
            title('Tsne','Fontsize',14)
        end
    else
        subplot(2,2,1)
        dendrogram(Z_pca)
        set(gca,'YTick', [],'box','on')
        title('Dendrogram of Cluster Means: PCA','Fontsize',14)
        subplot(2,2,2)
        dendrogram(Z_pm)
        set(gca,'YTick', [],'box','on')
        title('Dendrogram of Cluster Means: Path Metrics','Fontsize',14)
        subplot(2,2,3)
        dendrogram(Z_umap)
        set(gca,'YTick', [],'box','on')
        title('Dendrogram of Cluster Means: Umap','Fontsize',14)
        subplot(2,2,4)
        dendrogram(Z_tsne)
        set(gca,'YTick', [],'box','on')
        title('Dendrogram of Cluster Means: Tsne','Fontsize',14)
    end
    
else
    
    if strcmp(DataSet, 'CellMixSng')
        figure
        if MakeIndividualPlots == 0
            subplot(1,2,1)
        end
        dendrogram(Z_pca,'Labels',LabelNames)
        axis square
        set(gca,'YTick',[],'box','on')
        if ShowTitles == 1
            title('PCA (40 dims)','Fontsize',14)
        end
        if MakeIndividualPlots == 0
            subplot(1,2,2)
        else
            figure
        end
        dendrogram(Z_pm,'Labels',LabelNames)
        axis square
        set(gca,'YTick',[],'box','on')
        if ShowTitles == 1
            title('Path Metrics','Fontsize',14)
        end
    else
        figure
        subplot(1,2,1)
        dendrogram(Z_pca)
        set(gca,'YTick', [],'box','on')
        title('Dendrogram of Cluster Means: PCA','Fontsize',14)
        subplot(1,2,2)
        dendrogram(Z_pm)
        set(gca,'YTick', [],'box','on')
        title('Dendrogram of Cluster Means: Path Metrics','Fontsize',14)
    end    
    
end