# scRNAseq clustering and visualization with path-metrics
Here we provide source code for a clustering and visualization method for single-cell RNA-seq data sets, that exploits **path metrics**. Using this method, distances between cells are measured in a data-driven way which is both density
sensitive (decreasing distances across high density regions) and respects the underlying data geometry.
By combining path metrics with multidimensional scaling, a low dimensional embedding of the data is
obtained which respects both the global geometry of the data and preserves cluster structure.

# Requirements 

1. [Adjusted Rand index Matlab function version 1.0.0.0](https://www.mathworks.com/matlabcentral/fileexchange/49908-adjusted-rand-index)
2. [Munkres Assignment Algorithm Matlab function](https://www.mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm)
3. [Uniform Manifold Approximation and Projection (UMAP) Matlab function >= version 1.5.2](https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)

# Description of files
* To reproduce results of processed data sets:
  1.	Use the link in Processed data sets folder to access download data sets.
  2.	Run *RunPathMetrics.m* using Matlab version 2021a.


* To produce the simulated beta cells data set:
  1.	Download "GSM2230757_human1_umifm_counts.csv.gz" from [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2230757)
  2.	Run *Simulationscript.R*
  
  
* To process the data sets used in results:
  1.	Download the real data sets using the links mentioned on “processing and imputation of data.R”
  2.	Run *processing and imputation of data.R*
  3.  Scale data
   *  **For Basic scaling:** input in Matlab the csv file of the data set of interest produced in step 2, along with the true labels file. Then run *Code/ProcessData.m* after uncommenting the name of the data set and the line “Normalization = 'Basic';”.
   *  **For Linnorm scaling:** run *Linnorm transformation.R* for the data set of interest produced in a. Then input in Matlab the output csv file along with the true labels file. Then run *Code/ProcessData.m* after uncommenting the name of the data set and the line “Normalization = 'Linnorm';”. Now you can use the output to run *Code and data sets/Code_for_github/RunPathMetrics.m*.
   *  **For SCT scaling:** run *sctranform_for_PM_after_SAVER_imputation.R* for the data set of interest produced in a. Then input in Matlab the output csv file along with the true labels file. Then run *Code/ProcessData.m* after uncommenting the name of the data set and the line “Normalization = 'SCT';”. Now you can use the output to run *Code/RunPathMetrics.m*.


* Other files description:
  1.	*Seurat_clustering_analysis.R :* Apply Seurat clustering on data set produced by *processing and imputation of data.R*
  2.	*Seurat_clustering_analysis_default_res:* Apply Seurat clustering with resolution 0.8 on data set produced by *processing and imputation of data.R*
  3.	*RunSNNClustering.R:* Apply SNN clustering on data sets after Basic or SCT scaling (created as described in 3b.)
  4.	*clustering_evaluation.R:* Function that produces ARI, Entropy of cluster accuracy and Entropy of cluster Purity for a clustering.
  5.	*numerical_to_word_labels.xlsx:* For interpretation purposed we provide the correspondence of numerical labels and descriptive labels for the data sets used for the results.


