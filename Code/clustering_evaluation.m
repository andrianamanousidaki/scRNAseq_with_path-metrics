%Clustering evaluation function.
%Given predicted and true cluster labels as input,
%we get ECA,ECP,ARI, number of clusters as output.
%Remark:
%
% **The input must be a categorical vector. Use categorical() if
%   needed.
% **Also need addpath(genpath('cmccomb-rand_index-10a43e9'))
%

function f= clustering_evaluation(predicted_clusters, true_clusters)
   
   
    function f= entropy(x)
        
    freqs = histcounts(x)/sum(histcounts(x));
    freqs = freqs(freqs>0);
    f=-dot(freqs,log(freqs));
    end

    ARI = rand_index(predicted_clusters, true_clusters, 'adjusted');
    ECP=mean( arrayfun(@(x) entropy(predicted_clusters(true_clusters == x)), unique(true_clusters))); 
    ECA=mean( arrayfun(@(x) entropy(true_clusters(predicted_clusters == x)), unique(predicted_clusters)));
    
    number_pred_clusters= size(unique(predicted_clusters));
    number_pred_clusters=number_pred_clusters(:,1);
    
    f=table(ARI,ECP,ECA,number_pred_clusters);

end    