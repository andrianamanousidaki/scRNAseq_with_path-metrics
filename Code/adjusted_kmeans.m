function  klabels=adjusted_kmeans(U,k)

        pm_labels = kmeans(U,k,'replicates',20);% kmeans
        H=hist(pm_labels,1:1:k); % find size of kmeans cluster
        
        if min(H)<sqrt(size(U,1))/2 % if the minimum cluster size is less than 30 cells then
        
            num_tiny_clusters=length(find(H<sqrt(size(U,1))/2));
            [pm_labels, C] = kmeans(U,k+ num_tiny_clusters,'replicates',20); 
            pm_labels_original=pm_labels;
            H=hist(pm_labels,1:1:k+ num_tiny_clusters);
            [H2,idxH]=sort(H);
            tiny_cluster_idx = idxH(1:num_tiny_clusters);
            CD =squareform(pdist(C));% P CENTROID DISTANCE MATRIX 
            options = CD(:,tiny_cluster_idx);
            %if size(tiny_cluster_idx,2)>1
                for i=1:length(tiny_cluster_idx)
                    v=options(:,i);
                    merge_idx = find( options(:,i) == min(v(v>0)));
                    pm_labels(pm_labels==tiny_cluster_idx(i)) = merge_idx;     
                end

            remaining_labels = unique(pm_labels);
            for i=1:length(remaining_labels)
                pm_labels(pm_labels==remaining_labels(i))=i;
            end
            
        end 
        klabels=pm_labels;
end