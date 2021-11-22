function [Ap] = ConnectGraph(Ap, X, p)

    % Add in smallest edge connecting each pair of CC's:
    CC = conncomp(graph(Ap));
    NumberCC = size(unique(CC),2);
    for i=1:NumberCC
        CCindex{i} = find(CC==i);
    end
    for i=1:NumberCC
        for j=i+1:NumberCC
            SampleDistances = zeros(length(CCindex{i}),length(CCindex{j}));
            for s=1:length(CCindex{i})
                for q=1:length(CCindex{j})
                    SampleDistances(s,q) = norm(X(CCindex{i}(s),:)-X(CCindex{j}(q),:));
                    %SampleDistances(s,q) = Ap(CCindex{i}(s),CCindex{j}(q));
                end
            end
            MinDistance = min(min(SampleDistances));
            [s0, q0] = find(SampleDistances==MinDistance);
            Ap(CCindex{i}(s0),CCindex{j}(q0))=MinDistance^p;
            Ap(CCindex{j}(q0),CCindex{i}(s0))=Ap(CCindex{i}(s0),CCindex{j}(q0));
        end
    end
    
end

