function [ args ] = exclude_outliers( args )
%EXCLUDE_OUTLIERS Summary of this function goes here
%   Detailed explanation goes here

ampclass = args.sorted.amplitude_class;
unit_ids = args.sorted.id;
regfeatures = args.features.regular.features;
lgfeatures = args.features.large.features;

[cluster_ids,ci] = unique(unit_ids);
cluster_class = ampclass(ci);
num_clusters = length(cluster_ids);

for i=1:num_clusters
    % ***************outlier exclusion **************
    cl = cluster_ids(i); %cluster id    
    cli=find(unit_ids==cl); %index of all spikes in this cluster
    class=cluster_class(i); %amp class
    if class==1
        M = regfeatures(cli,:);
    else
        M = lgfeatures(cli,:);
    end
    
    % find N-dimensional distance from cluster mean
    % normalize by dividing by stdev (ie z-score)
    % next, sum the squares of these values, then log-transform to create a
    % distribution of distances from the cluster mean
    L = log(sum(bsxfun(@rdivide,bsxfun(@minus,M,mean(M)),std(M)) .^2, 2));
    L = L-median(L);
    MAD = median(abs(L));
    outliers = abs(L)>(5*MAD);
    
    %Set the cluster ID for these cases to zero (unsorted)
    args.sorted.id(cli(outliers)) = 0;
    
end



end

