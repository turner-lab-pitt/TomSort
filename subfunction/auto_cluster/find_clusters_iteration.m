function output = find_clusters_iteration( argstruct )
%FIND_CLUSTERS Summary of this function goes here
%   Detailed explanation goes here
output=argstruct;

if ~isfield(argstruct.spc,'results')
    output.clustering.error='Error in find_clusters: No spc results were found in the input';
    fprintf('Error in find_clusters: No spc results were found in the input');
    return;
end

spc_cluster = argstruct.spc.results(:,3:end);

% limit to 10 clusters
cluster_limit = 10;

ind=false(argstruct.detection.num_spikes, size(spc_cluster,1), cluster_limit);
for ii=1:size(ind,3)
    for jj=1:size(ind,2)
        ind(:,jj,ii)=spc_cluster(jj,:)==ii-1;
    end
end

features=argstruct.spc.input_features;

pdist = {};
pthr = {};
spc_temps = size(spc_cluster,1);
for ii=spc_temps:-1:1
    cluster_counts(ii,:)=reshape(sum(ind(:,ii,:)),1,[]);
    for cl = find(cluster_counts(ii,:)>100)
        c = ind(:,ii,cl);
        distribution = fitgmdist(features(c,:),1);
        gaussians{ii,cl}=distribution;
        prob=pdf(distribution,features);
        
        thr = prctile(prob(c,:), 1);
                
        pdist{ii, cl} = prob;
        pthr{ii, cl} = thr;
        
    end
end

new_cluster_threshold = 0.9; %proportion of spikes that need to be previously unclustered to call it a new cluster
same_cluster_threshold = 1; %proportion of spikes that were previously assigned to a given cluster to call the whole thing part of that cluster

clusters = zeros(argstruct.detection.num_spikes,1);
cluster_defs={};

for ii=size(pdist,1):-1:1
    for jj=1:size(pdist,2)        
        found_clusters=unique(clusters(clusters>0));
        
        dist=pdist{ii,jj};
        if(size(dist)==0), continue; end
        
        thr = pthr{ii,jj};
        
        this = dist>thr;
%         cluster_exists = false;
        if sum(this & clusters==0)/sum(this) > new_cluster_threshold
            cluster_index=length(found_clusters)+1;
            clusters(this & clusters==0) = cluster_index;
            cluster_defs{cluster_index} = gaussians{ii,jj};
        else
            
        end
        
        for kk=1:length(found_clusters)
            if  sum(this & clusters==kk) > same_cluster_threshold
                % merge immediately?

            end
        end
    end
end

for ii=1:max(clusters)
    mu(ii,:) = cluster_defs{ii}.mu;
    sigma(:,:,ii) = cluster_defs{ii}.Sigma;
end
gaussian_mixture=gmdistribution(mu,sigma); %add third param p for mixing proportions?

posterior_prob = posterior(gaussian_mixture, features);
[~,clusters_posterior]=max(posterior_prob,[],2);

output.clustering.initial = clusters;
output.clustering.posterior_probabilty = posterior_prob;
output.clustering.expanded = clusters_posterior;
output.clustering.multiple_possibile = sum(double(posterior_prob>0.99),2)==0;

expanded_clusters=output.clustering.expanded.*double(~output.clustering.multiple_possibile); 

mahal_matrix=zeros(size(features,1),length(cluster_defs));

for cl=unique(clusters(clusters>0))'
    cl_index = expanded_clusters==cl;
    X=features(cl_index,:);
    GM=fitgmdist(X, 1);
    mahal_matrix(:,cl)=mahal(GM,features);
end
mahal_index = sub2ind(size(mahal_matrix),...
    1:length(output.clustering.expanded),output.clustering.expanded');
output.clustering.mahal_distance=mahal_matrix(mahal_index)';


end

