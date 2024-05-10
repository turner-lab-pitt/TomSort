function output = find_clusters_wc( filename )
%FIND_CLUSTERS Summary of this function goes here
%   Detailed explanation goes here
if length(filename)>4 && strcmp(filename(end-3:end),'.mat')~=0
    filename = filename(1:end-4);
end
fndg = ['data_' filename '.dg_01.lab'];
fprintf('Loading SPC clustering results from %s... ',fndg);
dg = load(fndg);
fprintf('done\n');

fntimes = ['times_' filename '.mat'];
fprintf('Loading wave features from %s... ',fntimes);
spctimes = load(fntimes, 'ipermut', 'inspk', 'index', 'par', 'spikes', 'cluster_class');
fprintf('done\n');

output=struct;
%output.loaded=struct('dg',dg,'times',spctimes);

%use the ipermut index to remap the columns in dg back to the original
%spike index
%output.permutation_index=spctimes.ipermut;
spc_cluster(:,spctimes.ipermut)=dg(:,3:end);

% min_spc_cluster = 60;
% for ii=1:15,clustersize(ii,:)=sum(spc_cluster==ii,2);end
% num_clusters=max(sum(clustersize>min_spc_cluster));
output.parameters = spctimes.par;
output.wave_features=spctimes.inspk;
output.wave_shapes=spctimes.spikes;
output.wave_times=spctimes.cluster_class(:,2);

% limit to 10 clusters
cluster_limit = 10;
ind=false(size(spctimes.inspk,1), size(dg,1), cluster_limit);
for ii=1:cluster_limit
    for jj=1:size(dg,1)
        ind(:,jj,ii)=spc_cluster(jj,:)==ii-1;
    end
end


pdist = {};
pthr = {};
spc_temps = size(dg,1);
for ii=spc_temps:-1:1
    cluster_counts(ii,:)=reshape(sum(ind(:,ii,:)),1,[]);
    for cl = find(cluster_counts(ii,:)>100)
        c = ind(:,ii,cl);
        distribution = fitgmdist(spctimes.inspk(c,:),1);
        gaussians{ii,cl}=distribution;
        prob=pdf(distribution,spctimes.inspk);
        
        thr = prctile(prob(c,:), 1);
                
        pdist{ii, cl} = prob;
        pthr{ii, cl} = thr;
        
    end
end

new_cluster_threshold = 0.9; %proportion of spikes that need to be previously unclustered to call it a new cluster
same_cluster_threshold = 1; %proportion of spikes that were previously assigned to a given cluster to call the whole thing part of that cluster

clusters = zeros(size(spctimes.ipermut'));
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

posterior_prob = posterior(gaussian_mixture, spctimes.inspk);
[~,clusters_posterior]=max(posterior_prob,[],2);

output.cluster_centers = clusters;
output.cluster_posterior_probabilty = posterior_prob;
output.cluster_expand = clusters_posterior;
output.cluster_multiple_possibile = sum(double(posterior_prob>0.99),2)==0;

expanded_clusters=output.cluster_expand.*double(~output.cluster_multiple_possibile); 

mahal_matrix=zeros(size(spctimes.inspk,1),length(cluster_defs));

for cl=unique(clusters(clusters>0))'
    cl_index = expanded_clusters==cl;
    X=spctimes.inspk(cl_index,:);
    GM=fitgmdist(X, 1);
    mahal_matrix(:,cl)=mahal(GM,spctimes.inspk);
end
mahal_index = sub2ind(size(mahal_matrix),1:length(output.cluster_expand),output.cluster_expand');
output.cluster_mahal_distance=mahal_matrix(mahal_index)';


end

