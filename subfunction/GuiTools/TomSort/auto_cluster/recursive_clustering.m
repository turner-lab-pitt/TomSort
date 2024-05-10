function output = recursive_clustering(F, params, recursion)

spc_args=struct;
spc_args.params = params;

% limit to 10 clusters per SPC run
cluster_limit = 10;

% find number of spikes
num_spikes=size(F,1);

if(num_spikes<100)
    fprintf('Warning: less than 100 spikes, cluster_recursively returning immediately.\n');
    output= zeros(1,num_spikes);
    return;
end



ns_max = params.spc.max_spikes_abs;

pdist = {};
pthr = {};
gaussians = {};
clusters = zeros(num_spikes,1);
cluster_defs={};        

ns=min(ns_max,num_spikes);
ui = 1:num_spikes;

indices = ui(randperm(length(ui)));
indices = indices(1:ns);
spc_args.data = F(indices,:);

tic;
fprintf('Running SPC (recursion level %d) on %d spikes...',recursion,ns);
spc_results = spc_iteration(spc_args);
spc = ones(num_spikes,size(spc_results,1)) * -1;
spc(indices,:) = spc_results(:,3:end)';
fprintf('done in %0.1f minutes.\n',toc/60);

[classes, classes_orig] = assign_clusters;
fprintf('%d spikes assigned to %d clusters\n',num_spikes,length(unique(classes(classes>0))));
    
clust=unique(clusters(clusters>0))';
figure;hold on;title(sprintf('Recursion level %d cluster centers',recursion));
plot(F(:,1),F(:,2),'.k','markersize',0.5);
for ci=clust
    if sum(classes_orig==ci)>0
        plot(F(classes_orig==ci,1),F(classes_orig==ci,2),'.','markersize',0.5);
    end
end
figure;hold on;title(sprintf('Recursion level %d expanded clusters',recursion));
plot(F(:,1),F(:,2),'.k','markersize',0.5);
for ci=clust
    if sum(classes==ci)>0
        plot(F(classes==ci,1),F(classes==ci,2),'.','markersize',0.5);
    end
end
drawnow;
if length(clust)<2
    fprintf('Only one cluster found (recursion %d). Returning.\n',recursion);
else
    orig_classes = classes;
    classes=zeros(size(classes));
    for ci=clust
        fprintf('Recursively clustering class %d of %d.\n',ci,length(clust));        
        rc = recursive_clustering(F(orig_classes==ci,:), params, recursion+1);
        classes(orig_classes==ci)=rc + double(rc>0)*max(classes); %classes from each sub-clustering should come back 1-N
    end
    figure;hold on;title(sprintf('Recursion level %d after subclustering',recursion));
    plot(F(:,1),F(:,2),'.k','markersize',0.5);
    for ci=unique(classes(classes>0))'
        if sum(classes==ci)>0
            plot(F(classes==ci,1),F(classes==ci,2),'.','markersize',0.5);
        end
    end
end

output = classes;

    function [final, initial] = assign_clusters
        
        ind=false(size(spc,1), size(spc,2), cluster_limit);
        for ii=1:size(ind,3)
            for jj=1:size(ind,2)
                ind(:,jj,ii)=spc(:,jj)==ii-1;
            end
        end
        figure;
        plot(squeeze(log10(sum(ind))));
        title(sprintf('SPC, recursion %d',recursion));
        
        spc_temps = size(spc,2);
        for ii=spc_temps:-1:1
            cluster_counts(ii,:)=reshape(sum(ind(:,ii,:)),1,[]);
            if ii==spc_temps, continue; end
            for cl = find(cluster_counts(ii,  :)>params.spc.min_spc_cluster &...
                          sum(cluster_counts(ii:end,:)>params.spc.min_spc_cluster)>0);
                c = ind(:,ii,cl);
                distribution = fitgmdist(F(c,:),1);
                gaussians{end+1}=distribution;
                prob=pdf(distribution,F);

                thr = prctile(prob(c,:), 1);
                mthr=prctile(mahal(distribution,F(c,:)),99);
                mdist=mahal(distribution,F);
%                 pdist{end+1} = prob;
%                 pthr{end+1} = thr;
                pdist{end+1} = mdist;
                pthr{end+1}=mthr;
            end
        end

        new_cluster_threshold = 0.9; %proportion of spikes that need to be previously unclustered to call it a new cluster
        cindices=zeros(size(ind,1),1);
        for ii=1:length(pdist)                    
            found_clusters=unique(clusters(clusters>0));

            dist=pdist{ii};
            thr = pthr{ii};
            this = dist<thr; %assign all spikes less than 
            
            if sum(this & clusters==0)/sum(this) > new_cluster_threshold
                cluster_index=length(found_clusters)+1;
                clusters(this & clusters==0) = cluster_index;
                cluster_defs{cluster_index} = gaussians{ii};
            end
        end
        
        if max(clusters)==1
            c = ind(:,floor(spc_temps/3),1);
            distribution = fitgmdist(F(c,:),1);
            
            mthr=prctile(mahal(distribution,F(c,:)),99);
            mdist=mahal(distribution,F);
%                 pdist{end+1} = prob;
%                 pthr{end+1} = thr;
            initial=double(c);
            final=double(mdist<mthr);
        else
            for ii=1:max(clusters)
                mu(ii,:) = cluster_defs{ii}.mu;
                sigma(:,:,ii) = cluster_defs{ii}.Sigma;
            end
            gaussian_mixture=gmdistribution(mu,sigma); %add third param p for mixing proportions?

            posterior_prob = posterior(gaussian_mixture, F);
            [~,clusters_posterior]=max(posterior_prob,[],2);

            dist = pdf(gaussian_mixture,F);
            thr = prctile(pdf(gaussian_mixture,F(clusters>0,:)),1);
            initial = clusters;
            final = clusters_posterior;
        end
    end
end