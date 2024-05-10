function output = recursive_clustering(argstruct)

output=argstruct;
params=argstruct.params;

spc_args=struct;
spc_args.params = params;

% limit to 10 clusters per SPC run
cluster_limit = 10;

num_spikes=argstruct.detection.num_spikes;
F = argstruct.spc.input_features; %features
U = ones(num_spikes,1); %index of unassigned spikes - potential SPC inputs

ns_rel = num_spikes * params.spc.max_spikes_rel;
ns_abs = params.spc.max_spikes_abs;
ns_max = min(ns_rel,ns_abs);

pdist = {};
pthr = {};
gaussians = {};
clusters = zeros(num_spikes,1);
cluster_defs={};        

max_iterations = 10;
num_unclustered = num_spikes;

for iter = 1:max_iterations
    ns=min(ns_max,sum(U));
    ui = find(U);
    if(length(ui)<100), break; end
    
    indices = ui(randperm(length(ui)));
    indices = indices(1:ns);
    spc_args.data = F(indices,:);
    tic;
    fprintf('Running SPC (iteration %d) on %d spikes...',iter,ns);
    spc_results = spc_iteration(spc_args);
%     loaded=load('spcpass1.mat');
%     spc_results=loaded.spc_results;
%     indices=loaded.indices;
    spc = ones(num_spikes,size(spc_results,1)) * -1;
    spc(indices,:) = spc_results(:,3:end)';
    
    figure;hold on;title(['SPC Iteration ' num2str(iter)]);
    plot(F(U==1,1),F(U==1,2),'.k','markersize',0.5);
    drawnow;
    
    fprintf('done in %0.1f minutes.\n',toc/60);
    assign_clusters(spc);
    fprintf('%d spikes assigned to %d clusters so far\n',sum(U==0),length(unique(clusters))-1);
    
    if sum(U==1)==num_unclustered, break; end;
    num_unclustered = sum(U==1);
    
    clust=1:size(clusters,2);
    figure;hold on;title(['Clustering iteration ' num2str(iter)]);
    plot(F(:,1),F(:,2),'.k','markersize',0.5);
    for ci=clust
        if sum(clusters(:,ci)==ci)>0
            plot(F(clusters(:,ci)==ci,1),F(clusters(:,ci)==ci,2),'.','markersize',0.5);
        else
            %plot(F(clusters(:,ci)>0,1),F(clusters(:,ci)>0,2),'.k','markersize',0.5);
        end
    end
    drawnow;
end

    function assign_clusters(r)
        %only previously-unclustered spikes will be included in 'r'
        
        previously_found_clusters=unique(clusters(clusters>0));
        
        ind=false(size(r,1), size(r,2), cluster_limit);
        for ii=1:size(ind,3)
            for jj=1:size(ind,2)
                ind(:,jj,ii)=r(:,jj)==ii-1;
            end
        end
        
        inner_loop_offset=length(pdist);
        
        spc_temps = size(r,2);
        for ii=spc_temps:-1:1
            cluster_counts(ii,:)=reshape(sum(ind(:,ii,:)),1,[]);
            for cl = find(cluster_counts(ii,:)>100)
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
        same_cluster_threshold = 0.5; %proportion of spikes that were previously assigned to a given cluster to call the whole thing part of that cluster
        new_iter_threshold=0.9;
        

        for ii=inner_loop_offset+1:length(pdist)
                    
            found_clusters=unique(clusters(clusters>0));

            dist=pdist{ii};
            if(size(dist)==0), continue; end

            thr = pthr{ii};

            this = dist<thr;
%             this = dist>thr;% & U==1;
    %         cluster_exists = false;
            if sum(this & sum(clusters(:,length(previously_found_clusters)+1:end),2)==0)/sum(this) > new_cluster_threshold
                cluster_index=length(found_clusters)+1;
                clusters(this & sum(clusters(:,length(previously_found_clusters)+1:end),2)==0, cluster_index) = cluster_index;
                cluster_defs{cluster_index} = gaussians{ii};
            end

            
            
        end
        
        new_clusters = setdiff(found_clusters,previously_found_clusters)';
        merged=[];
        for nc=new_clusters
            merge_candidates=sum(clusters(clusters(:,nc)==nc,:)>0)/sum(clusters(:,nc)==nc);
            merge_candidates(nc)=0;%can't merge with itself...
            if sum(merge_candidates) < (1-new_iter_threshold)
                clusters(U==0,nc)=0; %don't re-classify previously-classified spikes
                c=clusters(:,nc)==nc;
                distribution = fitgmdist(F(c,:),1);
                mdist=mahal(distribution,F);
                mthr=prctile(mdist(c),99);
                clusters(U==1 & mdist<mthr,nc)=nc;%expand cluster again
                U(mdist<mthr)=0; %don't put spikes too close to the expanded cluster back into spc
                
                fprintf('Cluster #%d identified as new cluster with %d spikes\n',nc,sum(clusters(:,nc)>0));
            else
                merge_candidates=merge_candidates>same_cluster_threshold;
                
                if sum(merge_candidates)==1 %merge
                    merged(end+1)=nc;
                    merge_target=find(merge_candidates);
                    fprintf('Cluster #%d merged into #%d\n',nc,merge_target);
                    clusters(clusters(:,nc)==nc,merge_target)=merge_target;
                    c=clusters(:,merge_target)==merge_target;
                    distribution = fitgmdist(F(c,:),1);
                    mdist=mahal(distribution,F);
                    mthr=prctile(mdist(c),99);
                    U(mdist<mthr)=0; %don't put spikes too close to the expanded cluster back into spc

                else %not mergable due to overlap with multiple previous; set to nc+0.5 to invalidate
                    clusters(clusters(:,nc)==nc,nc)=nc+0.5;
                end
            end
            %remove merged columns and adjust cluster numbers
            
            
        end
        while ~isempty(merged)
            index=merged(end);
            clusters(:,index:end) = clusters(:,index:end) + (clusters(:,index:end)>0  * -1);
            clusters(:,index)=[];
            merged(end)=[];
        end
        
        U(sum(clusters,2) > 0) = 0;
    end
end