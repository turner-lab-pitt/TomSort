function [cluster_ids,ind] = assign_clusters(spc, cluster_limit, min_cluster)

    cluster_ids=zeros(1,size(spc,1));

    ind=false(size(spc,1), size(spc,2), cluster_limit);
    for ii=1:size(ind,3)
        for jj=1:size(ind,2)
            ind(:,jj,ii)=spc(:,jj)==ii-1;
        end
    end

%         figure; hold on; plot(squeeze(log10(sum(ind)))); plot([0 18],log10([60 60]),'-k')

    cluster_num = 0;
    cluster_temps=cell(0);
    spc_temps = size(spc,2);
    for ii=spc_temps:-1:1
        cluster_counts(ii,:)=reshape(sum(ind(:,ii,:)),1,[]);
        if ii==spc_temps, continue; end %don't worry about the highest temp
%             for cl = find(sum(cluster_counts(ii:ii+1, :)>min_cluster)==2);
        for cl = find(cluster_counts(ii, :)>min_cluster)
            c = ind(:,ii,cl); %c = indices of the spikes belonging to the spc grouping "cl"
            
            %find previously-id'd clusters that some of these belong to
            prev_id = unique( cluster_ids(cluster_ids > 0 & c') );
            
            if isempty(prev_id) %start a new cluster
                cluster_num = cluster_num+1; %increment cluster count
                cluster_ids(c) = cluster_num; %assign spikes to this cluster
                cluster_temps{cluster_num}{1}=c; %keep track of the tempertures this cluster exists at
            else % at least some spikes have previously been assigned to a cluster
                % figure out what to do with the new spc group
                
                % L = Index of which clusters this group is larger than
                L = false(1,length(prev_id));                
                % E = Index of which clusters this group is an expansion of
                % (>80% of previous spikes in this group now)
                E = false(1,length(prev_id));
                % Populate the L and E indices
                for jj=1:length(prev_id)
                    pid=prev_id(jj);
                    L(jj) = sum(c) > sum(cluster_ids==pid);
                    E(jj) = sum(c' & cluster_ids==pid)/sum(cluster_ids==pid)>0.8;
                end
                
                % If a single cluster satisfies these criteria, expand that
                % cluster.  
                if sum(L & E) == 1
                    pid = prev_id(L & E);
                    cluster_ids(c)=pid; %expand previously found cluster
                    cluster_temps{pid}{end+1}=cluster_ids==pid;%keep track of temp
                end
                % If multiple clusters satisfy them, ignore this
                % cluster because it is just merging 2+ clusters at this
                % lower temperature. 
                if sum(L & E) > 1
                end
                % If no clusters satisfy them, create a
                % new cluster out of the previously-unclustered spikes in
                % the group, without reassigning the others.
                if sum(L & E) == 0
                    cluster_num = cluster_num+1; %increment cluster count
                    cluster_ids(cluster_ids == 0 & c') = cluster_num; %assign spikes to this cluster
                    cluster_temps{cluster_num}{1} = cluster_ids ==cluster_num; %keep track of the tempertures this cluster exists at
                end
                
            end
        end
    end

    cluster_ids=zeros(1,size(spc,1));
    for ii=1:length(cluster_temps)
        len = length(cluster_temps{ii});
        ind = ceil(len/2);
        cluster_ids(cluster_temps{ii}{ind}) = ii;
    end
        
        
end