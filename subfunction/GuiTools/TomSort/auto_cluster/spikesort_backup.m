function output = spikesort(argstruct)

features = argstruct.spc.input_features;
argstruct.params.clustering.exclude_outliers = true;
argstruct.params.clustering.cocluster_thresholds = [25 50];

if size(argstruct.detection.wave_shapes,1) > 100
    argstruct = do_permutation_clustering(argstruct);
else
    info=struct;
    info.cluster_indices={};
    info.cluster_ids=[];
    argstruct.clustering = struct('id',zeros(size(argstruct.detection.wave_shapes,1),1),'info',info);
    
    argstruct.clustering.outliers.pre_merge_MAD=zeros(size(argstruct.detection.wave_shapes,1),1);
    argstruct.clustering.outliers.post_merge_MAD=zeros(size(argstruct.detection.wave_shapes,1),1);
    output = argstruct;
    return;
end
ids = argstruct.clustering.id;
argstruct.clustering.history.id=ids;
argstruct.clustering.history.event={'Initial clustering'};
cluster_defs = argstruct.clustering.info.cluster_indices;
cluster_def_ids=argstruct.clustering.info.cluster_ids;

if sum(ids==0)>100 && true
    % try reclustering the unsorted spikes
    % small clusters not picked up the first time because of splitting the
    % population into subsets should get identified now.
    % still allow outlier exclusion
    recluster_index=ids==0;
    fprintf('Reclustering %d unsorted spikes\n',sum(ids==0));
    new_cluster_struct = recluster(argstruct,ids==0);
    new_ids = new_cluster_struct.clustering.id;
    new_ids(new_ids>0) = new_ids(new_ids>0) + max(ids); %if there were 2 clusters before (1, 2) make these (3, 4, ...)
    ids(ids==0) = new_ids;
    nci = new_cluster_struct.clustering.info.cluster_indices;
    cluster_defs(recluster_index,end+1:end+size(nci,2))=nci;
    cluster_def_ids(end+1:end+size(nci,2))=new_cluster_struct.clustering.info.cluster_ids+ max(ids);
    argstruct.clustering.history.id(:,end+1)=ids;
    argstruct.clustering.history.event{end+1}='Unsorted reclustered';
    clear new_cluster_struct;
    fprintf('Finished reclustering unsorted spikes\n');
end


% identify clusters with bi- or multi-modal features using Hartigan's dip test
cluster_ids = unique(ids(ids>0));
h_score = zeros(1,length(cluster_ids));
for ii=1:length(cluster_ids)
    id = cluster_ids(ii);
    F = features(ids==id,:);
    h_score(ii) = hscore(F);
end
to_recluster = find(h_score > 2.5);
fprintf('Clusters with bi-/multi-modal features: %d.\n',length(to_recluster));
for ii=1:length(to_recluster)
    id = to_recluster(ii);
    numspikes = sum(ids==id);
    F=features(ids==id,:);
    cluster_def_ids(cluster_def_ids==id)=0;
    if numspikes > 100
        fprintf('Reclustering #%d (%d spikes)\n',id, numspikes);
        arg_copy = argstruct;
        arg_copy.params.clustering.exclude_outliers=false;%these were already called not outliers. don't exclude now.
        for jj=1:5
            %incrementally make it easier to find multiple clusters
            arg_copy.params.clustering.cocluster_thresholds=0.75 * arg_copy.params.clustering.cocluster_thresholds;
            numspikes = sum(ids==id);
            if numspikes > 100
                fprintf('Reclustering #%d (%d spikes)\n',id, numspikes);
                new_cluster_struct=recluster(arg_copy,ids==id);
                recluster_index=ids==id;
                new_ids=new_cluster_struct.clustering.id;
                cluster_ids=unique(new_ids(new_ids>0));
                new_hscores=zeros(1,length(cluster_ids));
                
                for kk=1:length(cluster_ids)
                    this_id = cluster_ids(kk);
                    subF = F(new_ids==this_id,:);
                    new_hscores(kk) = hscore(subF);
                    if new_hscores(kk) < 2.5
                        temp = new_ids;
                        temp(temp~=this_id) = 0;
                        temp(temp==this_id) = temp(temp==this_id) + max(ids); %if there were 2 clusters before (1, 2) make these (3, 4, ...)

                        nci = new_cluster_struct.clustering.info.cluster_indices;
                        ncid = new_cluster_struct.clustering.info.cluster_ids;
                        nci_sub=nci(:,ncid==this_id);
                        cluster_defs(recluster_index,end+1:end+size(nci_sub,2))=nci_sub;
                        cluster_def_ids(end+1:end+size(nci_sub,2))=ncid(ncid==this_id) + max(ids);


                        ids(recluster_index) = ids(recluster_index)+temp;
%                         argstruct.clustering.history.id(:,end+1)=ids;
%                         argstruct.clustering.history.event{end+1}=sprintf('Reclustered #%d',this_id);
                    end
                end
            end
            if max(new_hscores) < 2.5, break; end % Found unimodal populations. Stop iterating.
        end
%         new_ids(new_ids>0) = new_ids(new_ids>0) + max(ids); %if there were 2 clusters before (1, 2) make these (3, 4, ...)
%         ids(ids==id) = new_ids;
        argstruct.clustering.history.id(:,end+1)=ids;
        argstruct.clustering.history.event{end+1}=sprintf('Reclustered #%d',id);
        
        ids(ids==id) = 0; %anything that still hasn't been assigned to unimodal peaks, leave unclustered
        clear new_cluster_struct;
        fprintf('Finished reclustering #%d.\n',id);
    else
        fprintf('Too few spikes found to recluster #%d.\n',id);
    end
end

% %%%%%%% gaussian mixture model of unimodal clusters
% mu = zeros(length(cluster_def_ids),size(features,2));
% sigma = zeros(size(features,2),size(features,2),length(cluster_def_ids));
% for ii=1:length(cluster_def_ids)
%     M = features(cluster_defs(:,ii),:);
%     M_range = range(M);
%     const_column = M_range==0;
%     if sum(const_column) > 0
%         M_noise=bsxfun(@times,randn(size(M)).*0.1*min(M_range(M_range>0)), const_column);
%         M = M+M_noise;
%         features(cluster_defs(:,ii),:) = features(cluster_defs(:,ii),:)+M_noise;
%     end
%     gm = fitgmdist(M,1);
%     mu(ii,:)=gm.mu;
%     sigma(:,:,ii)=gm.Sigma;
% end
% gmm=gmdistribution(mu,sigma);
% P = posterior(gmm,features);
% P_ids=unique(cluster_def_ids(cluster_def_ids>0));
% assignment_prob = zeros(size(features,1),length(P_ids));
% for ii=1:length(P_ids)
% %     assignment_prob(:,ii) = sum(P(:,cluster_def_ids==P_ids(ii)),2);
% %     assignment_prob(:,ii) = max(P(:,cluster_def_ids==P_ids(ii)),[],2);
% 
% % weight the most likely cluster for each supercluster at twice the others
%     assignment_prob(:,ii) = max(P(:,cluster_def_ids==P_ids(ii)),[],2) + sum(P(:,cluster_def_ids==P_ids(ii)),2);
% end
% unit_mat = bsxfun(@eq,assignment_prob,max(assignment_prob,[],2));
% all_ids=sum(bsxfun(@times,unit_mat,1:size(unit_mat,2)),2);
% ids=all_ids;


% %%%%%%%% K-nearest-neighbor classifier to pull all data points into the
% identified clusters
MDL = fitcknn(features(ids>0,:),ids(ids>0),'numneighbors',10);
knnids = resubPredict(MDL);
newids = predict(MDL, features(ids==0,:));
oldids = ids;
ids(oldids==0) = newids;
ids(oldids~=0) = knnids;
all_ids = ids;

argstruct.clustering.history.id(:,end+1)=ids;
argstruct.clustering.history.event{end+1}=sprintf('KNN algorithm');
        
% %%%%%%%% Exclude outliers for merge detection testing %%%%%%%%%
pre_merge_MAD=zeros(size(all_ids));
% for cl=P_ids
%     cli=find(all_ids==cl);
%     if cl==0, continue; end
%     M = features(cli,:);
%     L = log(sum(bsxfun(@rdivide,bsxfun(@minus,M,mean(M)),std(M)) .^2, 2));
%     L = L-median(L);
%     MAD = median(abs(L));
%     outliers = abs(L)>(4*MAD);
%     ids(cli(outliers)) = 0;
%     pre_merge_MAD(cli) = abs(L)/MAD;
% end

% %%%%%%% Look for clusters that can be merged %%%%% %
break_out_of_merge_loop = false;
while break_out_of_merge_loop==false
    break_out_of_merge_loop = true; %if a merge is found, set this back to false.

    clusters = unique(ids(ids>0));
    [rr, cc] = find(triu(ones(length(clusters)),1));

    for ii=1:size(rr,1)
        id1 = ids==clusters(rr(ii));
        id2 = ids==clusters(cc(ii));
        all_id2=all_ids==clusters(cc(ii));
        
        F1 = features(id1,:);
        F2 = features(id2,:);
        
%         ax = mean(F2)-mean(F1);
%         pr1 = F1 * ax' / norm(ax);
%         pr2 = F2 * ax' / norm(ax);
%         pr = [pr1(:);pr2(:)];
%         pr_min=min(pr);
%         pr_max=max(pr);
%         idx=linspace(pr_min,pr_max,1000);
%         cdf1 = sum(bsxfun(@lt,pr1,idx)) / length(pr1);
%         cdf2 = sum(bsxfun(@lt,pr2,idx)) / length(pr2);
%         cmin=min(abs(cdf2 - (1-cdf1)));
%         mi=median(find(abs(cdf2 - (1-cdf1)) == cmin));
% 
%         opt_thr = interp1(1:length(idx),idx,mi);
%         
%         missed_ratio = sum(pr1>opt_thr)/length(pr1); % well-separated clusters will have very
        % little overlap. If this ratio is very small (e.g. <0.01 or 1% of
        % "mis-classified" spikes between the two clusters, they are
        % well-separated and shouldn't be merged. This could be made
        % smarter - maybe a d' value or something - in the future.
%         error_rate=max([sum(pr1>opt_thr)/sum(pr>opt_thr), sum(pr2<opt_thr)/sum(pr<opt_thr)]);
        
        raw_combined_h = hscore([F1;F2]);
        
        
%         resamp_combined_h = max(hs3);
        %resamp_to_orig_ratio = max(hs3./max([hs1; hs2]));
        if raw_combined_h < 2.5
            
            
            RD1 = zeros(5000,10);
            RD2 = RD1;
            if size(F1,1) >= size(RD1,1)
                RD1 = F1(randsample(size(F1,1), size(RD1,1)),:);
            else
                for jj=1:size(F1,2)
                    [E,X]=ecdf(F1(:,jj));
                    R=rand(size(RD1,1),1);
                    S=sum(bsxfun(@gt,E,R'));
                    RD1(:,jj)=X(S);
                    [E,X]=ecdf(F2(:,jj));
                    R=rand(size(RD1,1),1);
                    S=sum(bsxfun(@gt,E,R'));
                    RD2(:,jj)=X(S);
                end
            end
            if size(F2,1) >= size(RD2,1)
                RD2 = F2(randsample(size(F2,1), size(RD2,1)),:);
            else
                for jj=1:size(F2,2)
                    [E,X]=ecdf(F2(:,jj));
                    R=rand(size(RD2,1),1);
                    S=sum(bsxfun(@gt,E,R'));
                    RD2(:,jj)=X(S);
                end
            end

            RD = [RD1; RD2];
            SVMModel = fitcsvm(RD,[ones(length(RD1),1);ones(length(RD1),1)*-1],'KernelFunction','linear','Standardize',true);

            PROJ = bsxfun(@rdivide,bsxfun(@minus,[F1;F2],SVMModel.Mu),SVMModel.Sigma) * SVMModel.Beta + SVMModel.Bias;
            EXPECTED = ones(size(F1,1)+size(F2,1), 1)*-1;
            EXPECTED(1:size(F1,1)) = 1;
            
            AGREE = sign(PROJ)==sign(EXPECTED);
            TRUE_POS_F1 = sum(AGREE(1:size(F1,1)))/size(F1,1);
            TRUE_POS_F2 = sum(AGREE(end-size(F2,1)+1:end))/size(F2,1);
            FALSE_POS_F1 = sum(~AGREE(end-size(F2,1)+1:end))/size(F1,1);
            FALSE_POS_F2 = sum(~AGREE(1:size(F1,1)))/size(F2,1);
            
            DISTINCT = TRUE_POS_F1>0.99 && TRUE_POS_F2>0.99 && FALSE_POS_F1<0.05 && FALSE_POS_F2<0.05;
            
            if DISTINCT == false
    %         if max([raw_combined_h, resamp_combined_h]) < 3
                %found clusters that should be merged
                ids(id2) = clusters(rr(ii));
                all_ids(all_id2) = clusters(rr(ii));
                break_out_of_merge_loop = false; %continue trying to merge.
                fprintf('Merged #%d into #%d\n',clusters(cc(ii)),clusters(rr(ii)));
                argstruct.clustering.history.id(:,end+1)=ids;
                argstruct.clustering.history.event{end+1}=sprintf('Merged %d into %d\n',clusters(cc(ii)),clusters(rr(ii)));
                break; %break this loop and start over, since ids has changed now.
            end
            
            
            
        end
    end
%     argstruct.clustering.history.id(:,end+1)=ids;
%     argstruct.clustering.history.event{end+1}=sprintf('Merged units');
end

% %%%%%%% Sort units by spike count, and renumber them 1-N %%%%% %

orig_ids= all_ids;
ids = zeros(size(ids));
unique_ids= unique(orig_ids(orig_ids>0));
spike_counts = zeros(1,length(unique_ids));
for ii=1:length(spike_counts)
    spike_counts(ii) = sum(orig_ids == unique_ids(ii));
end
[~,order]=sort(spike_counts,'descend');
for ii=1:length(order)
    ids(orig_ids==unique_ids(order(ii)) ) = ii;
end

cluster_nums = unique(ids);
num_clusters=length(cluster_nums);
post_merge_MAD = zeros(size(ids));
for i=1:num_clusters
    % ***************outlier exclusion **************
    cl = cluster_nums(i);
    cli=find(ids==cl);
    if cl==0, continue; end
    M = features(cli,:);
    L = log(sum(bsxfun(@rdivide,bsxfun(@minus,M,mean(M)),std(M)) .^2, 2));
    L = L-median(L);
    MAD = median(abs(L));
    post_merge_MAD(cli) = abs(L)/MAD;
    
end

argstruct.clustering.id=ids;
argstruct.clustering.outliers.pre_merge_MAD=pre_merge_MAD;
argstruct.clustering.outliers.post_merge_MAD=post_merge_MAD;
output = argstruct;
end