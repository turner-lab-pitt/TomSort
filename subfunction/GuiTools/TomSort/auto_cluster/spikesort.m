function output = spikesort(argstruct)

features = argstruct.spc.input_features;
num_spikes=size(features,1);
argstruct.params.clustering.exclude_outliers = true;
argstruct.params.clustering.cocluster_thresholds = [20 45];
argstruct.params.clustering.hscore_cutoff = 2.0;
hscore_thr = argstruct.params.clustering.hscore_cutoff;

if ~isfield(argstruct,'clustering')
    argstruct.clustering = struct;
end

if num_spikes > argstruct.params.spc.min_cluster
    argstruct = do_permutation_clustering(argstruct);
else
    info=struct;
    info.cluster_indices={};
    info.cluster_ids=[];
    argstruct.clustering.id = zeros(num_spikes,1);
    argstruct.clustering.info = info;
    
    argstruct.clustering.outliers.pre_merge_MAD=zeros(num_spikes,1);
    argstruct.clustering.outliers.post_merge_MAD=zeros(num_spikes,1);
    output = argstruct;
    return;
end
ids = argstruct.clustering.id;
argstruct.clustering.history.id=ids;
argstruct.clustering.history.event={'Initial clustering'};
cluster_defs = argstruct.clustering.info.cluster_indices;
cluster_def_ids=argstruct.clustering.info.cluster_ids;

if sum(ids==0)>argstruct.params.spc.min_cluster
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
to_recluster = find(h_score > hscore_thr);
fprintf('Clusters with bi-/multi-modal features: %d.\n',length(to_recluster));
for ii=1:length(to_recluster)
    id = to_recluster(ii);
    numspikes = sum(ids==id);
    F=features(ids==id,:);
    cluster_def_ids(cluster_def_ids==id)=0;
    if numspikes > argstruct.params.spc.min_cluster
        fprintf('Reclustering #%d (%d spikes)\n',id, numspikes);
        arg_copy = argstruct;
        arg_copy.params.clustering.exclude_outliers=false;%these were already called not outliers. don't exclude now.
        for jj=1:6
            %incrementally make it easier to find multiple clusters
            arg_copy.params.clustering.cocluster_thresholds=0.75 * arg_copy.params.clustering.cocluster_thresholds;
            if jj==6
                arg_copy.params.clustering.cocluster_thresholds=[2 4];
            end
            numspikes = sum(ids==id);
            if numspikes > argstruct.params.spc.min_cluster
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
                    if new_hscores(kk) < hscore_thr
                        temp = new_ids;
                        temp(temp~=this_id) = 0;
                        temp(temp==this_id) = temp(temp==this_id) + max(ids); %if there were 2 clusters before (1, 2) make these (3, 4, ...)

                        nci = new_cluster_struct.clustering.info.cluster_indices;
                        ncid = new_cluster_struct.clustering.info.cluster_ids;
                        nci_sub=nci(:,ncid==this_id);
                        cluster_defs(recluster_index,end+1:end+size(nci_sub,2))=nci_sub;
                        cluster_def_ids(end+1:end+size(nci_sub,2))=ncid(ncid==this_id) + max(ids);


                        ids(recluster_index) = ids(recluster_index)+temp;

                    end
                end
            end
            if max(new_hscores) < hscore_thr, break; end % Found unimodal populations. Stop iterating.
        end
        
        argstruct.clustering.history.id(:,end+1)=ids;
        argstruct.clustering.history.event{end+1}=sprintf('Reclustered #%d',id);
        
        ids(ids==id) = 0; %anything that still hasn't been assigned to unimodal peaks, leave unclustered
        clear new_cluster_struct;
        fprintf('Finished reclustering #%d.\n',id);
    else
        fprintf('Too few spikes found to recluster #%d.\n',id);
    end
end

% %%%%%%%% K-nearest-neighbor classifier to pull all data points into the
% identified clusters; if tiny clusters are left, eliminate them and re-run
% until all spikes are classified into large-enough clusters
do_knn = true;
while do_knn == true
    do_knn = false;
    I = ids>0;
    if sum(I)==0, break; end
    F = features(I,:);
    MDL = fitcknn(F,ids(I),'numneighbors',10);
    knnids = resubPredict(MDL);
    newids = predict(MDL, features(ids==0,:));
    oldids = ids;
    ids(oldids==0) = newids;
    ids(oldids~=0) = knnids;
    ids=ids(:);
    for ID = unique(ids(ids>0))'
        if sum(ids==ID) < argstruct.params.spc.min_cluster
            ids(ids==ID)=0;
            do_knn = true;
        end
    end
end


all_ids = ids;
argstruct.clustering.history.id(:,end+1)=ids;
argstruct.clustering.history.event{end+1}=sprintf('KNN algorithm');

% %%%%%%%% Exclude outliers for merge detection testing %%%%%%%%%
pre_merge_MAD=zeros(size(all_ids));
for cl=unique(all_ids)'
    cli=find(all_ids==cl);
    if cl==0, continue; end
    M = features(cli,:);
    L = log(sum(bsxfun(@rdivide,bsxfun(@minus,M,mean(M)),std(M)) .^2, 2));
    L = L-median(L);
    MAD = median(abs(L));
    outliers = abs(L)>(4*MAD);
    ids(cli(outliers)) = 0;
    pre_merge_MAD(cli) = abs(L)/MAD;
end

% %%%%%%% Look for clusters that can be merged %%%%% %
break_out_of_merge_loop = false;
clusters = unique(ids(ids>0));
comparison_matrix = triu(ones(length(clusters)),1);
comparison_needed = ones(length(clusters));
centroids = zeros(length(clusters),size(features,2)-1);
for ii=1:length(clusters)
    centroids(ii,:) = median(features(ids==clusters(ii),1:end-1));
end
centroid_distances = pdist2(centroids,centroids,'euclidean');
    
while break_out_of_merge_loop==false
    break_out_of_merge_loop = true; %if a merge is found, set this back to false.
    [rr, cc] = find(comparison_matrix & comparison_needed);
    
    [~,centroid_sort]=sort(centroid_distances(sub2ind(size(centroid_distances),rr,cc)),'ascend');
    rr = rr(centroid_sort);
    cc = cc(centroid_sort);
    for ii=1:size(rr,1)
        id1 = ids==clusters(rr(ii));
        id2 = ids==clusters(cc(ii));
        all_id2=all_ids==clusters(cc(ii));
        
        F1 = features(id1,:);
        F2 = features(id2,:);
        
        if isempty(F1) || isempty(F2), continue; end
       
        raw_combined_h = hscore([F1;F2]);

        if raw_combined_h < hscore_thr % consider merging the cases                               
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

            index= id1 | id2;
            PROJ = bsxfun(@rdivide,bsxfun(@minus,features(index,:),SVMModel.Mu),SVMModel.Sigma) * SVMModel.Beta + SVMModel.Bias;


            EXPECTED = zeros(size(features,1), 1);
            EXPECTED(id1) = 1;
            EXPECTED(id2) = -1;
            
            AGREE = sign(PROJ)==sign(EXPECTED(index));
            
            TRUE_POS_F1 = sum(AGREE(1:size(F1,1)))/size(F1,1);
            TRUE_POS_F2 = sum(AGREE(end-size(F2,1)+1:end))/size(F2,1);
            FALSE_POS_F1 = sum(~AGREE(end-size(F2,1)+1:end))/size(F1,1);
            FALSE_POS_F2 = sum(~AGREE(1:size(F1,1)))/size(F2,1);
            proj_hscore = hscore(PROJ);
            
            DISTINCT = proj_hscore>hscore_thr || (TRUE_POS_F1>0.99 && TRUE_POS_F2>0.99 && FALSE_POS_F1<0.05 && FALSE_POS_F2<0.05);
            
            if DISTINCT
            
                comparison_needed(rr(ii),cc(ii))=0;
            else % DISTINCT == false: found clusters that should be merged
                ids(id2) = clusters(rr(ii));
                all_ids(all_id2) = clusters(rr(ii));
               
                break_out_of_merge_loop = false; %continue trying to merge.
                fprintf('Merged #%d into #%d\n',clusters(cc(ii)),clusters(rr(ii)));
%                 argstruct.clustering.history.id(:,end+1)=ids;
%                 argstruct.clustering.history.event{end+1}=sprintf('Merged %d into %d\n',clusters(cc(ii)),clusters(rr(ii)));
             
                comparison_needed(rr(ii),:)=1;
                comparison_needed(:,cc(ii))=0;
                
            end
        else
            comparison_needed(rr(ii),cc(ii))=0;
        end
    end
end
ids = all_ids;
argstruct.clustering.history.id(:,end+1)=ids;
argstruct.clustering.history.event{end+1}=sprintf('Merged subclusters\n');
             

% %%%%%%%%% Use an SVM to reclassify borderline cases
clusters = unique(ids(ids>0));
comparison_matrix = triu(ones(length(clusters)),1);
[rr, cc] = find(comparison_matrix);

for ii=1:size(rr,1)
    id1 = ids==clusters(rr(ii));
    id2 = ids==clusters(cc(ii));

    F1 = features(id1,:);
    F2 = features(id2,:);

    if isempty(F1) || isempty(F2), continue; end

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

    index= id1 | id2;
    PROJ = bsxfun(@rdivide,bsxfun(@minus,features(index,:),SVMModel.Mu),SVMModel.Sigma) * SVMModel.Beta + SVMModel.Bias;

    X = linspace(min(PROJ), max(PROJ),10000);
    KS1 = ksdensity(PROJ(id1(index)),X);
    KS2 = ksdensity(PROJ(id2(index)),X);
    DENS = ksdensity(PROJ,X);
    DENS(X<X(KS2==max(KS2))) = inf;
    DENS(X>X(KS1==max(KS1))) = inf;
    [~,DC] = min(DENS);
    CUTOFF = X(DC);


    split_ids = zeros(size(PROJ));
    split_ids(PROJ<CUTOFF) = clusters(cc(ii));
    split_ids(PROJ>CUTOFF) = clusters(rr(ii));
    ids(index)=split_ids;
    all_ids(index)=split_ids;
    
end

argstruct.clustering.history.id(:,end+1)=ids;
argstruct.clustering.history.event{end+1}=sprintf('SVM reclassifier\n');


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

% cluster_nums = unique(ids);
% num_clusters=length(cluster_nums);
% post_merge_MAD = zeros(size(ids));
% for i=1:num_clusters
%     % ***************outlier exclusion **************
%     cl = cluster_nums(i);
%     cli=find(ids==cl);
%     if cl==0, continue; end
%     M = features(cli,:);
%     L = log(sum(bsxfun(@rdivide,bsxfun(@minus,M,mean(M)),std(M)) .^2, 2));
%     L = L-median(L);
%     MAD = median(abs(L));
%     post_merge_MAD(cli) = abs(L)/MAD;
%     
% end

argstruct.clustering.history.id(:,end+1)=ids;
argstruct.clustering.history.event{end+1}=sprintf('Ordered by spike count\n');

argstruct.clustering.id=ids;
argstruct.clustering.info='Completed';%don't bother keeping all the extra data
output = argstruct;
end