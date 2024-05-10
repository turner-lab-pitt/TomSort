function output = permutation_clustering(F, params)

spc_args=struct;
spc_args.params = params;
exclude_outliers = isfield(params.clustering,'exclude_outliers') && params.clustering.exclude_outliers;
cocluster_thresholds = params.clustering.cocluster_thresholds;

% limit to 10 distinct clusters per SPC run (should be EXCEEDINGLY unlikely)
cluster_limit = 10;

% find number of spikes
num_spikes=size(F,1);

if(num_spikes<params.spc.min_cluster)
    fprintf('Warning: less than %d spikes, permutation_clustering not run.\n',params.spc.min_cluster);
    output=struct('id',zeros(1,num_spikes),'info','Too few spikes to cluster');
    return;
end

ns_max = params.spc.num_spikes_for_spc; 
ns = min(ns_max,num_spikes);
num_spc_per_permutation=floor(num_spikes/ns);
num_spc_reps = params.spc.num_spc_reps;

% how many permutations to do? at least 2 - one in temporal order, one
% randomized. Enough times that a minimum of num_spc_reps are done.
num_reps = max([2 ceil(num_spc_reps/num_spc_per_permutation)]);
index=zeros(1,num_spikes);
index(1:ceil(num_spikes/num_spc_per_permutation):end)=1;
index=cumsum(index);
ui = 1:num_spikes;

indices=cell(num_reps,num_spc_per_permutation);
for ii=1:num_reps
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%
 
%     if ii<=4
%         [~,rp] = sort(F(:,ii)); %sort spikes by feature values - small clusters show up better
%     elseif ii==5
%         rp = ui; %leave spikes ordered by time of event 
%     else
%         rp=randsample(ui,length(ui));
%     end
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ii==1
        rp = ui;
    else
        rp = randsample(ui,length(ui));
    end
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for jj=1:num_spc_per_permutation
        indices_to_include = rp(index==jj);
        indices{ii,jj} = randsample(indices_to_include,length(indices_to_include));
    end
end

num_spc = numel(indices);

base_filename = spc_args.params.filename;
min_cluster=params.spc.min_cluster;
clusters=nan(num_spc,num_spikes);
cluster_cell = cell(num_spc,1);
args_cell=cell(num_spc,1);

ts = tic;
fprintf('Running SPCx%d in parallel.\n',num_spc);
for ii=1:num_spc
    args_cell{ii}=spc_args;
    args_cell{ii}.params.filename=[base_filename '_' num2str(ii)];
    args_cell{ii}.data=F(indices{ii},:);
end

parfor ii=1:num_spc
    spc_results = spc_iteration(args_cell{ii});
    spc = spc_results(:,3:end)';
    [cluster_cell{ii},~]=assign_clusters(spc,cluster_limit,min_cluster);
end

fprintf('Parallel spc done in %0.1f minutes.\n',toc(ts)/60);

for ii=1:num_spc
    clusters(ii,indices{ii}) = cluster_cell{ii};
end

cluster_defs=cell(0);
cluster_ind = cell(0);
cluster_mat = zeros(num_spikes,0);
for ii=1:num_spc
    classes = unique(clusters(ii,clusters(ii,:)>0));
    for jj=classes
        
        M = F(clusters(ii,:)==jj,:);
        kk=0;
        while rank(M) < size(M,2) %add some noise...
            M(:,end-kk)=M(:,end-kk)+randn(size(M,1),1)*max(std(M(:,end-kk)),std(M(:)));
            kk=kk+1;
        end
        
        cluster_mat(:,end+1)=clusters(ii,:)==jj;
        try
            cluster_defs{end+1} = fitgmdist(M,1);
        catch
            error('SKIP this channel in spreadsheet!!!!!!')
%             disp('Would like to skip this channel, but cant skip due to error with choosing units');
%             cluster_defs{end+1} = fitgmdist(M,1,'Regularize',1e-6);%added large regularization 5/5/2019
        end
        cluster_ind{end+1} = [ii jj];
        
    end
end

% a few spc runs with single large clusters can overly dominate the results
% so exclude large cluster outliers.
spikes_per_cluster=sum(cluster_mat);
outlier_threshold = 5*median(spikes_per_cluster);
large_cluster_index = spikes_per_cluster > outlier_threshold;
cluster_defs = cluster_defs(~large_cluster_index);
dist=zeros(length(cluster_defs));
cluster_mat(:,large_cluster_index)=[];

for ii=1:length(cluster_defs)
    for jj=1:length(cluster_defs)
        dist(ii,jj)=mahal(cluster_defs{ii},cluster_defs{jj}.mu);
    end
end

coclusters = eye(size(dist))==1;

for ii=1:length(cluster_defs)
    cc = coclusters(ii,:);
    oldsum = 0;
    while oldsum < sum(cc)
        oldsum = sum(cc);
        for jj=find(cc)
            d1 = dist(jj,:);
            d2 = dist(:,jj)';
            dmin=min(d1,d2);
            dmax=max(d1,d2);
            c = dmin<cocluster_thresholds(1) & dmax<cocluster_thresholds(2);
    
            cc(c)=true;
        end
    end
    coclusters(ii,:)=cc;
end
superclusters=unique(coclusters,'rows');
if(max(sum(superclusters,1))>1)
    error('Something is wrong: An SPC cluster was assigned to multiple superclusters.')
end
cluster_ids = sum(bsxfun(@times,superclusters,(1:size(superclusters,1))'),1);

mu = zeros(length(cluster_defs),size(F,2));
sigma = zeros(size(F,2),size(F,2),length(cluster_defs));
for ii=1:length(cluster_defs)
    mu(ii,:)=cluster_defs{ii}.mu;
    sigma(:,:,ii)=cluster_defs{ii}.Sigma;
end
assignment_prob = zeros(num_spikes,max(cluster_ids));
if ~isempty(clusters(clusters>0))
    big_gmm=gmdistribution(mu,sigma);
    P = posterior(big_gmm,F);

    for ii=1:max(cluster_ids)
        assignment_prob(:,ii) = sum(P(:,cluster_ids==ii),2);
    end

    % %%%%% Exclude outliers if requested %%%%%% %
    if exclude_outliers
        PDF = -log(pdf(big_gmm,F));
        % Use modified z-score based on median for robustness to outliers
        % http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
        modified_z = (PDF-median(PDF))*0.6745 / median(abs(PDF-median(PDF))); % 
        outliers = modified_z > 5.5;
    else
        outliers = false(size(F,1),1);
    end
else
    outliers = true(size(F,1),1);
end
unit_mat = bsxfun(@eq,assignment_prob,max(assignment_prob,[],2));
% [~,order]=sort(sum(unit_mat(~outliers,:)),'descend'); 
%Don't sort by spike count!!!!! Leads to the cluster definition indices being scrambled.
% unit_mat = unit_mat(:,order);
% assignment_prob = assignment_prob(:,order);

unit = sum(bsxfun(@times,unit_mat,1:size(unit_mat,2)),2);
% if ~isequal(unique(unit(~outliers))', 1:length(unique(unit(~outliers))))
%     error('Problem with assigning unit id!');
% end
unit(outliers)=0;

info = struct;
info.spc0 = sum(cluster_mat,2)==0;
info.unitmat = unit_mat;
info.lowprob = max(assignment_prob,[],2)<0.95;
info.outliers = outliers;
info.assignment_prob=assignment_prob;
info.cluster_indices=cluster_mat==1;
info.cluster_ids = cluster_ids;


output=struct('id',unit,'info',info);

end