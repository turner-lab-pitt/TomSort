function output = select_features(inputfeatures, inputfeaturenames, index, timestamp)

output.num_spikes = sum(index);
output.features=zeros(size(inputfeatures,1),0);
output.sortable_features=output.features;
if output.num_spikes==0
    return
end

%get the subset of features to use for picking features
selectionfeatures = inputfeatures(index,:);
selection_timestamp = timestamp(index);
%normalize all features for output based on the interquartile range of the
%subset being selected for use with
iqr_selection = iqr(selectionfeatures); %added by HW 5/5/2019 for case when dividing by IQR of 0
iqr_selection(iqr_selection==0) = 1;
normalized_inputfeatures = bsxfun(@rdivide, inputfeatures, iqr_selection)*100;
%continue with the normalized selection features;
selectionfeatures = normalized_inputfeatures(index,:);

hscores=zeros(1,size(selectionfeatures,2));
histograms=zeros(size(selectionfeatures,2), 100);
parfor i=1:length(hscores)
    [~,hscores(i),alias_factor(i)]=hscore(selectionfeatures(:,i),'antialias');
    histograms(i,:)=histcounts(selectionfeatures(:,i),100);
end

%parameters for selecting number of features
hthr=0.001;
max_features=10;
min_features=10;


%sort features by Hartigan's Dip Test score: favor multi-modal features
[~,hindex]=sort(hscores,'descend');

ordered_hscores = hscores(hindex);
ordered_alias_factor=alias_factor(hindex);
ordered_histograms = histograms(hindex,:);
ordered_featurenames = inputfeaturenames(hindex);
ordered_features = selectionfeatures(:,hindex);

%Start with predefined features - PC1, PC2, amplitude
I1 = strcmp(ordered_featurenames,'PC1');
I2 = strcmp(ordered_featurenames,'PC2');
I3 = strcmp(ordered_featurenames, 'amplitude');
N =~(I1 | I2 | I3);
ORDER = [find(I1),find(I2),find(I3),find(N)];
ordered_hscores = ordered_hscores(ORDER);
ordered_alias_factor=ordered_alias_factor(ORDER);
ordered_histograms = ordered_histograms(ORDER,:);
ordered_featurenames = ordered_featurenames(ORDER);
ordered_features = ordered_features(:,ORDER);


%remove features with evidence of aliasing
aliased = struct('name',ordered_featurenames(ordered_alias_factor>1.5),'factor',ordered_alias_factor(ordered_alias_factor>1.5));
ordered_features(:,alias_factor>1.5) = [];
ordered_featurenames(alias_factor>1.5) = [];
ordered_hscores(ordered_alias_factor>1.5) = [];
ordered_histograms(ordered_alias_factor>1.5,:) = [];
ordered_alias_factor(ordered_alias_factor>1.5) = [];

%remove linearly dependent columns
if rank(ordered_features) < size(ordered_features,2)
    ms = min(100,size(ordered_features,1));
    [~,independent_features]=rref(ordered_features(1:ms,:));
    ordered_features = ordered_features(:,independent_features);
    ordered_featurenames = ordered_featurenames(independent_features);
    ordered_hscores = ordered_hscores(independent_features);
    ordered_histograms = ordered_histograms(independent_features,:);
    ordered_alias_factor = ordered_alias_factor(independent_features);
end

num_features=sum(ordered_hscores>hthr);
num_features = min(num_features,max_features);
num_features = max(num_features,min_features);

%Add the timestamp of the events as a "feature" (make it the Nth one selected)
ORDER = [1:num_features num_features:length(ordered_hscores)];
ordered_hscores = ordered_hscores(ORDER);
ordered_alias_factor=ordered_alias_factor(ORDER);
ordered_histograms = ordered_histograms(ORDER,:);
ordered_featurenames = ordered_featurenames(ORDER);
ordered_features = ordered_features(:,ORDER);


ordered_hscores(num_features) = 0;
ordered_alias_factor(num_features)=0;
ordered_histograms(num_features,:)=histcounts(timestamp,100);
ordered_featurenames(num_features)={'timestamp'};
ordered_features(:,num_features)=selection_timestamp';

output.selection.features=ordered_features;
output.selection.featurenames=ordered_featurenames;
output.selection.histograms=ordered_histograms;
output.selection.alias_factor=ordered_alias_factor;
output.selection.hscores = ordered_hscores;

included = ismember(inputfeaturenames,ordered_featurenames);
output.all.features = [normalized_inputfeatures(:,included) timestamp'];
output.all.featurenames = {inputfeaturenames{included} 'timestamp'};
[~,full_feature_order]=ismember(ordered_featurenames,output.all.featurenames);
reorderd_full_feature_matrix=output.all.features(:,full_feature_order);
output.features = reorderd_full_feature_matrix(:,1:num_features);
output.featurenames = ordered_featurenames(1:num_features);
output.sortable_features = ordered_features(:,1:num_features);
output.index=index;
end