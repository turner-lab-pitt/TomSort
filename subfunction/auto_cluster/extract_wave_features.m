function output = extract_wave_features(argstruct)
%Calculates the spike features

detection=argstruct.detection.all_spikes;
params=argstruct.params;
output=argstruct;


feat_before_peak = params.sorting.features_before_peak;
dim_red_feat = params.sorting.dimensionality_reduction_features;
points_before_peak=params.sorting.points_before_peak;

feature_index= (1:dim_red_feat) + points_before_peak-feat_before_peak;
%full_spikes=detection.wave_shapes;
WS = detection.wave_shapes(:,feature_index,:);
feature_spikes=reshape(WS,size(WS,1),[],1);
timestamp = detection.wave_index / detection.timeseries_length;% range [0 1] for plotting in feature space

class = detection.class;

%% Find features that may be useful for clustering
% -- wavelet coefficients: Determined individually for each waveform
% meaning that the whole matrix can be done at once
% -- PCA scores: Depend on which spikes are being analyzed. Do the
% different subsets independently, and apply the resulting coefficients to
% the whole matrix.

%Do wavelet decomposition for potential clusterable features

if isfield(argstruct,'features') && isfield(argstruct.features,'wavelets')
    wavelet_coeffs = argstruct.features.wavelets;
else
    tic;
    fprintf('Extracting wavelets... ');
    wavelet_coeffs = do_wavelet_decomposition(feature_spikes);
    output.features.wavelets = wavelet_coeffs;
    output.features.timing.wavedec=toc;
    fprintf(' finished in %0.1f seconds.\n',output.features.timing.wavedec);
end



%Do principal component decomposition for potential clusterable featuers
tic;
fprintf('Beginning PCA...');
%regular amplitude subset PCA
regular_spikes = class==1;
large_spikes = class==2;
% [PCAcoeffs,PCAscores] = pca(feature_spikes(regular_spikes,:));
reg_coeffs = pca(feature_spikes(regular_spikes,:));
regPCAscores = bsxfun(@minus,feature_spikes,mean(feature_spikes(regular_spikes,:))) / reg_coeffs';
lg_coeffs = pca(feature_spikes(large_spikes,:));
lgPCAscores = bsxfun(@minus,feature_spikes,mean(feature_spikes(large_spikes,:))) / lg_coeffs';

output.features.timing.pca=toc;
fprintf(' finished in %0.1f seconds.\n',output.features.timing.pca);

%basic spike features - amplitude
AMP = detection.wave_amplitude;



largespike_features=[wavelet_coeffs lgPCAscores AMP];
regularspike_features =[wavelet_coeffs regPCAscores AMP];

allfeatureclasses=[...
    repmat({'WC'},1,size(wavelet_coeffs,2)),...
    repmat({'PC'},1,size(regPCAscores,2)),...
    'amplitude'
    ];
allfeatureindex=[...
    1:size(wavelet_coeffs,2),...
    1:size(regPCAscores,2)...
    ];
allfeaturenames = allfeatureclasses;
for i=1:length(allfeatureclasses)-1
    allfeaturenames{i}=[allfeatureclasses{i} num2str(allfeatureindex(i))];
end

tic;
fprintf('Beginning feature selection...');
%feature selection
%features for sorting large spikes
lgFeatSel= select_features(largespike_features,allfeaturenames,large_spikes,timestamp);
%features for sorting regular spikes
regFeatSel = select_features(regularspike_features,allfeaturenames,regular_spikes,timestamp);

output.features.regular = regFeatSel;
output.features.large = lgFeatSel;
fprintf(' done in %0.1f seconds.\n',toc);


