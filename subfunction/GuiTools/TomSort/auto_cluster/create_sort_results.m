function A = create_sort_results(args)
A = args; % return the original argument structure plus newly created fields.

amplitudes = A.detection.all_spikes.wave_amplitude;
spike_shapes = A.detection.all_spikes.wave_shapes;
timestamps = A.detection.all_spikes.wave_index / A.params.sr;
ampclass = A.sorted.amplitude_class;
unit_ids = A.sorted.id;
regfeatures = A.features.regular.features;
lgfeatures = A.features.large.features;



% Order the units by mean peak amplitude
unique_unit_ids = unique(unit_ids(unit_ids>0));
if isempty(unique_unit_ids) %none were successfully clustered
    unit_ids(:)=1;
    unique_unit_ids=1;
end
mean_spike_shape = nan(length(unique_unit_ids), size(spike_shapes,2));
for ii=1:length(unique_unit_ids)
    unit_index = unit_ids == unique_unit_ids(ii);
    mean_spike_shape(ii,:) = mean(spike_shapes(unit_index,:));
end

[unit_amplitude,peak_index] = max(abs(mean_spike_shape),[],2);

[~,unit_order]=sort(unit_amplitude,'descend');

[~,ordered_unit_ids] = ismember(unit_ids,unit_order);

% Renumber the highest unit ID to -1, indicating the noise floor "cluster"
ordered_unit_ids(ordered_unit_ids==max(unit_order)) = -1;
classes = unique(ordered_unit_ids);
class_descriptions = {};
class_featureset = zeros(size(classes));
for ii=1:length(classes)
    if classes(ii) == -1
        class_descriptions{ii} = 'Noise floor';
    elseif classes(ii) == 0
        class_descriptions{ii} = 'Unclustered';
    elseif classes(ii) > 0
        class_descriptions{ii} = sprintf('Sorted unit %d',classes(ii));
        class_featureset(ii) = ampclass(find(ordered_unit_ids==classes(ii),1,'first'));
    else
        class_descriptions{ii} = sprintf('Unknown class %d',classes(ii));
    end
end



% Create a "units" field in the output structure, and make it an array of
% structures describing the units
unit_labels = unique(ordered_unit_ids(ordered_unit_ids>0));
%Start with an empty struct array
A.units = struct('label',{},'spikes',{},'mean_spike',{},'timestamps',{});
for ii=1:length(unit_labels)    
    unit = struct;
    idx = ordered_unit_ids==unit_labels(ii);
    unit.label = unit_labels(ii);
    unit.spikes = spike_shapes(idx,:);
    unit.mean_spike = mean(unit.spikes,2);
    unit.timestamps = timestamps(idx);
    A.units(ii) = unit;
end

S = struct;

S.spikes = spike_shapes;
S.timestamps = timestamps;
S.spike_class = ordered_unit_ids;
S.classes = classes;
S.class_featureset = class_featureset;
S.class_descriptions = class_descriptions;
S.num_units = length(A.units);
S.spike_features.set1 = regfeatures;
S.spike_features.set2 = lgfeatures;

A.sort_results = S;



end