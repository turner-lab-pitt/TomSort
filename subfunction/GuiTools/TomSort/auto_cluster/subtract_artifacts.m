function [ cont_out, artifacts_out ] = subtract_artifacts( cont, fs, art, n, window_start, window_end )
%SUBTRACT_PERIEVENT_MEAN Summary of this function goes here
% Inputs
%	cont - data vector (continuous samples)
%	fs - sampling frequency (or timestamp vector)
%	art - artifact timestamps
%   n - number of artifacts in the moving average (must be odd)
%   window_start - perievent start time relative to event (seconds)
%   window_end - perievent end time relative to event (seconds)
%
%
% Output
%	cont_out - artifact-subtracted signal; same size as input cont
%
% TMP 2018-04
% 

if isscalar(fs)
    artifact_index = round(art(:) * fs); %column vector of artifact indices
else
    artifact_index = find(ismember(fs, art));
    ts = fs;
    fs = 1/mean(diff(fs));
end
window = round(window_start*fs) : round(window_end * fs); %row vector of relative indices
artifact_matind = bsxfun(@plus, artifact_index,window); %2D matrix of linear indices

artifacts = nan(size(artifact_matind));
valid = artifact_matind>0 & artifact_matind<=length(cont);
[r,c]=find(valid);
v=sub2ind(size(valid),r,c);
artifacts(v) = cont(artifact_matind(v));
hanwind = hanning(length(window));
if isempty(artifacts)
    artifacts_out = [];
    cont_out = cont;
elseif size(artifacts,1) <= n
    avg = bsxfun(@times, mean(artifacts), hanwind');
    artifacts_out = repmat(mean(artifacts),size(artifacts,1),1);
    artifact_vector = zeros(size(cont));
    artifact_vector(artifact_matind(:)) = avg(:);
    cont_out = cont - artifact_vector;
else
    padding_length = round( (n-1)/2 );
    padded_artifacts = [flipud(artifacts(1:padding_length,:)); artifacts; artifacts(end:-1:end-padding_length,:)];

    
    indices_to_average = bsxfun(@plus, (1:length(art))', (0:n-1));

    a = reshape(indices_to_average,1,length(art),[]);

    r = bsxfun(@plus, repmat((1:length(art))',1,length(window),n), reshape((0:n-1),1,1,n));
    c = repmat(1:length(window),length(art),1,n);
    linear_ind = sub2ind(size(padded_artifacts),r,c);

    artifact_3D = padded_artifacts(linear_ind);

    moving_avg = bsxfun(@times, mean(artifact_3D,3), hanwind');

    artifact_vector = zeros(size(cont));
    artifact_vector(artifact_matind(:)) = moving_avg(:);
    cont_out = cont - artifact_vector;

    artifacts_out = moving_avg;
end



end

