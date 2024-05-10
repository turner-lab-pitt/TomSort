function [ output ] = sort_channel( chan_file,flip_flag)
%SORT_CHANNEL Run spike sorting algorithm on a single channel datafile
%   file: string or MatFile pointing to the data file

%% check input argment
if ~exist('flip_flag','var')
    flip_flag = 0;
end

%% file loading
if ischar(chan_file) && exist(chan_file,'file')==2 % <-- filename was given
    file = matfile(chan_file);
elseif isa(chan_file,'matlab.io.MatFile')
    file = chan_file;
else
    error('Input argument must be a MatFile object or a path to a .mat file');
end
%file is now a MatFile object.

%%
params = struct;
params.filtering = struct('type','none');

input = struct;
% Rename data and samplerate variables to match the processing script's
% expected format

fields = fieldnames(file);
if ismember('signal',fields)
    input.data = file.signal;
elseif ismember('hp_cont',fields)
    input.data = file.hp_cont;
else
    error('Input file must have `hp_cont` or `signal` variable')
end
if flip_flag
    input.data = input.data * (-1);
end
if ismember('samplerate',fields)
    input.sr = file.samplerate;
else
    error('Input file must have `samplerate` variable')
end
input.filename = file.Properties.Source;

fprintf('Sorting channel %s\n',input.filename);
sorted = process_channel(input,params);
output = sorted.results{1};



end

