function TomSort2Nex(fname,folder,nexfname,nexfolder,fs)

% Input Arguments
%       fname: TomSort result data file
%      folder: folder where the TomSort results are stored
%    nexfname: Nex file name
%   nexfolder: folder to store the nex file
%          fs: sampling rate of the extracellular data

if ~exist('fname','var') || ~exist('folder','var')
    [fname,folder] = uigetfile('*.mat','select file');
end
if ~exist('nexfname','var') || isempty(nexfname)
    nexfname = strrep(fname,'.mat','.nex');
end
if ~exist('sfolder','var') || isempty(nexfolder)
    nexfolder = folder;
end
s = load([folder fname]);
if ~exist('fs','var')
    fs = 30000;
end
unitlabel = 'abcdefghijklmnopqrstuvwxyz';

% Add 1 to n_units so as to include noise
n_units = s.SorterOutput.sort_results.num_units+1;

% start new nex file data
nexFile = nexCreateFileData(fs);

% Step through units 
for u=1:n_units
	
	% get logical indices of each unit type
	if u==1
		% First save unsorted + noise
		wave_inds = s.SorterOutput.sort_results.spike_class < 1;
		unit_name = 'U_wf';
	else
		wave_inds = s.SorterOutput.sort_results.spike_class == u-1;
		unit_name = [unitlabel(u-1),'_wf'];
	end
	
	waveTs = s.SorterOutput.sort_results.timestamps(wave_inds)';
	waves = s.SorterOutput.sort_results.spikes(wave_inds,:)';
	preThresT = 0.001;	% I don't know what this should be
	[nPointsWave,~] = size(waves);
	chan_n = str2num(fname( (regexp(fname,'Ch_')+3):(regexp(fname,'\.','once')-1)));
	unit_n = u-1;
	
	% add waveform with pre-threshold time, num. points in wave, wire and unit numbers
	nexFile = nexAddWaveform(nexFile, fs, waveTs, waves, unit_name, ...
		preThresT, nPointsWave, chan_n, unit_n);
	
end

writeNexFile(nexFile, [nexfolder nexfname]);

