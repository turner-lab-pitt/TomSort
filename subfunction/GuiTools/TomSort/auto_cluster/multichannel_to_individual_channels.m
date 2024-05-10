function [filenames] = multichannel_to_individual_channels(multichannel_file, data_path)
%% constants

%data_path = 'C:\Users\thomasp\Documents\MATLAB\EphysChannelData\';
%% file loading
if ischar(multichannel_file) && exist(multichannel_file,'file')==2 % <-- filename was given
    file = matfile(multichannel_file);
elseif isa(multichannel_file,'matlab.io.MatFile')
    file = multichannel_file;
else
    error('Input argument must be a MatFile object or a path to a .mat file');
end

datafile_orig = file.Properties.Source;
fprintf('Processing file %s\n',datafile_orig);

filename_components = strsplit(datafile_orig,filesep);
filename = filename_components{end};
matches = regexp(filename,'([a-z]_.*)_(\d)_chans_(\d+)_(\d+)','tokens','ignorecase');
if ~isempty(matches)
    name_components = matches{1};
    name_base = [name_components{1} '_session_' name_components{2} '_chan_'];
    session_filename = [name_components{1} '_' name_components{2} '.mat'];
else
    
    %Test for the case where there is no session number in the string
    matches = regexp(filename,'([a-z]_.*)_chans_(\d+)_(\d+)','tokens','ignorecase');
    if isempty(matches)
        error('The filename pattern does not match the expected pattern.');
    end
    m = matches{1};
    
    matches{1}={m{1}, '1', m{2:end}} ;
    name_components = matches{1};
    name_base = [name_components{1} '_session_' name_components{2} '_chan_'];
    session_filename = [name_components{1} '.mat'];
end
session_filecomponents=filename_components;
session_filecomponents{end}=session_filename;
ch = str2double(name_components(3:4));
channels = ch(1):ch(2);

file_details = whos(file);
hp_cont_index = strcmp('hp_cont',{file_details.name});
if ~any(hp_cont_index)
    error('The file does not have the expected field: hp_cont');
end
hp_cont_details = file_details(hp_cont_index);
hp_cont_size = hp_cont_details.size;
hp_cont_bytes= hp_cont_details.bytes;
bytes = hp_cont_bytes;
prefixes = {'','kilo','mega','giga'};
for ii=1:length(prefixes)
    prefix = prefixes{ii};
    if bytes/1024 < 1
        break;
    end
    bytes = bytes/1024;
end

%check to make sure we know what channels these are
[num_channels,dimension] = min(hp_cont_size);
if length(channels) ~= num_channels
    error('hp_cont contains %d channels of data. The filename indicated %d channels (%d:%d)',...
        num_channels,length(channels),channels(1),channels(end));
end

channels_to_process = 1:num_channels;
for ii=1:num_channels
    fname = [data_path name_base num2str(channels(ii)) '.mat'];
    if exist(fname,'file')
        channels_to_process(channels_to_process==ii)=[];
    end    
end

if isempty(channels_to_process)
    fprintf('The individual channel data files already exist in %s. Skipping this file.\n',data_path);
else
    fprintf('Loading data from matfile. Size of hp_cont is %d by %d (%0.2f %sbytes)\n',...
        hp_cont_size(1),hp_cont_size(2),bytes,prefix);
    load_start = tic;
    hp_cont = file.hp_cont;
    samplerate = file.samplerate; %#ok<*NASGU>
    task_codes = struct;
    stim=struct('gpi',[],'bcx',[],'cp',[]);
    
    fields = fieldnames(file);
        
    if ismember('tasks',fields)
        task_codes=file.tasks;
    else
        session_matfile=matfile(strjoin(session_filecomponents,filesep));
        task_codes=session_matfile.tasks;
    end
    if ismember('gpis',fields), stim.gpi = file.gpis; end
    if ismember('bcxs',fields), stim.bcx = file.bcxs; end
    if ismember('cps',fields), stim.cp = file.cps; end
    
    load_length = toc(load_start);
    fprintf('Load completed in %.1f seconds\n',load_length);
    fprintf('Saving data to individual files in %s...\n',data_path);
    filenames = {};
    save_start = tic;
    for ii=channels_to_process
        fname = [data_path filesep name_base num2str(channels(ii))];
        if(dimension == 1)
            signal = hp_cont(ii,:);
        else
            signal = hp_cont(:,ii);
        end
        
        [signal,stim_artifacts] = subtract_gpi_bcx_cp_artifacts(signal,samplerate,stim.gpi,stim.bcx,stim.cp);
        [signal, code_artifacts] = subtract_artifacts(signal,samplerate,task_codes.E_SYSTEMCODES_START,11,-0.04,0.04);
        % Save the files using the -v6 flag, to save lots of time by not
        % compressing the file.
        save(fname,'-v6','signal','samplerate','datafile_orig','stim_artifacts','code_artifacts');
        filenames{end+1} = fname;
    end
    save_length = toc(save_start);
    fprintf('Writing files completed in %.1f seconds\n',save_length);
end



end