function [ gui ] = SpikeSorter( )
%SPIKESORTER Summary of this function goes here
%   Detailed explanation goes here

gui = figure('Name','SpikeSorter','resize','off');

%defaults:
results_path='S:\Monkey_data\SpikeSorting\SortResults\';
data_path = 'S:\Monkey_data\SpikeSorting\EphysChannelData\';
init_path = 'S:\Monkey_data\';

%% set up layout of the widget
layout.base.height=420;
layout.lineheight = 15;
layout.uigettarget=[20 layout.base.height-50 100 15];
layout.uigetdir = [20 layout.uigettarget(2)-20 100 15];
layout.dirinfo = [20 layout.uigetdir(2)-(5+layout.lineheight) 600 layout.lineheight];
layout.fileinfo = [20 20 320 layout.dirinfo(2)-(5+20)];
layout.gobutton = [...
    sum(layout.fileinfo([1 3]))+20,... %x
    sum(layout.fileinfo([2 4])) - 20,...%y
    100, 15]; %w, h

%% create ui elements

tg = uitabgroup(gui);
sortTab = uitab(tg,'title','Files to sort');
multichannelTab = uitab(tg,'title','Convert multi-channel files to single channels');

%% sortable files

uiSortDirname = uicontrol(sortTab,'style','text',...
    'position',layout.dirinfo,...
    'horizontalalignment','left','string','No data selected');

uiSortFileinfo = uicontrol(sortTab,'style','listbox',...
    'position',layout.fileinfo,...
    'fontname','monospaced');
    
uicontrol(sortTab,'style','pushbutton','string','Select data',...
    'callback',@selectDirectory,...
    'position',layout.uigetdir,...
    'userdata',struct('dirname',uiSortDirname,'listbox',uiSortFileinfo,'filetype','sortable'));
uicontrol(sortTab,'style','text',...
    'position',layout.uigetdir + [layout.uigetdir(3)+10 0 250 0],...
    'horizontalalignment','left','string','Must contain variables: signal, samplerate, datafile_orig');    
uicontrol(sortTab,'style','pushbutton','string','Spikesort!',...
    'position',layout.gobutton,...
    'callback',@sortItems);

uicontrol(sortTab,'style','pushbutton','string','Target directory',...
    'callback',@pickSortTarget,...
    'position',layout.uigettarget);
uiresultspath = uicontrol(sortTab,'style','text',...
    'position',layout.uigettarget + [layout.uigettarget(3)+10 0 250 0],...
    'horizontalalignment','left','string',results_path);

%% multichannel files
uiMultichannelDirname = uicontrol(multichannelTab,'style','text',...
    'position',layout.dirinfo,...
    'horizontalalignment','left','string','No data selected');

uiMultichannelFileinfo = uicontrol(multichannelTab,'style','listbox',...
    'position',layout.fileinfo,...
    'fontname','monospaced');

uicontrol(multichannelTab,'style','pushbutton','string','Select data',...
    'callback',@selectDirectory,...
    'position',layout.uigetdir,...
    'userdata',struct('dirname',uiMultichannelDirname,'listbox',uiMultichannelFileinfo,'filetype','multichannel'));

uicontrol(multichannelTab,'style','text',...
    'position',layout.uigetdir + [layout.uigetdir(3)+10 0 250 0],...
    'horizontalalignment','left','string','Must contain variables: hp_cont (matrix), samplerate');    

uicontrol(multichannelTab,'style','pushbutton','string','Split multichannel file(s)',...
    'position',layout.gobutton,...
    'callback',@processMultichannel);

uicontrol(multichannelTab,'style','pushbutton','string','Target directory',...
    'callback',@pickMultichannelTarget,...
    'position',layout.uigettarget);
uidatapath = uicontrol(multichannelTab,'style','text',...
    'position',layout.uigettarget + [layout.uigettarget(3)+10 0 250 0],...
    'horizontalalignment','left','string',data_path);

% Make the listbox less clunky: https://undocumentedmatlab.com/blog/smart-listbox-editbox-scrollbars
try  % graceful-degradation for future compatibility
   % Get the Java scroll-pane container reference
   jScrollPane = findjobj(uiSortFileinfo); 
   % Modify the scroll-pane's scrollbar policies
   % (note the equivalent alternative methods used below)
   set(jScrollPane,'VerticalScrollBarPolicy',javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);   %VERTICAL_SCROLLBAR_AS_NEEDED=20
   %jScrollPane.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);  %HORIZONTAL_SCROLLBAR_AS_NEEDED=30
catch
   % Never mind...
end

%%
filepath = 'S:\Monkey_data\';

%% Callbacks
    function pickSortTarget(~,~,~)
        p = uigetdir(results_path,'Select folder for saving sort results');
        if ~isempty(p) && ischar(p)
%             disp(['Why are we here? ' num2str(~isempty(p) && ischar(p))])
            results_path=p;
            set(uiresultspath,'string',p);
        end
    end
    function pickMultichannelTarget(~,~,~)
        p = uigetdir(data_path, 'Select folder for saving channel data');
        if ~isempty(p) && ischar(p)
            data_path=p;
            set(uidatapath,'string',p);
        end
    end
    function selectDirectory(btn,~,~)
        userdata = get(btn,'userdata');
        dirname = userdata.dirname;
        listbox = userdata.listbox;
        path = uigetdir2(init_path,'Select files or folders');
        %disp(path);
        if isempty(path) || isnumeric(path)
            set(dirname,'string','No directory/file chosen')
            set(listbox,'string',{});
            return;
        elseif length(path)==1
            set(dirname,'string',sprintf('Selected: %s',path{1}))
        else
            set(dirname,'string',sprintf('%d files/directories selected',length(path)));
        end
        
        %fileinfo = cell.empty();
        data = listbox.get('userdata');
        data.filedata = {}; %clear the old file info
        set(listbox,'userdata',data);
        set(listbox,'string',{}); %clear the old list
        for ii=1:length(path)
            p = path{ii};
            if ~isdir(p)
                [filetype,f,valid] = validateFile(p);
                if valid && strcmp(userdata.filetype,filetype)
                    addFileToList(f,listbox);
                end
            else %isdir==true
                dirContents = dir([p '\*.mat']); % list only matfiles
                for jj=1:length(dirContents)
                    fp = [p filesep dirContents(jj).name];
                    [filetype,f,valid] = validateFile(fp);                    
                    if valid && strcmp(userdata.filetype,filetype)
                        addFileToList(f,listbox);
                    end
                end
            end
        end
        
        %set(uiFileinfo,'string',fileinfo);
    end
    function sortItems(~,~,~)
        data = get(uiSortFileinfo,'userdata');
        if isfield(data,'filedata')
            files = data.filedata;
            for ii=1:length(files)
                file = files{ii};
                %disp(file);
                handleFileInput('sortable',file);
                string = get(uiSortFileinfo,'string');
                string(1)=[];
                set(uiSortFileinfo,'string',string);
                drawnow;
                fd = files;
                fd(1)=[];
                data.filedata=fd;
                set(uiSortFileinfo,'userdata',data);
            end
        end
    end
    function processMultichannel(~,~,~)
        data = get(uiMultichannelFileinfo,'userdata');
        if isfield(data,'filedata')
            files = data.filedata;
            for ii=1:length(files)
                file = files{ii};
                %disp(file);
                handleFileInput('multichannel',file);
                string = get(uiMultichannelFileinfo,'string');
                string(1)=[];
                set(uiMultichannelFileinfo,'string',string);
                drawnow;
                fd = files;
                fd(1)=[];
                data.filedata=fd;
                set(uiMultichannelFileinfo,'userdata',data);
            end
        end
    end


    function handleFileInput(filetype,matobj)
        
        switch filetype
            case 'multichannel'
                multichannel_to_individual_channels(matobj, data_path);
            case 'sortable'
                current_pwd=pwd;
                cd('S:\data_Pearce\MATLAB\GuiTools\SpikeSorter\working_directory');
                do_sorting_and_save(matobj, results_path);
                cd(current_pwd);
            otherwise
                warning('Unrecognized file type');
        end
    end
end
function do_sorting_and_save(matobj, results_path)
output = sort_channel(matobj);
inputfile = strsplit(matobj.Properties.Source,filesep);
inputfile = inputfile{end};
resultsfile = matfile( [results_path filesep inputfile],'Writable',true );
resultsfile.SorterOutput=output;
resultsfile.numUnits = output.sort_results.num_units;
keyfile = matfile( [results_path filesep 'key.mat'],'Writable',true);
try
    results_cell = keyfile.results_cell;
catch E
    results_cell = {};
end
results_cell{end+1} = struct('filename',inputfile,'num_units',output.sort_results.num_units);
keyfile.results_cell = results_cell;

% dirstatus = mkdir([results_path filesep 'sorted']);
% if dirstatus==1
%     status = 0;
%     tries_left = 20;
%     while status == 0 && tries_left>0
%         status = movefile(matobj.Properties.Source,[ results_path filesep 'sorted']);
%         tries_left = tries_left-1;
%         if status==0
%             warning('There was a problem moving the file. Trying again in 500 ms');
%             pause(0.5);
%         end
%     end
% end

fprintf('Finished sorting channel %s\n*******\n\n',matobj.Properties.Source);

end
function var_info = checkForVariable(matobj,name)
    details = whos(matobj);
    var_info = 0;
    for ii=1:length(details)
        var = details(ii);
        if strcmp(name, var.name)
            var_info = var;
            return;
        end
    end
end

function [filetype,matobj,valid] = validateFile(filename)
    filetype = 'Not recognized';
    valid=false;
    matobj = matfile(filename);

    % If the file contains 'hp_cont' and 'samplerate', it might be a
    % multichannel
    hp_cont = checkForVariable(matobj,'hp_cont');
    samplerate = checkForVariable(matobj,'samplerate');
    if ~isstruct(samplerate)
        ts = checkForVariable(matobj,'ts');%Intan data files may have ts instead of samplerate
        if isstruct(ts)
            samplerate = 1 / mean(diff(matobj.ts));
            matobj.Properties.Writable=true;
            matobj.samplerate = samplerate;
        end
    end
    % If hp_cont is a vector already, this file is 'sortable'
    % rather than a multichannel file to split into single-channel files
    if isstruct(hp_cont) && ~isvectorFromSize(hp_cont.size) && (isstruct(samplerate) || isstruct(ts))        
        filetype = 'multichannel';
        valid=true;
        return;
    end

    % If the file contains ('signal' OR 'hp_cont (vector only)') AND ('samplerate')
    % it is a file that can be used for the sorter.
    % If it came from a larger file, the 'datafile_orig' variable indicates
    % the source
    signal = checkForVariable(matobj,'signal');
    hp_cont = checkForVariable(matobj,'hp_cont');
    sr = checkForVariable(matobj,'samplerate');
%     dfo= checkForVariable(matobj,'datafile_orig');
    if (isstruct(signal) || isstruct(hp_cont)) && isstruct(sr)
        filetype = 'sortable';
        valid=true;
        return;
    end
      

end

%% Helpers
    function l = addFileToList(f,l)
        variables = whos('-file',f.Properties.Source);
        MB = sum([variables.bytes])/(1024*1024);
        
        filename = strsplit(f.Properties.Source,filesep);
        filename = filename{end}(1:end-4);
        filelist = get(l,'string');            
        filelist{end+1} = sprintf('%30s: %0.2f MB',filename,MB);
        data = get(l,'userdata');            
        data.filedata{end+1}=f;                        
        set(l,'string',filelist);
        set(l,'userdata',data);
        drawnow;
        
    end


function [pathname] = uigetdir2(start_path, dialog_title)
% Pick a directory with the Java widgets instead of uigetdir

import javax.swing.JFileChooser;

if nargin == 0 || isempty(start_path) || (isnumeric(start_path)&&startpath == 0) % Allow a null argument.
    start_path = pwd;
end

jchooser = javaObjectEDT('javax.swing.JFileChooser', start_path);

jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
jchooser.setMultiSelectionEnabled(true);
if nargin > 1
    jchooser.setDialogTitle(dialog_title);
end

status = jchooser.showOpenDialog([]);

if status == JFileChooser.APPROVE_OPTION
    jFile = jchooser.getSelectedFiles();
    for ii=1:length(jFile)
        pathname{ii} = char(jFile(ii).getPath());
    end
elseif status == JFileChooser.CANCEL_OPTION
    pathname = {};
else
    error('Error occured while picking file.');
end

end

function isvec = isvectorFromSize(sz)

isvec =  ndims(sz)==2 && ismember(1,sz);

end