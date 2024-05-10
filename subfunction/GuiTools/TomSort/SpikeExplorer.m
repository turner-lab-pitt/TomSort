function gui = SpikeExplorer()

figtitle = 'SpikeExplorer';

%%Create the figure;
gui = figure('name',figtitle,...
    'resizeFcn',@resizeCallback,...
    'visible','off');

%% Constants
startingFilepath = 'S:\Monkey_data\';

%% "Global" variables
currentlyLoaded = struct('matfile',false,'SorterOutput',false,'baseFilename','');
editorData = struct;
editorData.waveHandle = 0;
editorData.featureHandle=0;

editor_features=0;
editor_waveshapes=0;
editor_isi=0;
editor_text=0;
editor_reclassify = 0;

S = struct; %sort results
clusters=[];
history = {};
wave_shapes = [];
cluster_nums = [];
cluster_types = [];
cluster_classes = [];
num_clusters = [];
peak_amplitude=[];
wave_times =[]; 
wave_features = [];
large_features = [];
wsX = [];
wsY = [];

sort_quality=[];
spike_selector=0;

isibinwidth=0.5;
isiedges=0:isibinwidth:50.5;

redraw_required=0;

%% Create tabs for different user interface functionalties
tgroup=uitabgroup('parent',gui, 'SelectionChangedFcn',@tabChangedCallback);
fileTab = uitab('Title','Select File','parent',tgroup);
summaryTab=uitab('Title','Summary','Parent',tgroup);
summaryPanel=uipanel(summaryTab,'Position',[1 1 98 98]/100, 'Units','Normalized');
waveTab=uitab('Title','Wave features','Parent',tgroup);
largewaveTab=uitab('Title','Large spike features','Parent',tgroup);
thresholdTab = uitab('Title','Threshold','Parent',tgroup);
editorTab = uitab('Title','Sort Editor','Parent',tgroup);


%% set up layout of the widget
layout.base.height=420;
layout.lineheight = 15;
layout.uigetdir = [20 layout.base.height-50 100 15];
layout.dirinfo = [20 layout.uigetdir(2)-(5+layout.lineheight) 400 layout.lineheight];
layout.fileinfo = [20 20 320 layout.dirinfo(2)-(5+20)];
layout.savebutton = [.8 .05 .15 .05];
layout.viewresults = [...
    sum(layout.fileinfo([1 3]))+20,... %x
    sum(layout.fileinfo([2 4])) - 20,...%y
    100, 15]; %w, h
layout.toworkspace = layout.viewresults + [0 -25 0 0];



%% Set up file tab
uicontrol(fileTab,'style','pushbutton','string','Select a data directory',...
    'callback',@selectDirectory,...
    'position',layout.uigetdir,...
    'userdata',struct('fixYpos',layout.base.height-layout.uigetdir(2)));
uiDirname = uicontrol(fileTab,'style','text',...
    'position',layout.dirinfo,...
    'horizontalalignment','left','string','Please select file(s) to process',...
    'userdata',struct('fixYpos',layout.base.height-layout.dirinfo(2)));

% "Listbox" control for 
uiFileinfo = uicontrol(fileTab,'style','listbox',...
    'position',layout.fileinfo,...
    'fontname','monospaced',...
    'userdata',struct('fixTop',...
        layout.base.height-sum(layout.fileinfo([2 4]))));
    
uicontrol(fileTab,'style','pushbutton','string','View Results',...
    'position',layout.viewresults,...
    'callback',@itemChosen,...
    'userdata',struct('fixYpos',layout.base.height-layout.viewresults(2)));

uicontrol(fileTab,'style','pushbutton','string','Export to workspace',...
    'position',layout.toworkspace,...
    'callback',@exportToWorkspace,...
    'userdata',struct('fixYpos',layout.base.height-layout.toworkspace(2)));

%% Add save button to Summary tab
uicontrol(summaryTab,'style','pushbutton','string','Save Changes',...
    'units','normalized','position',layout.savebutton,...
    'callback',@saveChanges)


%% Create Editor tab UI
setupEditorTab();


%% sort editor user interface setup
    function setupEditorTab()
        editor_features=axes(editorTab,'pos',[.5 .5 .5 .5],'units','normalized','xtick',[],'ytick',[]);
        editor_waveshapes=axes(editorTab,'pos',[0 0 .5 .5],'units','normalized','xtick',[],'ytick',[]);
        editor_isi=axes(editorTab,'pos',[.5 0 .5 .5],'units','normalized','xtick',[],'ytick',[]);
        editor_text=uicontrol(editorTab,'style','text','units','pixels','pos',[5 layout.base.height-45 200 15],...
            'userdata',struct('fixYpos',45),...
            'string','Use the Summary tab to select a unit','horizontalalignment','left');
        
        reclassify_panel = uipanel(editorTab,'units','pixels','pos',[5 layout.base.height-95 250 50],...
            'userdata',struct('fixYpos',95));
        uicontrol(reclassify_panel,'style','text','units','pixels','pos',[5 28 100 15],...
            'string','Reclassify entire cluster as:','horizontalalignment','left');
        editor_reclassify = uicontrol(reclassify_panel,'style','popup','string',{'-----'},...
            'pos',[5 10 100 15]);
        uicontrol(reclassify_panel,'style','pushbutton','string','Reclassify',...
            'units','pixels','pos',[115 2 60 25],...
            'callback',@reclassify_cluster);
        
        spikeselect_panel = uipanel(editorTab,'units','pixels','pos',[5 layout.base.height-165 250 60],...
            'userdata',struct('fixYpos',165));
        uicontrol(spikeselect_panel,'style','text','units','pixels','pos',[5 38 100 15],...
            'string','Select spikes:','horizontalalignment','left');
        selector_action = uicontrol(spikeselect_panel,'style','popup','string',{'Add spikes', 'Remove spikes'},...
            'pos',[5 10 100 15]);
        uicontrol(spikeselect_panel,'style','pushbutton','string','Apply to selection',...
            'units','pixels','pos',[115 2 40 25],...
            'callback',{@test_selection,selector_action});
%         uicontrol(spikeselect_panel,'style','pushbutton','string','Confirm',...
%             'units','pixels','pos',[155 2 50 25],...
%             'callback',{@confirm_selection,selector_action});

        uicontrol(spikeselect_panel,'style','pushbutton','string','Point cloud',...
            'units','pixels','pos',[155 32 70 25],...
            'callback',@startImpoly);
        
        uicontrol(spikeselect_panel,'style','pushbutton','string','Waveforms',...
            'units','pixels','pos',[80 32 70 25],...
            'callback',@startImline);
        
        editorData.plotted=false;
    end

    function test_selection(src,evt,action_dropdown)
        if ~isvalid(spike_selector)
            error('Invalid action: you must use a Selection Tool first')
        end
        plotted = editorData.plotted;
        editing_id = plotted(1);
        editing_features=plotted(2:3);
        feat = wave_features(:,editing_features);
        
        if isa(spike_selector,'imline')%selected spikes from waveform plot
            line = spike_selector.getPosition();            
            wx = reshape(wsX',[],1);
            wy = reshape(wsY',[],1);
            [~, ~, ii] = polyxpoly(line(:,1),line(:,2),wx,wy);
            
            [~,spike_row] = ind2sub(size(wsX'),ii(:,2));
            selection = false(size(clusters));
            selection(unique(spike_row)) = true;
            
        else %selected spikes from feature scatter plot
            poly = spike_selector.getPosition();
            
            selection = inpolygon(feat(:,1),feat(:,2),poly(:,1),poly(:,2));
        end
        
        delete(spike_selector)
        
        wx = reshape(wsX(selection,:)',[],1);
        wy = reshape(wsY(selection,:)',[],1);
        plot(editor_waveshapes,wx,wy,'-r');
        plot(editor_features,feat(selection,1),feat(selection,2),...
                    '.','color','r','markersize',.5);
                
        index = get(action_dropdown,'value');
        opts = get(action_dropdown,'string');
        action = opts{index};
        %disp(['Testing ' action]);
        switch action
            case 'Remove spikes'
                remove_spikes(selection,editing_id);
            case 'Add spikes'
                add_spikes(selection,editing_id);
        end
        confirmEditDialog();
        setupSortEditor();
        if strcmp(action,'Add spikes') && editing_id <= 0 && redraw_required           
            tgroup.SelectedTab = summaryTab;
            drawGraphs();
        end
    end
    function add_spikes(selection,editing_id)
        history{end+1} = clusters;
        if editing_id <= 0
            % create a new cluster
            new_cluster_id =  max(1, max(cluster_nums)+1);
            clusters(clusters<=0 & selection) = new_cluster_id;
            cluster_nums(end+1)=new_cluster_id;
            cluster_types(end+1)=1;
            num_clusters = length(cluster_nums);
        else
            % add spikes to the existing cluster
            clusters(clusters<=0 & selection) = editing_id;        
        end
        
        setupSortEditor(true);
    end
    function remove_spikes(selection,editing_id)
        history{end+1} = clusters;
        clusters(clusters==editing_id & selection) = -1;        
        setupSortEditor(true);
    end
    function confirmEditDialog()
        quest='Proceed with these changes?';
        title='Confirm changes';
        btn1='Yes';
        btn2='No (undo)';
        defbtn='Yes';
        answer = questdlg(quest,title,btn1,btn2,defbtn);
        switch answer
            case 'Yes'
                redraw_required=1;
                return;
            case 'No (undo)'
                undo();
                return;
            case ''
                undo();
                return;
        end
    end
    function undo()

        clusters=history{end};
        history(end)=[];
        setupSortEditor(true);

    end
%     function confirm_selection(src,evt,action_dropdown)
%         index = get(action_dropdown,'value');
%         opts = get(action_dropdown,'string');
%         action = opts{index};
%         disp(['Confirming ' action]);
%     end
%%

% Make the edit less clunky: https://undocumentedmatlab.com/blog/smart-listbox-editbox-scrollbars
try  % graceful-degradation for future compatibility
   % Get the Java scroll-pane container reference
   jScrollPane = findjobj(uiFileinfo);
 
   % Modify the scroll-pane's scrollbar policies
   % (note the equivalent alternative methods used below)
   set(jScrollPane,'VerticalScrollBarPolicy',javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);   %VERTICAL_SCROLLBAR_AS_NEEDED=20
   %jScrollPane.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);  %HORIZONTAL_SCROLLBAR_AS_NEEDED=30
catch
   % Never mind...
end

% now that all elements are created, turn the gui on
set(gui,'visible','on');

%% Callbacks
    function saveChanges(src,eventdata)
        disp('save changes called');
        classes={};
        all_classified = true;
        for i=1:length(sort_quality)
            dropdown=sort_quality(i);
            opts = get(dropdown,'string');
            dropdown_index = get(dropdown,'value');
            quality = opts{dropdown_index};
            if strcmpi(quality,'select sort quality')
                all_classified=false;
            end
            userdata=get(dropdown,'userdata');
            cluster_id=userdata.cluster;
            class=struct;
            class.quality=quality;
            class.dropdown_index=dropdown_index;
            class.cluster_id=cluster_id;
            classes{i}=class;
        end
        answer = inputdlg('Enter your initials','Sign off on the sort results',[1 60]);
%         uiwait(answer);
        if isempty(answer)
            %Cancel clicked; silently do nothing
        elseif isempty(answer{1})
            errordlg('Results not saved: You must enter your initials','Warning!!!')
        else
            filename = strrep([currentlyLoaded.baseFilename '_' answer{1} '.mat'],'\','/');
            msg=sprintf('Saving %s\n',filename);
            fprintf(msg);
            h = msgbox(msg);
            S.saved_results.clusters = clusters;        
            S.saved_results.cluster_nums = cluster_nums;
            S.saved_results.cluster_types = cluster_types;
            S.saved_results.classes = classes;
            S.saved_results.base_filename = currentlyLoaded.baseFilename;
            S.saved_results.finalized = all_classified;
            S.saved_results.signed_by = struct('initials',answer{1},'timestamp',now);
            mf=matfile(filename,'writable',true);
            mf.SorterOutput=S;
            mf.numUnits = sum(cluster_nums>0);
            delete(h);
            if all_classified
                for i=1:length(classes)
                    class=classes{i};
                    class.sorter_file=filename;
                    if any(strcmpi(class.quality,{'A','B','C','Fiber'}))
                        save_unit_matfile(class);
                    end
                end
            end
            tgroup.SelectedTab = fileTab;
        end
%         S = struct; %sort results
%         clusters=[];
%         history = {};
%         wave_shapes = [];
%         cluster_nums = [];
%         cluster_types = [];
%         num_clusters = [];
%         peak_amplitude=[];
%         wave_times =[]; 
%         wave_features = [];
%         large_features = [];
%         wsX = [];
%         wsY = [];
    end

    function save_unit_matfile(class)
        cluster_id = class.cluster_id;
        quality = class.quality;
        index = clusters==cluster_id;      
        timestamps = S.sort_results.timestamps(index);
        sorter_file=class.sorter_file;
        fileparts = strsplit(sorter_file,'/');
        fname=fileparts{end}(1:end-4);
        fname = [fname '_unit_' num2str(cluster_id) '.mat'];
        unit_path = 'S:/Monkey_data/Thal_Ctx/Data/Units/';
        mf = matfile([unit_path fname],'writable',true);
        mf.sorter_file = sorter_file;
        mf.cluster_id = cluster_id;
        mf.quality = quality;
        mf.timestamps = timestamps;
    end

    function tabChangedCallback(src, eventdata)
        if strcmp(eventdata.OldValue.Title,'Sort Editor') && redraw_required
            drawGraphs();
        end
        if strcmp(eventdata.NewValue.Title, 'Sort Editor')
            setupSortEditor();
        end
    end

    function selectClusterToEdit(src,eventdata)
        if isa(editorData.waveHandle,'matlab.graphics.axis.Axes') && isvalid(editorData.waveHandle)
            set(editorData.waveHandle,'xcolor','k','ycolor','k','color','w');
        end
        editorData.waveHandle = src;
        set(editorData.waveHandle,'xcolor','r','color',[0.8 1 0.8]);
        disp(['Selected unit # ' num2str(src.UserData.cluster)])
    end
    function selectDimensionsToEdit(src,eventdata)
        if isa(editorData.featureHandle,'matlab.graphics.axis.Axes') && isvalid(editorData.featureHandle)
            set(editorData.featureHandle,'color','w');
        end
        editorData.featureHandle = src;
        set(editorData.featureHandle,'color',[0.8 1 0.8]);
        disp(['Selected dimensions ' num2str(src.UserData.dims)])
    end

    function itemChosen(~,~,~)
        index = get(uiFileinfo,'value');
        data = get(uiFileinfo,'userdata');
        if isfield(data,'filedata')
            files = data.filedata;
            file = files{index};
            tgroup.SelectedTab = summaryTab;
            loadSortResults(file);
            
            drawGraphs();
        end
    end



    function loadSortResults(matfileObj)
        % First clear the old tabs
        delete(allchild(summaryPanel));
        delete(allchild(waveTab));
        delete(allchild(largewaveTab));
        delete(allchild(thresholdTab));       
        cla(editor_features);
        cla(editor_waveshapes);
        cla(editor_isi);
%         setupEditorTab();
        drawnow;
        
        S = matfileObj.SorterOutput;
        currentlyLoaded = struct('matfile',matfileObj,'SorterOutput',S);

        wave_shapes = S.sort_results.spikes;
        if isfield(S,'history')
            history = S.history;
        else
            history = {};
        end
        if isfield(S,'saved_results')
            clusters = S.saved_results.clusters;
            cluster_nums = S.saved_results.cluster_nums;
            cluster_types = S.saved_results.cluster_types;
            cluster_classes = S.saved_results.classes;
            currentlyLoaded.baseFilename = S.saved_results.base_filename;
        else
            clusters = S.sort_results.spike_class;
            cluster_nums = S.sort_results.classes;
            cluster_types = S.sort_results.class_featureset;
            cluster_classes = [];
            currentlyLoaded.baseFilename = matfileObj.Properties.Source(1:end-4);
        end
%         clusters = S.sort_results.spike_class;
%         cluster_nums = S.sort_results.classes;
%         cluster_types = S.sort_results.class_featureset;

        num_clusters = length(cluster_nums);
        peak_amplitude=wave_shapes(:,S.params.sorting.points_before_peak);
        wave_times = S.sort_results.timestamps * 1000; %convert from seconds to milliseconds
        wave_features = S.sort_results.spike_features.set1;
        large_features = S.sort_results.spike_features.set2;

        wsX = repmat([1:size(wave_shapes,2) nan],size(wave_shapes,1),1);
        wsY = [wave_shapes nan(size(wave_shapes,1),1)];
    end

    function drawGraphs()
        % First clear the old tabs
        delete(allchild(summaryPanel));
        delete(allchild(waveTab));
        delete(allchild(largewaveTab));
        delete(allchild(thresholdTab));
        sort_quality=[];
        drawnow;
        
        detection_window=(1:S.params.sorting.dimensionality_reduction_features)+...
            S.params.sorting.points_before_peak - S.params.sorting.features_before_peak;

        
        % Plot the spike shapes, ISI histograms, and amplitude histograms
        vmax=max(max(wave_shapes(clusters>0,detection_window)))*1.1;
        vmin=min(min(wave_shapes(clusters>0,detection_window)))*1.1;
        if isempty(vmax) || vmax-vmin==0
            vmax=max(wave_shapes(:))*1.1;
            vmin=min(wave_shapes(:))*1.1;
        end
        %vbinedges=linspace(vmin,vmax,100);
        
        w=100/(sum(cluster_nums>0)+1);

        for i=1:num_clusters
            cl = cluster_nums(i);

            if cl<-1
                continue;
            elseif cl==-1
                cli = find(clusters==0|clusters==-2);
                clr = find(clusters==-1);
                xpos=100-w;
            elseif cl==0
                continue;
            else
                cn = find(unique(cluster_nums(cluster_nums>0))==cl);
                xpos=(cn-1)*w;
                cli=find(clusters==cl);        
            end
            
            %If there are lots of spikes in the unit, randomly sample
            %for plotting so it isn't so slow
            numspikes = length(cli);
            fullcli = cli;
            if numspikes > 25000
                rs = randsample(numspikes,25000);
                cli=cli(rs);
            end
            
            aWaves=axes('parent',summaryPanel,'position',[xpos 75 w 20]/100,'units','normalized');
            box on;
            aWaves.ButtonDownFcn = @selectClusterToEdit;
            
            hold(aWaves,'on');    
            aAmp=axes('parent',summaryPanel,'position',[xpos 35 w 20]/100,'units','normalized');
            hold(aAmp,'on');
            aISI=axes('parent',summaryPanel,'position',[xpos 55 w 20]/100,'units','normalized');
            hold(aISI,'on');
%             aDropdown = uicontrol('parent',summaryPanel,'style','popup','string',{'test','foo','bar'},...
%                 'position',[xpos 25 w 5]/100,'units','normalized');
            hCatPanel=uipanel('parent',summaryPanel,'position',[xpos+2 22 w-4 8]/100,'units','normalized');
            sort_quality(end+1)=uicontrol('parent',hCatPanel,'style','popup',...
                'position',[5 15 120 15],'units','pixels',...
                'string',{'Select Sort Quality','A','B','C','Fiber','Background','Noise'},...
                'userdata',struct('cluster',cl));
            
            if ~isempty(cluster_classes)
                for class_i=1:length(cluster_classes)
                    class=cluster_classes{class_i};
                    if class.cluster_id==cl
                        set(sort_quality(end),'value',class.dropdown_index);
                        break;
                    end
                end
            end
            
            if cl>0, 
                set(aWaves,'ColorOrderIndex',cn);
            end
            lh = plot(aWaves,reshape(wsX(cli,:)',[],1),reshape(wsY(cli,:)',[],1),'linestyle','-');
            if cl==-1
                set(lh,'color','k');
                lh = plot(aWaves,reshape(wsX(clr,:)',[],1),reshape(wsY(clr,:)',[],1),'linestyle','-');
                set(lh,'color',[0.7 0.7 0.7]);
                plot(aWaves, mean(wave_shapes(clr,:)),':k','linewidth',2)
            end

            plot(aWaves, mean(wave_shapes(cli,:)),'-k','linewidth',2)
            if cl>0
                title(aWaves,[num2str(cl) ': ' num2str(numspikes)]);
            elseif cl<=0
                titlestring = sprintf('Noise floor (%d events)',sum(clusters==-1));
                unclust = sum(clusters==0|clusters==-2);
                if unclust>0
                    titlestring = sprintf('%d not clustered; %s',unclust,titlestring);
                end
                title(aWaves,titlestring);
            end
            %if cl==-1, cli = union(cli,clr); end
            if cl==-1, fullcli = clr; end

            aWaves.UserData = struct('cluster',cl,'color',get(lh,'color'));

            [pa,pa_edges] = histcounts(peak_amplitude(fullcli),100);
            bar(aAmp,pa_edges(1:end-1)+0.5*diff(pa_edges(1:2)),pa,1,'EdgeColor','none','FaceColor','k'); 
            isi=diff(wave_times(fullcli));
            isihist=histcounts(isi,isiedges);
            pct_violations=sum(isi<1.5)/length(isi)*100;
            pct_color = 'k';
            if pct_violations > 0.5, pct_color='r';end
            rectangle(aISI,'position',[0 0 1.5 max(isihist)*1.1],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
            bar(aISI, isiedges(1:end-1)+0.5*isibinwidth, isihist, 1, 'EdgeColor','none','FaceColor','k');
            bar(aISI, isiedges(isiedges<1.5)+0.5*isibinwidth, isihist(isiedges<1.5),1,'EdgeColor','none','facecolor','r');
            text(aISI,50,max(max(isihist),1),sprintf('%0.2f%% <1.5ms',pct_violations),'horizontalalignment','right','color',pct_color);
            set(aWaves,'xlim',[1 size(wave_shapes,2)],'ylim',[vmin vmax],'xtick',[],'ytick',[]);
            set(aAmp,'xtick',[],'ytick',[],'xlim',[min(peak_amplitude) max(peak_amplitude)]);
            set(aISI,'xlim',[0 50],'ylim',[0 max(max(isihist)*1.1,1)],'xtick',[],'ytick',[]);

        end

        
        inputs = struct;
        inputs.clusters = clusters;
        inputs.cluster_types = cluster_types;
        inputs.wave_features = wave_features;
        inputs.large_features = large_features;
        inputs.waveTab = waveTab;
        inputs.largewaveTab = largewaveTab;
        inputs.axesButtonDownFcn = @selectDimensionsToEdit;
        multidimensional_feature_plot(inputs);
        

        
        % Plot the amplitude threshold information
        amplitude_sort_plot(S.detection.amplitude_info,thresholdTab);
%         setupEditorTab();
        redraw_required = 0;
    end


    function selectDirectory(~,~,~)
        path = uigetdir2(startingFilepath,'Select a folder containing sort results');
        %disp(path);
        if isempty(path)
            set(uiDirname,'string','No directory/file chosen')
%             set(uiFilename,'string',{});
            return;
        elseif length(path)==1
            set(uiDirname,'string',sprintf('Directory: %s',path{1}))
        else
            set(uiDirname,'string',sprintf('%d files/directories selected',length(path)));
        end
        
        %fileinfo = cell.empty();
        data = uiFileinfo.get('userdata');
        data.filedata = {}; %clear the old file info
        set(uiFileinfo,'userdata',data);
        set(uiFileinfo,'string',{});
        set(uiFileinfo,'value',1);
        for ii=1:length(path)
            p = path{ii};
            if ~isdir(p)
                f = matfile(p);
                addFileToList(f,uiFileinfo);
            else %isdir==true
                dirContents = dir([p '\*.mat']); % list only matfiles
                for jj=1:length(dirContents)
                    filepath = [p filesep dirContents(jj).name];
                    f = matfile(filepath);                    
                    addFileToList(f,uiFileinfo);
                end
            end
        end
        
        %set(uiFileinfo,'string',fileinfo);
    end

    function exportToWorkspace(~,~,~)
        selection = get(uiFileinfo,'value');
        userdata = get(uiFileinfo,'userdata');
        mf = userdata.filedata{selection};
        if isa(currentlyLoaded.matfile,'matlab.io.MatFile') && isequal(mf, currentlyLoaded.matfile)
            assignin('base','SorterOutput',currentlyLoaded.SorterOutput);
        else
            assignin('base','SorterOutput',mf.SorterOutput);
        end
    end
    function resizeCallback(~,~,~)
        children = findobj(gui);
        figsize = get(gui,'position');
        for ii=1:length(children)
            child = children(ii);
            data = get(child,'userdata');
            if isstruct(data)
                pos = get(child,'pos');
                if isfield(data,'fixYpos')                    
                    pos(2) = figsize(4)-data.fixYpos;
                end
                if isfield(data,'fixTop')
                    pos(4) = figsize(4)-(pos(2) + data.fixTop);
                end
                set(child,'pos',pos);
                
            end
        end
    end

    function setupSortEditor(force_redraw)
        if nargin == 0
            force_redraw=false;
        end
        %disp('setup sort editor called')
        if ~isa(editor_features,'matlab.graphics.axis.Axes') || ~isvalid(editor_features)
            return;
        end
        %plot(editor_features,[0 1],[0 1],'--r')
        if ~isa(editorData.waveHandle,'matlab.graphics.axis.Axes') || ~isvalid(editorData.waveHandle)
            cla(editor_features);
            cla(editor_waveshapes);
            cla(editor_isi);
%             ssetupEditorTab();
        else
            plot_features =...
                isa(editorData.featureHandle,'matlab.graphics.axis.Axes')...
                && isvalid(editorData.featureHandle);
        
            w = get(editorData.waveHandle,'userdata');           
            
            
            editing = clusters == w.cluster;
            noise = clusters <= 0;
            others = ~editing & ~noise;
            
            other_clusters = cluster_nums(cluster_nums ~= w.cluster & cluster_nums>0);
            other_cluster_names=arrayfun(@(s)sprintf('Unit #%d',s), other_clusters,'uni',false);
            other_cluster_names{end+1}='Noise';
            other_clusters(end+1)=-1;
            set(editor_text, 'string', sprintf('Editing unit #%d',w.cluster));
            reclassify_userdata = get(editor_reclassify,'userdata');
            reclassify_userdata.others = other_clusters;
            reclassify_userdata.current = w.cluster;
            set(editor_reclassify,'string',other_cluster_names,'userdata',reclassify_userdata,'value',1);
            
            if plot_features
                f = get(editorData.featureHandle,'userdata');
                dims = f.dims;
            else
                dims = [0 0];
            end
            
            to_plot = [w.cluster, dims];
            
            %If the user hasn't picked a new unit or new features, don't
            %replot - saves time
            if isequal(to_plot,editorData.plotted) && ~force_redraw
                return;
            end
           
            cla(editor_features);
            if plot_features
                feat = wave_features;
                if strcmp(f.type,'large'), feat = large_features; end
                
                hold(editor_features,'on');
                
                plot(editor_features,feat(others,dims(1)),feat(others,dims(2)),...
                    '.','color',[.8 .8 .8],'markersize',.5);
                plot(editor_features,feat(noise,dims(1)),feat(noise,dims(2)),...
                    '.','color',[.6 .6 .6],'markersize',.5);
                plot(editor_features,feat(editing,dims(1)),feat(editing,dims(2)),...
                    '.','color',w.color,'markersize',.5);
            end
            
            % Plot wave shapes
            editorPlotWaveshapes(editing, w)
            
            %Plot ISI histogram
            editorPlotISIHistograms(editing)
            
            %keep track of the unit and features currently plotted
            editorData.plotted = to_plot;
        end
        
    end
    function editorPlotWaveshapes(cluster_idx, waveData)
        % Plot wave shapes
        cla(editor_waveshapes);
        hold(editor_waveshapes,'on');
        others = find(~cluster_idx);
        this = find(cluster_idx);
        
        %Limit to 25k waveshapes for plotting speed
        if length(others)>25000
            rs = randsample(length(others),25000);
            others = others(rs);
        end
        if length(this)>25000
            rs = randsample(length(this),25000);
            this = this(rs);
        end
        wx = reshape(wsX(others,:)',[],1);
        wy = reshape(wsY(others,:)',[],1);
        plot(editor_waveshapes,wx,wy,'-','color',[.7 .7 .7]);
        wx = reshape(wsX(this,:)',[],1);
        wy = reshape(wsY(this,:)',[],1);
        plot(editor_waveshapes,wx,wy,'-','color',waveData.color);
    end
    function editorPlotISIHistograms(cluster_idx)
        %Plot ISI histogram
        cla(editor_isi);
        hold(editor_isi,'on');
        isi=diff(wave_times(cluster_idx));
        isihist=histcounts(isi,isiedges);
        pct_violations=sum(isi<1.5)/length(isi)*100;
        pct_color = 'k';
        if pct_violations > 0.5, pct_color='r';end
        rectangle(editor_isi,'position',[0 0 1.5 max(isihist)*1.1],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
        bar(editor_isi, isiedges(1:end-1)+0.5*isibinwidth, isihist, 1, 'EdgeColor','none','FaceColor','k');
        bar(editor_isi, isiedges(isiedges<1.5)+0.5*isibinwidth, isihist(isiedges<1.5),1,'EdgeColor','none','facecolor','r');
        text(editor_isi,50,max(max(isihist),1),sprintf('%0.2f%% <1.5ms',pct_violations),'horizontalalignment','right','color',pct_color);
        set(editor_isi,'xlim',[0 50],'ylim',[0 max(max(isihist)*1.1,1)],'xtick',[],'ytick',[]);
            
    end

    function startImpoly(src,evt)
        disp('Starting point-cloud selection tool')
        if ~isequal(spike_selector,0) && isvalid(spike_selector)
            delete(spike_selector)
        end
        spike_selector = impoly(editor_features);
    end
    function startImline(src,evt)
        disp('Starting spike selection tool')
        if ~isequal(spike_selector,0) && isvalid(spike_selector)
            delete(spike_selector)
        end
        spike_selector = imline(editor_waveshapes);
    end
        
%% Helpers
    function l = addFileToList(f,l)
        variables = who('-file',f.Properties.Source);
        if ismember('numUnits',variables)
            filename = strsplit(f.Properties.Source,filesep);
            filename = filename{end}(1:end-4);
            filelist = get(l,'string');            
            filelist{end+1} = sprintf('%30s: %d units',filename,f.numUnits);
            data = get(l,'userdata');            
            data.filedata{end+1}=f;                        
            set(l,'string',filelist);
            set(l,'userdata',data);
            drawnow;
        end
    end

    function link(src,~,enlarged_axis,threeD_axis)
        cla(enlarged_axis,'reset');
        data=get(src,'children');
        copyobj(data,enlarged_axis);
        set(enlarged_axis,'xtick',[],'ytick',[]);

        cla(threeD_axis,'reset');
        hold(threeD_axis,'on');
        set(threeD_axis,'xtick',[],'ytick',[]);
        X=[];
        Y=[];
        for ii=1:length(data)
            X=[X data(ii).XData];
            Y=[Y data(ii).YData];
        end
        [~,Xedges,Yedges] = histcounts2(X,Y);
        for ii=1:length(data)
            d = data(ii);
            histogram2(threeD_axis, d.XData, d.YData,Xedges,Yedges,'FaceColor',d.Color,'EdgeColor','none','FaceAlpha',0.5);
        end
    end

    function reclassify_cluster(evt,src)
        dropdown_value = get(editor_reclassify,'value');
        cluster_ids = get(editor_reclassify,'userdata');
        cluster_names=get(editor_reclassify,'string');
        current_id = cluster_ids.current;
        new_id = cluster_ids.others(dropdown_value);
        new_name = cluster_names{dropdown_value};
        fprintf('Reclassifying %d as %s (%d).',current_id,new_name, new_id);
        history{end+1}=clusters;
        clusters(clusters== cluster_ids.current) = new_id;
        idx=cluster_nums==current_id;
        cluster_nums(idx)=[];
        cluster_types(idx)=[];
        num_clusters = length(cluster_nums);
        drawGraphs();
        tgroup.SelectedTab = summaryTab;
    end
            
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