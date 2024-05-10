function output = process_channel(args, params)
% args (first param) can be a struct with fields sr(sampling rate),
% data(timeseries) and filename(string).
    f = struct;
    if isstruct(args)
        f = args;        
    elseif ischar(args)
        filename = args;
        if length(filename)>4 && strcmp(filename(end-3:end),'.mat')==true
            filename = filename(1:end-4);
        end
        if ~exist([filename '.mat'],'file')
            error([filename '.mat not found'])
        end
        fprintf('Loading timeseries (varname: data) and sampling rate (varname: sr) from %s.mat\n',filename);

        filecontents=load([filename '.mat'],'sr','data');
        f.data = filecontents.data;
        f.sr = filecontents.sr;
        fprintf(' done\n');
    end
    if ~exist('params','var')
        params=struct;
    end
    if ~iscell(params)
        params = {params};
    end
    
    if isfield(f,'loadfromfile') && exist(f.loadfromfile, 'file')
        fromfile = load(f.loadfromfile);
        f.results = fromfile.f.results;
    end
    if ~isfield(f,'results')
        f.results={};
    end
    
    
    default_params.sr=f.sr;  
    default_params.filtering = struct('type','elliptical','fmin',500,'fmax',3000);
    default_params.detection=struct('voltagesign','neg','threshold',3,'artifacts',80,'amplitude_sort',1);
    default_params.sorting=struct('points_before_peak',20,'points_after_peak',44,'features_before_peak',8,'dimensionality_reduction_features',32);
    default_params.interpolate='n';
    default_params.spc.min_temp=0.02;
    default_params.spc.max_temp=0.20;
    default_params.spc.temp_step=0.01;
    default_params.spc.cycles=100;
    default_params.spc.k_nearest_neighbors=11;
    default_params.spc.min_cluster=60;
    default_params.spc.num_spikes_for_spc=15000;
    default_params.spc.num_spc_reps = 7;
    current_time=datetime('now','TimeZone','local','Format','d-MMM-y-HH-mm');
    current_timevec = datevec(current_time);
    
    default_params.filename=sprintf('clustering_%s-%0.3f',current_time,current_timevec(end));
    for ii=1:length(params)
        args=struct;
        args.data=double(f.data);
        params{ii} = merge_structures(default_params, params{ii});
        args.params=params{ii};
        
        fprintf('Filtering raw data using filter type: %s...\n',args.params.filtering.type);
        args = filter_signal(args);
        fprintf('finished in %0.1f seconds.\n', args.timing.filtering);
        args = detect_spikes(args);

        args = extract_wave_features(args);
        features = args.features;
        %Do spikesorting on large spikes
        num_spikes = args.detection.large_amplitude_subset.num_spikes;
        large_ids = ones(num_spikes,1); %default to calling them all one cluster
        large_clustering = false;
        if features.large.num_spikes > args.params.spc.min_cluster
            sorterInput = args;
            sorterInput.spc.input_features = args.features.large.sortable_features;
            sorterOutput = spikesort(sorterInput);
            large_clustering = sorterOutput.clustering;
            large_ids = large_clustering.id;
        end
        
        %Do spikesorting on regular spikes
        num_spikes = args.detection.regular_amplitude_subset.num_spikes;
        reg_ids = ones(num_spikes,1); %default to calling them one cluster
        reg_clustering = false;
        if features.regular.num_spikes > args.params.spc.min_cluster
            sorterInput = args;
            sorterInput.spc.input_features = args.features.regular.sortable_features;
            sorterOutput = spikesort(sorterInput);
            reg_clustering = sorterOutput.clustering;
            reg_ids = reg_clustering.id;
        end
        
        num_spikes = args.detection.all_spikes.num_spikes;
        all_ids = zeros(num_spikes,1);
        large_spike_index = args.detection.all_spikes.class==2;
        all_ids(large_spike_index) = large_ids;
        max_large_id = max(large_ids);
        reg_spike_index=args.detection.all_spikes.class==1;
        if ~isempty(max_large_id)
            reg_ids(reg_ids>0) = reg_ids(reg_ids>0)+max_large_id;
        end
        all_ids(reg_spike_index) = reg_ids;
        
        args.sorted.id = all_ids;
        args.sorted.amplitude_class= args.detection.all_spikes.class;
        args.clustering.regular = reg_clustering;
        args.clustering.large = large_clustering;
        % exclude outliers
        args = exclude_outliers(args);
        
        %create the output data structure
        args = create_sort_results(args);
        
        args.data='finished processing'; % clear the memory of the original large time series
        args.time = datestr(datetime('now'));
        
        
        f.results{end+1}=args;        
        
    end
    
    if isfield(f,'savetofile')
        tofile = f.savetofile;
        save(tofile,'f');
    end
    
    
    
    output=f;
end