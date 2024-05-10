function output = recluster(args,subset)
        
    args.params.subset=subset;
    args.detection.all_spikes.class = double(subset);
    args = extract_wave_features(args);
    args.spc.input_features = args.features.regular.sortable_features;
    args = do_permutation_clustering(args);
    
%     args = find_clusters(args);
    % args = run_spc(args); do this recursively within find_clusters
    
    
    output=args;
end