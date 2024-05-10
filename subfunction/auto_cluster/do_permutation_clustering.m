function output = do_permutation_clustering(argstruct)
output=argstruct;


M = argstruct.spc.input_features;

if ~isfield(output, 'clustering')
    output.clustering = struct;
end
% try
output.clustering = permutation_clustering(M,argstruct.params);
% catch ME
%     fprintf('permutation_clustering encountered an error: ');
%     error(ME.message);
    
%     output.clustering.error = true;
%     output.clustering.exception = ME;
   
% end
end