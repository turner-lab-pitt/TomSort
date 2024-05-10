function output = do_recursive_clustering(argstruct)
output=argstruct;

M = argstruct.spc.input_features;

output.spc.clusters = recursive_clustering(M,argstruct.params,1);

end