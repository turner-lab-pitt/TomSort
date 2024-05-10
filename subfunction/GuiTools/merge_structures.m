function S = merge_structures( S1, S2 )
%MERGE_STRUCTURES Summary of this function goes here
%   Detailed explanation goes here
S = repmat(S1,size(S2));


for j=1:length(S)
    fields = fieldnames(S2(j));
    for i = 1:length(fields),
        val = S2(j).(fields{i});
        if isstruct(val) && isfield(S(j), fields{i}) && isstruct(S(j).(fields{i}))
            S(j).(fields{i}) = merge_structures(S(j).(fields{i}), val);
        else
            S(j).(fields{i}) = S2(j).(fields{i});
        end
    end

end

