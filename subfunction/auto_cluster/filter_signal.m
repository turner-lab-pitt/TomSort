function output = filter_signal(argstruct)
output = argstruct;
x = double(argstruct.data);

params = argstruct.params;

sr = params.sr;

filters = params.filtering;

xf = zeros(length(filters),length(x));
for ii=1:length(filters)
    filter_def = filters(ii);
    filter_type = filter_def.type;
    tic;
    switch filter_type
        case 'none'
            xf(ii,:) = x;
        case 'elliptical'
            fmin = filter_def.fmin;
            fmax = filter_def.fmax;
            % HIGH-PASS FILTER OF THE DATA
            %Checks for the signal processing toolbox
            if exist('ellip','file')
                cutoffs = [fmin fmax]*2/sr;
                cutoffs(cutoffs<0)=0;
                cutoffs(cutoffs>=1) = [];
                if length(cutoffs)==2
                    [b,a] = ellip(2,0.1,40,cutoffs);
                else
                    [b,a] = ellip(2,0.1,40,cutoffs,'high');
                end
                xf(ii,:) = filtfilt(b, a, x);
            else
                error('detect_spikes requires the signal processing toolbox for this filter');
            end
        case 'wavelet'
            xf(ii,:) = wavelet_HP(x,'db4',5);
        otherwise
            if exist(filter_type,'file')
                try
                    Flt = eval(filter_type);
                    filtered = filter(Flt.Numerator, 1, x);
                    filtered(1:floor(length(Flt.Numerator)/2))=[];
                    xf(ii,1:length(filtered)) = filtered;
                catch Exception
                    error(Exception);
                end
            end
    end
    
end

output.timing.filtering=toc;
num_samples = min(length(x),10000);
output.filtering.timeseries.unfiltered = x(1:num_samples);
output.filtering.timeseries.filtered = xf(:,1:num_samples);
output.filtering.completed = true;
output.data = xf;

end