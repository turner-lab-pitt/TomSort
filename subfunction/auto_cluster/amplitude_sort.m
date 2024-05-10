function amplitude_info = amplitude_sort(amp_above_thr,make_plot)
if ~exist('make_plot','var')
    make_plot=false;
end
amplitude_info=struct;
amplitude_info.class = ones(size(amp_above_thr)); %initialize the spikes to all be class 1
amplitude_info.Threshold1 = NaN;
amplitude_info.Threshold2 = NaN;

original_num_spikes = length(amp_above_thr);

A = abs(amp_above_thr);
%thr = min(thr);
%T = thr;
thr = 0;
T = 0;
A = A(A>T) - T; %subtract minimum threshold voltage and ignore negative values
A=[-1*A(:);A(:)];%invert signal for smoother kernel smoothing;
% look for the different segments of the distribution
% initial descending part:

bandwidth = 2;
[Y0,X0]=ksdensity(A,0:1:max(A)*2,'width',bandwidth);

%remove the mirrored half of the distribution
Y0(X0<1)=[];
X0(X0<1)=[];
A(A<0)=[];

%normalize to area of 1
norm_factor = sum(Y0);
Y0=Y0/norm_factor; 

%number of spikes used for fitting the density function
num_spikes_in_Y0 = length(A);

peak_finding_signal = Y0;

while true
    %find peaks; if too many are found, increase the smoothing factor and
    %try again.
    [pks,locs,w,~] = findpeaks(peak_finding_signal);    
    if diff(peak_finding_signal(1:2)) >0 
        %reject first peak since the distribution started with a positive slope  
        pks(1)=[];
        locs(1)=[];
        w(1)=[];
    end
    peak_inds = locs( pks.*w*length(A) > 30 );
    
    [~,~,~,valley_inds]=extrema(peak_finding_signal);
    valley_inds=sort(valley_inds);

    if length(peak_inds)<4
        %we've found few enough peaks; break out of the loop
        break
    else
        %too many peaks. Smooth more.
        bandwidth = bandwidth * 2;
        peak_finding_signal = ksdensity(A,0:1:max(A)*2,'width',bandwidth);        
    end
    
end

peak=[];
valley=[];

%determine which points to use for the exponential curve fit.
%Ignore the first couple of points; the smoothing procedure blunts the
%empirical distribution
if isempty(peak_inds)
    % no peaks found; fit all of the data with an exponential
    index_to_fit = 3:length(X0);    
else
    % found some peaks; fit only the data before the first peak to see if
    % that looks like a true noise type of distribution
    peak = peak_inds(1);
    valley = max(valley_inds(valley_inds<peak));
    %valley is the index of the valley before the first peak
    %Fit only the data up to the valley
    index_to_fit = 3:valley;    
end
%Fit an exponential to the selected portion of the data
F = fit(X0(index_to_fit)',Y0(index_to_fit)','exp1');
Yinit = F(X0)';
    
    
% I is the index of the valley or the index after which the fit distribution
% drops below 0.01% of the initial value.
I = min( [ find(Yinit>0.0001*Yinit(1),1,'last'), valley] );

% Calculate the error of the exponential fit. If it is small enough, this
% is considered to be noise and some can be excluded
PctOff = 100*sum(max(Y0(3:I)-Yinit(3:I),zeros(size(3:I))))/sum(Y0(1:I));
PctUnder = 100*sum(max(Yinit(3:I)-Y0(3:I),zeros(size(3:I))))/sum(Yinit(1:I));


if  max([PctOff PctUnder/5]) < 5
    % error < 5%: the left part of the distribution is essentially noise
    % exclude some noise by increasing the effective threshold
    
    % 
    % the goal is to estimate the number of non-noise events, and use that
    % to inform how much noise to keep.
        
    
    
    %Estimate true event count as the deviation of the exponential fit from
    %the empirical distribution.
    estimated_real_events = sum(max(Y0-Yinit, 0));
    
    %Figure out where to set a threshold to get all the estimated real
    %events plus min(500, 25% of real events,peak-based estimate) so the noise will cluster
    %well in the next steps
    min_noise_spikes=500;
    prctile_threshold = find(cumsum(Y0)<(1-estimated_real_events*1.25),1,'last');
    abs_threshold = find(cumsum(Y0)<(1-estimated_real_events-min_noise_spikes/num_spikes_in_Y0),1,'last');
    peak_based_threshold = length(Y0);
    
    % If a peak was detected in the distribution, determine where a good
    % cutoff would be based on the peak location
    if ~isempty(peak_inds)
        
        peak = peak_inds(1);
        % first look for a natural break - a trough in the density function
        % that dips effectively to zero.
        split_ind_last = find(Y0(1:peak)<1e-8, 1, 'last');
        split_ind_first = find(Y0(1:peak)<1e-8, 1, 'first');
        if ~isempty(split_ind_last)
            % Found a natural break.
            peak_based_threshold = ( X0(split_ind_last) + X0(split_ind_first) )/2;
%             amplitude_info.class(abs(amp)>split_amp) = 1;
%             amplitude_info.thr = split_amp;
        else
            % no natural break found.
            % One option: The threshold (as a noise estimate)
            % provides information about how much spread is expected in a
            % real distribution of spike amplitudes. Set the new threshold
            % at peak-threshold to get all the spikes
            % Second option: use the location of the valley
            % Final option: use the point where the exponential fit drops
            % to effectively zero.
            % Strategy: use the lower of these three cutoffs.
            thr_estimate = max( (X0(peak) - thr), 0);
            valley_estimate = X0(valley);
            exp_estimate = X0(find(Yinit < 1e-4,1,'first'));
            estimated_amp = min([thr_estimate, valley_estimate, exp_estimate]);
            
            %lower the estimate a bit so some noise remains but don't lower
            %it beyond the original threshold
            peak_based_threshold = max(estimated_amp-(0.5*thr), 0); 
            % for clustering and comparing putative units to noise floor 
            
            %but don't let it go lower than the initial threshold
%             split_amp = max([split_amp thr]);
            
            %Add spikes back in as class 1
%             amplitude_info.class(abs(amp)>split_amp)=1;
%             amplitude_info.thr=split_amp;
                    

        end
    end
    
    %set the adjusted threshold as the lowest of these three estimates
    [threshold_adjustment,adjustment_type] = min([prctile_threshold, abs_threshold, peak_based_threshold]);
    %don't lower it beyond the original
    final_threshold = thr + max(threshold_adjustment, 0);
    
    % Anything below this new threshold is noise: set to class 0
    amplitude_info.class(abs(amp_above_thr)<final_threshold)=0;
    amplitude_info.Threshold1=final_threshold;
    switch adjustment_type
        case 1
            amplitude_info.Threshold1Type = '25% of putative spikes';
        case 2
            amplitude_info.Threshold1Type = '500 noise events';
        case 3
            amplitude_info.Threshold1Type = 'peak-based adjustment';
    end
    
else 
    % The exponential was not a good fit. Don't touch the original
    % threshold.
    amplitude_info.Threshold1 = thr;
    amplitude_info.Threshold1Type = 'original';
end

% before returning, see if there are large-amplitude spikes that 
% should be treated separately
if ~isempty(peak_inds)
    %figure out where to draw the cutoff.
    %first, look for the lowest-amplitude peak that is separable from
    %the noise
    first_low_point = find(Y0<1e-8, 1, 'first');
    peak = peak_inds(peak_inds > first_low_point);
    if ~isempty(peak)
        % A clean split in amplitudes was detected.
        % find the best location to split, and assign those spikes to the
        % large spike subset.
        
        peak = peak(1);
        %now find the cutoff to the left of the peak
        split_ind_last = find(Y0(1:peak)<1e-8, 1, 'last');
        split_ind_first = find(Y0(1:peak)<1e-8, 1, 'first');

        split_amp = ( X0(split_ind_last) + X0(split_ind_first) )/2 + thr;
        amplitude_info.class(abs(amp_above_thr)>split_amp) = 2;
        amplitude_info.Threshold2 = split_amp;
    %else: no spikes should be assigned to the large group
    end
    
end
    
amplitude_info.X=X0;
amplitude_info.Y=Y0;
amplitude_info.Fit=F(X0);
amplitude_info.PctOff=PctOff;
amplitude_info.PctUnder=PctUnder;
amplitude_info.InitialThr = T;
amplitude_info.Peak = peak;
amplitude_info.Valley = valley;
amplitude_info.PeakFindingDist = peak_finding_signal;
    
if make_plot
    amplitude_sort_plot(amplitude_info);
end

end