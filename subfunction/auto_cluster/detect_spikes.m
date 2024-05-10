function output = detect_spikes(argstruct)
% Detect spikes with amplitude thresholding. Uses median estimation.

start_time=tic;
fprintf('Beginning detect_spikes... ');
output=argstruct;
x = argstruct.data;
params = argstruct.params;
plot_waves = false;

sr = params.sr;
window_before = params.sorting.points_before_peak;
window_after = params.sorting.points_after_peak;
detect = params.detection.voltagesign;
stdmin = params.detection.threshold;
stdmax = params.detection.artifacts;

% output.detection=struct;
% output.timing=struct;

xf = reshape(argstruct.data,1,[]);
xf_detection = xf(1,:);

%calculate a threshold on a second-by-second basis to 

samples_per_second = floor(sr);
overrun = mod(length(xf),samples_per_second);
if overrun > 0
    padding_length= samples_per_second - overrun;
    padding = xf(end:-1:length(xf)-padding_length+1);
    xf = [xf padding];    
end
R = reshape(xf,samples_per_second,[]);


noise_std_detect = reshape(repmat( median(abs(R),1)/0.6745, size(R,1), 1), 1, []);
noise_std_detect(length(x)+1 : end)=[];

%noise_std_detect = median(abs(xf_detection))/0.6745;
% noise_std_sorted = median(abs(xf))/0.6745;

thr = stdmin * noise_std_detect;        %thr for detection is based on detect settings.
thrmax = stdmax * noise_std_detect;     %thrmax for artifact removal is based on sorted settings.

output.detection.noise_std = noise_std_detect(1:samples_per_second:end);
output.detection.threshold = thr(1:samples_per_second:end);


%LOCATE SPIKE TIMES
switch detect
    case 'pos'
        xf_detect = xf_detection;
    case 'neg'
        xf_detect = xf_detection * -1;        
    case 'both'
        xf_detect = abs(xf_detection);
end



%define threshold crossing as first sample with amp > threshold after one
%sample with amp < threshold
thr_cross = find(xf_detect>=thr & ...
                [nan xf_detect(1:end-1)]<thr );% & ...
%                 [nan nan xf_detect(1:end-2)]<thr ); %this would be for
%                 two points below threshold
%remove threshold crossings that were in the first or last few samples
thr_cross(thr_cross<=window_before | thr_cross>=length(x)-window_after) = [];
            
%take the absolute value of xf_detect in case it isn't abs already
xf_detect=abs(xf_detect);

%detect_window_ms defines number of milliseconds after each threshold
%crossing to look for the peak amplitude, for alignment purposes
detect_window_ms=0.5;
post_window=round(sr/1000*detect_window_ms);

fprintf('Iterating over threshold crossing events to find peaks...');
tic;
for i=1:length(thr_cross)
    cross=thr_cross(i);    
    pre_window=0;
    post_window_inner=post_window;
    while(true)
        window=cross-pre_window:cross+post_window_inner;
        wave=xf_detect(window);    
        [~,mi]=max(wave.*double(sign(xf(window))==sign(xf(cross))));

        if mi==1
            pre_window=pre_window+1;
        elseif mi==length(wave)
            post_window_inner=post_window_inner+1;
        else
            break;
        end
        
    end
    thr_cross(i)=thr_cross(i)+(mi-1)-pre_window;
end

%reject duplicate spike times - this can occur if noise at threhold caused
%duplicate crossings which were both aligned to the same peak within the
%detection window.
thr_cross=unique(thr_cross);
output.timing.find_peaks=toc;
fprintf('finished in %0.1f seconds.\n', output.timing.find_peaks);

%suppress short-ISI detections with opposing voltages, as these are very
%likely to represent large spikes and afterhyperpolarizations. Also if two
%threshold crossings => two peaks are detected in a very short window and
%the voltage does not return to <thr/2 between, merge the peaks.

tic;


short_isi = find(diff(thr_cross)<window_after)+1;
fprintf('Applying duplicate-detection algorithm to %d pairs of spikes...',length(short_isi));
%work backwards
for s=length(short_isi):-1:1
    i=short_isi(s);
    crossing=thr_cross(i);
    if(isnan(crossing)), continue; end %already rejected this crossing
    dup_candidate_indices= i-1:-1:max(i-window_after, 1);
    dup_candidates=thr_cross(dup_candidate_indices);
    dup_candidates(dup_candidates<(crossing-window_after) | isnan(dup_candidates) )=[];
%     dup_candidates=fliplr(thr_cross(thr_cross<crossing&thr_cross>));
    for j=1:length(dup_candidates)
        % if a negative and positive peak are within the same detection
        % window, try to determine which one represents the true peak and
        % reject the other
        if sign(xf_detection(crossing)) ~= sign(xf_detection(dup_candidates(j)))
            if xf_detect(crossing) >= xf_detect(dup_candidates(j))
                thr_cross(thr_cross==dup_candidates(j))=nan;
            else
                thr_cross(i)=nan;
                break;
            end
        % if two negative or two positive threhold crossings are detected within
        % the same detection window, and the voltage does not return below
        % threshold/2 between them, reject one peak and set the
        % detected peak time to the mean of the two
        elseif min(xf_detect(dup_candidates(j):crossing)) > thr(crossing)/2
            thr_cross(thr_cross==dup_candidates(j))=nan;
            thr_cross(i)=round(mean([crossing dup_candidates(j)]));
        end        
    end
end
output.timing.short_isi=toc;
fprintf('finished in %0.1f seconds. Rejected %d peaks.\n', output.timing.short_isi, sum(isnan(thr_cross)));

spike_indices=thr_cross(~isnan(thr_cross));
wave_indices=bsxfun(@plus,spike_indices', -1*(window_before+1):(window_after+2));
%if a spike happened too early in the data stream to have the full
%peri-event window, exclude it.
spike_indices(wave_indices(:,1)<=0)=[];
wave_indices(wave_indices(:,1)<=0,:)=[];
% spikes = xf_detection(wave_indices);

mds = permute(xf,[3 2 1]);
mdspikes=mds(bsxfun(@plus,repmat(wave_indices,1,1,size(mds,3)), permute((1:size(mds,3))-1,[3 1 2])*size(mds,2)));
spikes = mdspikes;
switch params.interpolate
    case 'n'
        spikes(:,end-1:end,:)=[];       %eliminates borders that were introduced for interpolation
        spikes(:,1:2,:)=[];
    case 'y'
        %Does interpolation
        error('Interpolation not supported');
%         tic;
%         fprintf('Interpolate spike shapes for optimal alignment...');
%         interpolation_factor=10;
%         ls = size(spikes,2);
%         alignment=window_before+1;
%         interpolated = spline(1:ls, spikes, 1:(1/interpolation_factor):ls);
%         si=sign(diff(interpolated,1,2));
%         [zd,zi]=min(abs(bsxfun(@times,2:size(interpolated,2)-1,double(si(:,1:end-1)~=si(:,2:end))) - alignment*interpolation_factor),[],2);
%         zi=zi+1;
%         zi(zd>interpolation_factor*1.5)=alignment*interpolation_factor;%if the algorithm fails, leave it as it was
%         interp_index=bsxfun(@plus,(-window_before+1:window_after)*interpolation_factor,zi);
%         spikes=interpolated( bsxfun(@plus,(1:size(interpolated,1))', ((interp_index-1)*size(interp_index,1)) ));
%         output.timing.interpolation = toc;
%        
%         fprintf('finished in %0.1f seconds.\n', output.timing.interpolation);
end
% spikes = mdspikes;
%remove too-large spikes as likely artifactual
amplitude=spikes(:,window_before,1);
amp_above_thr = abs(amplitude) - abs(thr(spike_indices))';
artifacts = abs(amplitude) > thrmax(spike_indices)';

output.detection.artifacts=struct;
output.detection.artifacts.wave_orig_index = find(artifacts);
output.detection.artifacts.wave_shapes=spikes(artifacts,:,:);

spike_indices(artifacts)=[];
spikes(artifacts,:,:)=[];
amplitude(artifacts)=[];
amp_above_thr(artifacts)=[];

%amplitude_sort: uses a smarter algorithm to pick a threshold, excluding
%noise while preserving non-noise events. Gives each event a class: 0 =
%noise, 1=regular spike, 2=large spike
amp_info = amplitude_sort(amp_above_thr);
if ~params.detection.amplitude_sort
    amp_info.class = ones(length(spike_indices),1);
end
class0 = amp_info.class==0;
class1 = amp_info.class==1;
class2 = amp_info.class==2;

output.detection.amplitude_info=amp_info;
output.detection.pre_amplitude_sort.threshold=thr;

output.detection.regular_amplitude_subset.timeseries_length = length(x);
output.detection.regular_amplitude_subset.wave_index = spike_indices(class1);
output.detection.regular_amplitude_subset.wave_shapes=spikes(class1,:,:);
output.detection.regular_amplitude_subset.wave_amplitude = amplitude(class1);
output.detection.regular_amplitude_subset.wave_peak_sign = sign(amplitude(class1));
output.detection.regular_amplitude_subset.num_spikes = length(spike_indices(class1));
ns1=sum(class1);

output.detection.noise_exclusion.wave_index = spike_indices(class0);
output.detection.noise_exclusion.wave_shapes=spikes(class0,:,:);
output.detection.noise_exclusion.wave_amplitude = amplitude(class0);
output.detection.noise_exclusion.wave_peak_sign = sign(amplitude(class0));
output.detection.noise_exclusion.num_spikes = length(spike_indices(class0));
output.detection.noise_exclusion.timeseries_length = length(x);
ns0=sum(class0);

output.detection.large_amplitude_subset.wave_index = spike_indices(class2);
output.detection.large_amplitude_subset.wave_shapes=spikes(class2,:,:);
output.detection.large_amplitude_subset.wave_amplitude = amplitude(class2);
output.detection.large_amplitude_subset.wave_peak_sign = sign(amplitude(class2));
output.detection.large_amplitude_subset.num_spikes = length(spike_indices(class2));
output.detection.large_amplitude_subset.timeseries_length = length(x);
ns2 = sum(class2);

output.detection.all_spikes.wave_index = spike_indices(class1|class2);
output.detection.all_spikes.wave_shapes=spikes(class1|class2,:,:);
output.detection.all_spikes.wave_amplitude = amplitude(class1|class2);
output.detection.all_spikes.wave_peak_sign = sign(amplitude(class1|class2));
output.detection.all_spikes.num_spikes = length(spike_indices(class1|class2));
output.detection.all_spikes.class = amp_info.class(class1|class2);
output.detection.all_spikes.timeseries_length = length(x);

fprintf('Detected %d spikes (lg: %d/reg: %d/noise: %d), t=%0.2fs.\n',...
    length(spike_indices),ns2,ns1,ns0,toc(start_time));

if plot_waves
    figure;
    hold on;
    t = (1:length(xf))/params.sr;
    num_to_plot=100000;
    plot(t(1:num_to_plot),xf_detection(1:num_to_plot),'-b');
    plot(t(1:num_to_plot),xf(1:num_to_plot),'-r');
    set(gca,'xlim',[t(1) t(num_to_plot)]);
    plot([t(1) t(num_to_plot)], [thr thr], '--k');    
    plot([t(1) t(num_to_plot)], [-thr -thr], '--k');
    
    wi=spike_indices(spike_indices<=num_to_plot);
    wt=t(bsxfun(@plus,-19:44, wi'))';
    wv=spikes(spike_indices<num_to_plot,:)';
    wt(end+1,:)=nan;
    wv(end+1,:)=nan;
    plot(wt(:),wv(:),'-m','linewidth',2);
    
    tt=t([wi;wi]);
    tt(end+1,:)=nan;
    tv=repmat([-120;80;nan],1,length(wi));
    plot(tt(:),tv(:),'-k')
end

   
end