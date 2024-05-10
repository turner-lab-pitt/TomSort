function amplitude_sort_plot( ampinfo, target )
%AMPLITUDE_SORT_PLOT Summary of this function goes here
%   Detailed explanation goes here
X0 = ampinfo.X;
Y0 = ampinfo.Y;
Yinit = ampinfo.Fit';
Ypf = ampinfo.PeakFindingDist;
peak = ampinfo.Peak;
valley = ampinfo.Valley;
thr = ampinfo.InitialThr;
OFF = ampinfo.PctOff;

if nargin==1
    target = figure;    
end
axes(target);
hold on;
plot(Ypf,'displayname','Smoothed amplitude distrubution');
%plot(X0,Y0,'displayname','Peak amplitudes');
plot(X0,Yinit,'displayname','Exponential model');


% plot(X0(peak),Y0(peak),'ok');
% plot(X0(valley),Y0(valley),'or');


ylimit=get(gca,'ylim');
plot([1 1]*(ampinfo.Threshold1-thr),ylimit,'-k',...
    'displayname','Threshold')
if ~isnan(ampinfo.Threshold2)
    plot([1 1]*(ampinfo.Threshold2-thr),ylimit,'--k',...
    'displayname','Threshold for large-amplitude spikes')
end
%title(sprintf('Exponential model off by %0.2f%%',OFF));
ts = 'Threshold';
switch ampinfo.Threshold1Type
    case '25% of putative spikes'
        ts = 'Noise floor should be ~25% of total events'; 
    case '500 noise events'
        ts = 'Likely all noise: threshold to include ~500 events';
    case 'peak-based adjustment'
        ts = 'Threshold based on peak in amplitude distribution';
    case 'original'
        ts = 'Left part of distribution not well fit; threshold not adjusted';
    otherwise
        ts = sprintf('Threshold method: %s',ampinfo.Threshold1Type);
end


title(ts);
legend('show');

end

