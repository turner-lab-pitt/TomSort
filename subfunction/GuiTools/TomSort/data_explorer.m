function output = data_explorer( argstruct, varargin )
if length(varargin)==1
    figtitle = varargin{1};
else
    figtitle = 'Data explorer';
end

%% Create index that selects which spikes to view
if isfield(argstruct.params, 'subset')
    subset = argstruct.params.subset;
else
    subset = true(size(argstruct.detection.all_spikes.wave_shapes,1),1);
end

%% If an error occurred, set all spikes to class 0 (unsorted)
if isfield(argstruct.clustering,'error') && argstruct.clustering.error>0
    if ~isfield(argstruct.clustering,'id')
        argstruct.clustering.id=ones(size(subset))*-2;
    end
end


if ~isfield(argstruct.clustering,'large_amplitude_subset')    
    argstruct.clustering.large_amplitude_subset.id=[];
end
if ~isfield(argstruct.clustering,'regular_amplitude_subset')    
    argstruct.clustering.regular_amplitude_subset.id=[];
end
% if ~isfield(argstruct.spc,'large_amplitude_subset')    
%     argstruct.spc.large_amplitude_subset.input_features=[];
% end
% if ~isfield(argstruct.spc,'regular_amplitude_subset')    
%     argstruct.spc.regular_amplitude_subset.input_features=[];
% end
argstruct.spc.large_amplitude_subset.input_features = argstruct.features.large.features;
argstruct.spc.regular_amplitude_subset.input_features = argstruct.features.regular.features;
LAS = argstruct.clustering.large_amplitude_subset;
RAS = argstruct.clustering.regular_amplitude_subset;

%% Convert back to original form (clustering struct contains all data)
argstruct.clustering = RAS;
argstruct.clustering.large_amplitude_subset = LAS;

clusters=argstruct.sort_results.spike_class(:);
spc_clusters=clusters;
clusters=[clusters; argstruct.clustering.large_amplitude_subset.id(:)+max(clusters)];
spc_clusters=clusters;
clusters=[clusters; -1*ones(argstruct.detection.noise_exclusion.num_spikes, 1)];
if isfield(argstruct,'spc')
    argstruct.spc.input_features = [argstruct.spc.regular_amplitude_subset.input_features;...
        argstruct.spc.large_amplitude_subset.input_features];
    wave_features = argstruct.spc.input_features;
else
    wave_features=[];
end

wave_shapes=argstruct.detection.all_spikes.wave_shapes(subset,:,1);
large_shapes = argstruct.detection.large_amplitude_subset.wave_shapes(:,:,1);
wave_shapes = [wave_shapes;large_shapes;argstruct.detection.noise_exclusion.wave_shapes(:,:,1)];

argstruct.detection.wave_index=argstruct.detection.all_spikes.wave_index(:);
wave_times=[argstruct.detection.wave_index(subset);...
    argstruct.detection.large_amplitude_subset.wave_index(:);...
    argstruct.detection.noise_exclusion.wave_index(:)];
wave_times=wave_times/argstruct.params.sr*1000; %in ms

% cluster_expand=argstruct.clustering.expanded;
peak_amplitude=wave_shapes(:,argstruct.params.sorting.points_before_peak);

% spike_filter = argstruct.clustering.mahal_distance < 25 &...
%                argstruct.clustering.multiple_possibile==false;
           
% spike_filter=ones(size(spike_filter));

% clusters=cluster_expand .* double(spike_filter);

numspikes = size(wave_features,1);

% colors={'k','b','r','g','c','y','m','b','r','g','c','y','m','b','r','g','c','y','m','b','r','g','c','y','m'};

cluster_nums = unique(clusters);
if cluster_nums(1)~=0, cluster_nums=[0; cluster_nums(:)]'; end %make sure a "no cluster" category is present
num_features = size(wave_features,2);
num_feature_plots=num_features-1;
num_clusters=length(cluster_nums);

f = figure('name',figtitle);


tgroup=uitabgroup('parent',f);
tab1=uitab('Title','Channel Summary','Parent',tgroup);
tab2=uitab('Title','Wave features','Parent',tgroup);
% tab3=uitab('Title','First pass','Parent',tgroup);
panel=uipanel(tab1,'Position',[1 1 98 98]/100, 'Units','Normalized');



for i=1:num_clusters
    % ***************outlier exclusion **************
    cl = cluster_nums(i);
    if cl<=0, continue; end
    cli=find(clusters==cl);
    
    M = argstruct.spc.input_features(cli,:);
    L = log(sum(bsxfun(@rdivide,bsxfun(@minus,M,mean(M)),std(M)) .^2, 2));
    L = L-median(L);
    MAD = median(abs(L));
    outliers = abs(L)>(5*MAD);
    clusters(cli(outliers)) = 0;
    scli = find(spc_clusters==cl);
    spc_clusters(scli(outliers))=0;
end



detection_window=(1:argstruct.params.sorting.dimensionality_reduction_features)+argstruct.params.sorting.points_before_peak - argstruct.params.sorting.features_before_peak;

vmax=max(max(wave_shapes(clusters>0,detection_window)))*1.1;
vmin=min(min(wave_shapes(clusters>0,detection_window)))*1.1;
if isempty(vmax) || vmax-vmin==0
    vmax=max(wave_shapes(:))*1.1;
    vmin=min(wave_shapes(:))*1.1;
end
vbinedges=linspace(vmin,vmax,100);
isibinwidth=0.5;
isiedges=0:isibinwidth:50.5;
w=100/sum(cluster_nums>=0);

for i=1:num_clusters
    cl = cluster_nums(i);
    
    if cl<0
        continue;
    elseif cl==0
        cli = find(clusters==0|clusters==-2);
        clr = find(clusters==-1);
        xpos=100-w;
    else
        cn = find(unique(cluster_nums(cluster_nums>0))==cl);
        xpos=(cn-1)*w;
        cli=find(clusters==cl);        
    end
    aWaves=axes('parent',panel,'position',[xpos 75 w 20]/100,'units','normalized');
    hold(aWaves,'on');    
    aAmp=axes('parent',panel,'position',[xpos 35 w 20]/100,'units','normalized');
    hold(aAmp,'on');
    aISI=axes('parent',panel,'position',[xpos 55 w 20]/100,'units','normalized');
    hold(aISI,'on');
   
    wsX = repmat([1:size(wave_shapes,2) nan],size(wave_shapes,1),1);
    wsY = [wave_shapes nan(size(wave_shapes,1),1)];
    if cl>0, 
        set(aWaves,'ColorOrderIndex',cn);
    end
    lh = plot(aWaves,reshape(wsX(cli,:)',[],1),reshape(wsY(cli,:)',[],1),'linestyle','-');
    if cl==0
        set(lh,'color','k');
        lh = plot(aWaves,reshape(wsX(clr,:)',[],1),reshape(wsY(clr,:)',[],1),'linestyle','-');
        set(lh,'color',[0.7 0.7 0.7]);
        plot(aWaves, mean(wave_shapes(clr,:)),':k','linewidth',2)
    end
    
    plot(aWaves, mean(wave_shapes(cli,:)),'-k','linewidth',2)
    if cl>0
        title(aWaves,[num2str(cl) ': ' num2str(length(cli))]);
    elseif cl==0
        title(aWaves,sprintf('%d unsorted; %d rejected',sum(clusters==0|clusters==-2),sum(clusters==-1)));
    end
    if cl==0, cli=union(cli,clr); end
    
    
    [pa,pa_edges] = histcounts(peak_amplitude(cli),100);
    bar(aAmp,pa_edges(1:end-1)+0.5*diff(pa_edges(1:2)),pa,1,'EdgeColor','none','FaceColor','k'); 
    isi=diff(wave_times(cli));
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

clusters=spc_clusters;
cluster_nums = unique(clusters);
if isempty(cluster_nums)||cluster_nums(1)~=0, cluster_nums=[0; cluster_nums(:)]'; end %make sure a "no cluster" category is present
num_features = size(wave_features,2);
num_feature_plots=num_features-1;
num_clusters=length(cluster_nums);

for i=1:num_features
    for j=num_features:-1:i+1
        aWaves=axes('position',[((i-1)/num_feature_plots) (1-(j-1)/num_feature_plots) 1/num_feature_plots 1/num_feature_plots],'units','normalized','parent',tab2);
        hold on        
        for k=1:length(cluster_nums)
            cn=cluster_nums(k);
            if k==2, set(aWaves,'ColorOrderIndex',1); end
            lh = plot(wave_features(clusters==cn,i), wave_features(clusters==cn,j),'.','markersize',0.5);
            if cn==0, set(lh,'color','k'); end
            
        end
        set(aWaves,'ButtonDownFcn', @link, 'xtick',[],'ytick',[],'xticklabels','','yticklabels','' )
        set(get(aWaves,'children'),'HitTest','off')
                
    end
end

enlarged_axis = axes('position',[ceil(num_feature_plots/3) 2*ceil(num_feature_plots/3) floor(num_feature_plots/3) floor(num_feature_plots/3)]/num_feature_plots, 'units', 'normalized','parent',tab2);
threeD_axis = axes('position',[2*ceil(num_feature_plots/3) ceil(num_feature_plots/3) floor(num_feature_plots/3) 2*floor(num_feature_plots/3)]/num_feature_plots, 'units', 'normalized','parent',tab2);


function link(src,~)
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

output=[];

    function draw_waveshapes(~,~,ax)
        replace_proportion=1;
        ud=get(ax,'userdata');
        lines=findobj(ax,'type','line');
        current_index=ud(end,end) + floor(replace_proportion*length(lines));
        nrows=size(ud,1)-1;
        for li=1:length(lines)
            if current_index>nrows, current_index=1;end
            set(lines(li),'ydata',ud(current_index,:))
            current_index=current_index+1;
        end
        ud(end,end)=current_index;
        set(ax,'userdata',ud);
%         disp(['test' num2str(current_index)])
    end
            
end

