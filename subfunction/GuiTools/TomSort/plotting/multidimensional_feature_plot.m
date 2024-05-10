function [ output_args ] = multidimensional_feature_plot( inputs )
%MULTIDIMENSIONAL_FEATURE_PLOT Summary of this function goes here
%   inputs: struct
%      clusters: vector of cluster ID
%      wave_features: array of features to plot
%      waveTab: ui component to plot regular size features
%      largeWaveTab: ui component to plot large features

clusters = inputs.clusters;
cluster_types = inputs.cluster_types;
wave_features = inputs.wave_features;
large_features = inputs.large_features;
waveTab = inputs.waveTab;
largewaveTab = inputs.largewaveTab;
axesButtonDownFcn = inputs.axesButtonDownFcn;

% Plot the features of the regular-amplitude spikes
        cluster_nums = unique(clusters);
        num_features = size(wave_features,2);
        num_feature_plots=num_features-1;
        %num_clusters=length(cluster_nums);

        enlarged_axis = axes('position',[ceil(num_feature_plots/3) 2*ceil(num_feature_plots/3) floor(num_feature_plots/3) floor(num_feature_plots/3)]/num_feature_plots, 'units', 'normalized','parent',waveTab);
        threeD_axis = axes('position',[2*ceil(num_feature_plots/3) ceil(num_feature_plots/3) floor(num_feature_plots/3) 2*floor(num_feature_plots/3)]/num_feature_plots, 'units', 'normalized','parent',waveTab);

        for i=1:num_features
            for j=num_features:-1:i+1
                aWaves=axes('position',[((i-1)/num_feature_plots) (1-(j-1)/num_feature_plots) 1/num_feature_plots 1/num_feature_plots],...
                    'units','normalized','parent',waveTab);
                hold on
                flag=true;
                for k=1:length(cluster_nums)
                    cn=cluster_nums(k);
                    if cn>0 && flag
                        set(aWaves,'ColorOrderIndex',1);
                        flag=false;
                    end
                    lh = plot(wave_features(clusters==cn,i), wave_features(clusters==cn,j),'.','markersize',0.5);
                    if cn==-1, set(lh,'color',[.7 .7 .7]); end
                    if cn==0, set(lh,'color','k'); end
                end
                set(aWaves,'ButtonDownFcn', {@link,enlarged_axis,threeD_axis,axesButtonDownFcn},...
                         'xtick',[],'ytick',[],'xticklabels','','yticklabels','' )
                set(get(aWaves,'children'),'HitTest','off')
%                 aWaves.ButtonDownFcn = axesButtonDownFcn;
                aWaves.UserData = struct('dims',[i,j],'type','regular');
            end
        end


        % If there was a separate class of large amplitude spikes, plot
        % those features
        if any(cluster_types==2)
            enlarged_axis = axes('position',[ceil(num_feature_plots/3) 2*ceil(num_feature_plots/3) floor(num_feature_plots/3) floor(num_feature_plots/3)]/num_feature_plots,...
                'units', 'normalized','parent',largewaveTab);
            threeD_axis = axes('position',[2*ceil(num_feature_plots/3) ceil(num_feature_plots/3) floor(num_feature_plots/3) 2*floor(num_feature_plots/3)]/num_feature_plots,...
                'units', 'normalized','parent',largewaveTab);

            for i=1:num_features
                for j=num_features:-1:i+1
                    aWaves=axes('position',[((i-1)/num_feature_plots) (1-(j-1)/num_feature_plots) 1/num_feature_plots 1/num_feature_plots],...
                        'units','normalized','parent',largewaveTab);
                    hold on        
                    flag=true;
                    for k=1:length(cluster_nums)
                        cn=cluster_nums(k);
                        if cn>0 && flag
                            set(aWaves,'ColorOrderIndex',1);
                            flag=false;
                        end
                        lh = plot(large_features(clusters==cn,i), large_features(clusters==cn,j),'.','markersize',0.5);
                        if cn==-1, set(lh,'color',[.7 .7 .7]); end
                        if cn==0, set(lh,'color','k'); end
                    end
                    set(aWaves,'ButtonDownFcn', {@link,enlarged_axis,threeD_axis,axesButtonDownFcn},...
                        'xtick',[],'ytick',[],'xticklabels','','yticklabels','' )
                    set(get(aWaves,'children'),'HitTest','off')
%                     aWaves.ButtonDownFcn = axesButtonDownFcn;
                    aWaves.UserData = struct('dims',[i,j],'type','large');
                end
            end

        end
        
end

function link(src,evt,enlarged_axis,threeD_axis,cb)
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
    cb(src,evt);
end