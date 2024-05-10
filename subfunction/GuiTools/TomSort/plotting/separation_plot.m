function separation_plot(C, id1, id2, target)


i1 = C.sort_results.spike_class==id1;
i2 = C.sort_results.spike_class==id2;
both = i1|i2;

F = C.features.regular.features(both,:);
F1 = C.features.regular.features(i1,:);
F2 = C.features.regular.features(i2,:);
ids = C.sort_results.spike_class(both);


% Build equal-sized populations of spike features for the two classes of
% spikes
RD1 = zeros(5000,10);
RD2 = RD1;
if size(F1,1) >= size(RD1,1)
    RD1 = F1(randsample(size(F1,1), size(RD1,1)),:);
else
    for jj=1:size(F1,2)
        [E,X]=ecdf(F1(:,jj));
        R=rand(size(RD1,1),1);
        S=sum(bsxfun(@gt,E,R'));
        RD1(:,jj)=X(S);
        [E,X]=ecdf(F2(:,jj));
        R=rand(size(RD1,1),1);
        S=sum(bsxfun(@gt,E,R'));
        RD2(:,jj)=X(S);
    end
end
if size(F2,1) >= size(RD2,1)
    RD2 = F2(randsample(size(F2,1), size(RD2,1)),:);
else
    for jj=1:size(F2,2)
        [E,X]=ecdf(F2(:,jj));
        R=rand(size(RD2,1),1);
        S=sum(bsxfun(@gt,E,R'));
        RD2(:,jj)=X(S);
    end
end

RD = [RD1; RD2];

% Project the points onto an axis defined by the cluster means in feature
% space
ax = mean(F2)-mean(F1);
pr1 = F1 * ax' / norm(ax); %projection of population 1
pr2 = F2 * ax' / norm(ax); %projection of population 2
pr = [pr1(:);pr2(:)]; %concatenate
pr_min=min(pr);
pr_max=max(pr);
idx=linspace(pr_min,pr_max,1000);
cdf1 = sum(bsxfun(@lt,pr1,idx)) / length(pr1);
cdf2 = sum(bsxfun(@gt,pr2,idx)) / length(pr2);
auc = trapz(cdf2,cdf1);
cmin=min(abs(cdf2 - (1-cdf1)));
mi=median(find(abs(cdf2 - (1-cdf1)) == cmin));

opt_thr = interp1(1:length(idx),idx,mi);

dprime = (mean(pr2)-mean(pr1))/sqrt(0.5*(std(pr2)+std(pr1)));

total_percent_misses = sum(pr1>opt_thr)/length(pr1)*100 + sum(pr2<opt_thr)/length(pr2)*100;

do_svm = total_percent_misses > 1 && true;

if do_svm

    SVMModel = fitcsvm(RD,[ones(length(RD1),1);ones(length(RD1),1)*-1],'KernelFunction','linear','Standardize',true);

    G = bsxfun(@rdivide,bsxfun(@minus,RD,SVMModel.Mu),SVMModel.Sigma) * SVMModel.Beta + SVMModel.Bias;

    G=G*-1;
    g1=G(1:length(RD1));
    g2=G(length(RD1)+1:end);
    idx=linspace(min(G),max(G),1000);
    gcdf1=sum(bsxfun(@lt,g1,idx))/length(g1);
    gcdf2=sum(bsxfun(@lt,g2,idx))/length(g2);
    
    gauc = trapz(gcdf2,gcdf1);

end


if nargin<4 
    target = figure;
end
% subplot(1,2,2);
% hold on;
% plot(cdf2,cdf1,'-k');
% if do_svm
%     plot(gcdf2,gcdf1,':k');
% else
%     gauc=0;
% end
% title(sprintf('AUC: %0.1f, %0.1f',auc*100,gauc*100));


subplot(1,2,1);
hold on;

if ~do_svm
    x=linspace(pr_min,pr_max,100);
    hc1=histcounts(pr1,x)/length(pr1);
    hc2=histcounts(pr2,x)/length(pr2);
    xp=(x(1:end-1)+diff(x))-opt_thr;
    plot(xp,hc1);
    plot(xp,hc2);
    plot([0 0],[0 1]*max([hc1 hc2]),'-k')
    %title(sprintf('Axis: %0.2f%% misses',sum(pr1>opt_thr)/length(pr1)*100));
    title('Projection onto axis between cluster centers');
else %do_svm==true
    
    x=linspace(min(G),max(G),100);
    xp2=x(1:end-1)+diff(x);
    % xp2 = xp2 * range(xp)/range(xp2);
    hcg1=histcounts(g1,x)/length(g1);
    hcg2=histcounts(g2,x)/length(g2);
    plot(xp2,hcg1);
    plot(xp2,hcg2);
    % plot([0 0],[0 1]*max([hcg1 hcg2]),'-k')
    %title(sprintf('SVM: %0.2f, %0.2f errors',sum(g1>0)/length(g1)*100,sum(g2<0)/length(g2)*100));
    title('Projection onto SVM axis')
end

%     subplot(2,2,4);
%     hold on;
%     plot(histcounts(RD(:,1)));
%     plot(histcounts(RD(:,2)));
%     title('Features for dip test');


% subplot(2,2,4);
% hold on;
% plot(F1(:,1),F1(:,2),'.','markersize',2);
% plot(F2(:,1),F2(:,2),'.','markersize',2);
% title(sprintf('H-score 1: %0.1f (%d spikes)\nH-score 2: %0.1f (%d spikes)\nCombined: %0.1f, ratio: %0.1f',...
%     max(h1),size(F1,1),max(h2),size(F2,1),max(h0),max(h_ratio)))

end