function [h,scores,alias_factor] = hscore(F,varargin)
antialias=false;
alias_factor=false;
all = false;
if ~isempty(varargin)
    if strcmp(varargin{1},'antialias')
        antialias=true;
        all=true;
    end
    if strcmp(varargin{1},'all')
        all=true;
    end
end

if isvector(F)
    all = true;
end
if ~all
    F(:,end)=[];
end
Forig=F;


% LIMIT = 100000;
% if size(F,1) >= LIMIT
%     F = F(randsample(size(F,1), LIMIT),:);
% end

% add some noise to de-alias



if antialias
    F = F + bsxfun(@times, randn(size(F)), std(F)/10);
end


x = zeros(1,size(F,2));
num_distributions = 3;
if size(F,1)<5000, num_distributions = 10; end
if size(F,1)<1000, num_distributions = 20; end
if size(F,1)>1

    hnorm=zeros(num_distributions,size(F,2));
    for ii=1:length(x)
        x(ii) = HartigansDipTest(F(:,ii));
        for jj=1:size(hnorm,1)
            hnorm(jj,ii)=HartigansDipTest(icdf('Normal',rand(size(F,1),1),mean(F(:,ii)),std(F(:,ii))));
        end

    %     hnorm(2,ii)=HartigansDipTest(icdf('Normal',rand(size(F,1),1),mean(F(:,ii)),std(F(:,ii))));
    %     hnorm(3,ii)=HartigansDipTest(icdf('Normal',rand(size(F,1),1),mean(F(:,ii)),std(F(:,ii))));
    end
    scores = x./mean(hnorm);

else
    scores = zeros(1,size(F,2));
end
if antialias
    F=Forig;
    x = zeros(1,size(F,2));
    num_distributions = 3;
    if size(F,1)<5000, num_distributions = 10; end
    if size(F,1)<1000, num_distributions = 20; end

    hnorm=zeros(num_distributions,size(F,2));
    for ii=1:length(x)
        x(ii) = HartigansDipTest(F(:,ii));
        for jj=1:size(hnorm,1)
            hnorm(jj,ii)=HartigansDipTest(icdf('Normal',rand(size(F,1),1),mean(F(:,ii)),std(F(:,ii))));
        end

    %     hnorm(2,ii)=HartigansDipTest(icdf('Normal',rand(size(F,1),1),mean(F(:,ii)),std(F(:,ii))));
    %     hnorm(3,ii)=HartigansDipTest(icdf('Normal',rand(size(F,1),1),mean(F(:,ii)),std(F(:,ii))));
    end
    alias_scores = x./mean(hnorm);
    alias_factor = alias_scores./scores;
end

h=max(scores);

end