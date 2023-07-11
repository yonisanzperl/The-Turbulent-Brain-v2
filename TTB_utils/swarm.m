function [out,stat_out] = swarm(xx, gg,options)
% swarm Function to get a pretty bee swarm plot 
%   swarm(xx, gg)
%       xx is a cell array with the data, each column is a group.
%       gg is a cell array with the name of each group (each column of xx).
%
%   swarm(xx, gg, 'PARAM1', val1, 'PARAM2', val2, ...) specifies optional
%   parameter name/value pairs.
%       'tlt' Is a character string and specifies the title of the plot.
%   
%       'stat_test' Is a character string specifying the type of
%       statistical test to apply to the data. Default is ranksum 
%
%       'fdr' is a boolean to apply multiple comparison correction if N > 2
%       with the Benjamini-Hochberg method. The default is true.
% 
%       'printPvals' is a boolean to display pvalues. The default is false.
%       
%       'overlay_style' is a character string and specifies what to overlay
%       on the beeswarm plot. 
%       (false default | 'box' | 'sd' | 'ci' | 'median' | 'boxplot')
%
% elvira.delagua@upf.edu

    arguments
       xx cell
       gg cell
       options.tlt string = ''
       options.fdr logical = true 
       options.printPvals logical = false
       options.overlay_style string = ''
       options.dot_size double = 0.0;
       options.alpha double = 0.0;
       options.stat_test string = 'ranksum'
       options.name string = ''
       options.display_on_plot logical = true

    end
    
    % if xx is not a column cell array, collapse it into one
    if size(xx, 1) > 1
        temp = xx;
        clear xx;
        for i = 1:size(temp, 2)
            xx{1, i} = [temp{:, i}];
        end
    end
    
    % Put data (xx) into a column array x and groups gg into grouping
    % variables g paired with each data point in x
    x = []; g = [];
    for i=1:length(xx)
        x = vertcat(x, xx{i}(:));
        g = vertcat(g, repmat({gg{i}}, length(xx{i}),1));
    end
    
    % If we are printing pValues, write the title of the plot as a header
    % first:
    if options.printPvals
        fprintf('\n\n'); 
        disp(strcat('          <strong>', upper(options.tlt), '</strong>'))
    end
    fID=fopen(options.name,'w');
    % 
    p=[]; groups = {}; cont=1;
    for i=1:(length(xx)-1)
        for j = (i+1):length(xx)
            
            if strcmp(options.stat_test, 'permutation')
                stat = permutationTest(xx{i}, xx{j}, 10000);
            end
            if strcmp(options.stat_test, 'ranksum')
                stat = ranksum(xx{i}, xx{j});
            end
             if strcmp(options.stat_test, 'signrank')
                stat = signrank(xx{i}, xx{j});
             end
             if strcmp(options.stat_test, 'ttest')
                [~,stat] = ttest(xx{i}, xx{j});
             end
             if strcmp(options.stat_test, 'ttest2')
                [~,stat] = ttest2(xx{i}, xx{j});
            end           
            
            
                stat_out(i,j)=stat;
            if options.printPvals           
                fprintf('P-val of %s test, %s vs. %s: %s \n', options.stat_test, gg{i}, gg{j}, num2str(stat));
                fprintf(fID,'%d) P-val of %s test, %s vs. %s: %s \n',cont, options.stat_test, gg{i}, gg{j}, num2str(stat));

            end
            
            p = horzcat(p, stat); % array of p-values
            groups{end+1} = {gg{i}, gg{j}}; % cell array with group names of each pair-wise comparison
            cont=cont+1;
        end
    end
    
    % Optionally correct for multiple comparisons
    if options.fdr && length(groups) > 1
        ifdr = sort(FDR_benjHoch(p, 0.05));
        p = p(ifdr);
        groups = groups(ifdr);
        fprintf(fID,'Significant comparison after FDR correction %i \n', ifdr);
    end

    % If the overlay options are boxplot or compact boxplot and alpha is in
    % the default option, make beeswarm plot transparency 0.5
    if (strcmp(options.overlay_style, 'boxplot') || strcmp(options.overlay_style, 'compact boxplot')) ...
            && (options.alpha == 0)
        options.alpha = 0.5;
    elseif options.alpha == 0 % else make beeswarm plot transparency 0.25
        options.alpha = 0.25; 
    end

    % If dot_size has been left default, do not pass it as an argument to
    % the beeswarm function. The beeswarm function will choose an
    % appropiate size automatically
    if options.dot_size == 0
        beeswarm(celllbls2nums(g), x, corral_style='random', overlay_style=options.overlay_style, ...
            MarkerFaceAlpha=options.alpha); hold on;
    else
        beeswarm(celllbls2nums(g), x, corral_style='random', overlay_style=options.overlay_style, ...
            MarkerFaceAlpha=options.alpha, dot_size=options.dot_size); hold on;
    end

    % Set settings for each overlay_style
    if strcmp(options.overlay_style, 'boxplot')
        h = boxplot(x, g, 'symbol', ''); hold on;
        set(h,'LineWidth',0.75); set(h, 'Color', 'k');
    elseif strcmp(options.overlay_style, 'compact boxplot')
        h = boxplot(x, g, 'symbol', '', 'PlotStyle','compact', 'MedianStyle','line'); hold on;
        set(h,'LineWidth',1); set(h, 'Color', 'k');
     elseif strcmp(options.overlay_style, 'violin')
        h = violinplot(cell2mat(xx));%, g, 'symbol', '', 'PlotStyle','compact', 'MedianStyle','line'); hold on;
        %set(h,'LineWidth',1); set(h, 'Color', 'k');
    end

    xticks([1:1:length(unique(gg))]);
    xticklabels(gg);

    hold on;

    %--------------------------------------------------------------------------
  if options.display_on_plot
    % Now plot asterisk bars where p <= 0.05
    groups_sigstar = groups(p<=0.05); stats = p(p<=0.05);
    
    xlocs=nan(length(groups_sigstar),2); %matrix that will store the indices 
    xtl=get(gca,'XTickLabel'); 
    
    for ii=1:length(groups_sigstar)
        grp=groups_sigstar{ii};
    
        if isnumeric(grp)
            xlocs(ii,:)=grp; %Just store the indices if they're the right format already
    
        elseif iscell(grp) %Handle string pairs or string/index pairs
    
            if isstr(grp{1})
                a=strmatch(grp{1},xtl);
            elseif isnumeric(grp{1})
                a=grp{1};
            end
            if isstr(grp{2})
                b=strmatch(grp{2},xtl);
            elseif isnumeric(grp{2})
                b=grp{2};
            end
    
            xlocs(ii,:)=[a,b];
        end
    
        %Ensure that the first column is always smaller number than the second
        xlocs(ii,:)=sort(xlocs(ii,:));
    
    end
    
    %If there are any NaNs we have messed up. 
    if any(isnan(xlocs(:)))
        error('Some groups were not found')
    end
    
    % Sort sig bars from shortest to longest
    [~,ind]=sort(xlocs(:,2)-xlocs(:,1),'ascend');
    xlocs=xlocs(ind,:);groups_sigstar=groups_sigstar(ind);
    stats=stats(ind);
    
    %Add the sig bar lines and asterisks 
    holdstate=ishold;
    hold on
    
    H=ones(length(groups_sigstar),2); %The handles will be stored here
    
    y=ylim;
    yd=myRange(y)*0.05; %separate sig bars vertically by 5% 
    
    for ii=1:length(groups_sigstar)
        thisY=findMinY(xlocs(ii,:))+yd;
        H(ii,:)=makeSignificanceBar(xlocs(ii,:),thisY,stats(ii));
    end
    
    % Add downward ticks on the ends of each line
    yd=myRange(ylim)*0.01; %Ticks are 1% of the y axis range
    for ii=1:length(groups_sigstar)
        y=get(H(ii,1),'YData');
        y(1)=y(1)-yd;
        y(4)=y(4)-yd;   
        set(H(ii,1),'YData',y)
    end
  else
     holdstate=ishold;
  
  end
    
    
    %--------------------------------------------------------------------------

    title(options.tlt);
    set(gca,'fontname','times')  % Set it to times
    set(gca, 'FontSize', 12); hold on;
    
    % Return cell array with group label pairs and corresponding pvals
    if ~isempty(groups)
        for i = 1:length(groups)
            pair = groups{i};
            groups_out{i} = [pair{1} ' vs. ' pair{2}];
        end
        out = [groups_out; num2cell(p)];
    else
        out = {};
    end
    
    % Return hold state to whatever it was before we started
    if ~holdstate
        hold off
    elseif holdstate
        hold on
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Internal functions

function H=makeSignificanceBar(x,y,p)
    %makeSignificanceBar produces the bar and defines how many asterisks we get for a 
    %given p-value


    if p<=1E-3
        stars='***'; 
    elseif p<=1E-2
        stars='**';
    elseif p<=0.05
        stars='*';
    elseif isnan(p)
        stars='n.s.';
    else
        stars='';
    end
            
    x=repmat(x,2,1);
    y=repmat(y,4,1);

    H(1)=plot(x(:),y,'-k','LineWidth',1.5,'Tag','sigstar_bar');

    %Increase offset between line and text if we will print "n.s."
    %instead of a star. 
    if ~isnan(p)
        offset=0.005;
    else
        offset=0.02;
    end

    starY=mean(y)+myRange(ylim)*offset;
    H(2)=text(mean(x(:)),starY,stars,...
        'HorizontalAlignment','Center',...
        'BackGroundColor','none',...
        'Tag','sigstar_stars');

    Y=ylim;
    if Y(2)<starY
        ylim([Y(1),starY+myRange(Y)*0.05])
    end


end %close makeSignificanceBar



function Y=findMinY(x)
    % The significance bar needs to be plotted a reasonable distance above all the data points
    % found over a particular range of X values. So we need to find these data and calculat the 
    % the minimum y value needed to clear all the plotted data present over this given range of 
    % x values. 
    %
    % This version of the function is a fix from Evan Remington
    oldXLim = get(gca,'XLim');
    oldYLim = get(gca,'YLim');

    axis(gca,'tight')
    
    %increase range of x values by 0.1 to ensure correct y max is used
    x(1)=x(1)-0.1;
    x(2)=x(2)+0.1;
    
    set(gca,'xlim',x) %Matlab automatically re-tightens y-axis

    yLim = get(gca,'YLim'); %Now have max y value of all elements within range.
    Y = max(yLim);

    axis(gca,'normal')
    set(gca,'XLim',oldXLim,'YLim',oldYLim)

end %close findMinY


function rng=myRange(x)
    %replacement for stats toolbox range function
    rng = max(x) - min(x);
end %close myRange



function nums=celllbls2nums(g)
    % turn column cell array with labels to cell array with numbers
    gg = unique(g, 'stable');
    n = [1:1:length(gg)]'; 
    nums = zeros(size(g));
    for i = 1:length(gg)
        idx = strfind(g, gg{i});
        idx = find(not(cellfun('isempty', idx)));
        nums(idx) = n(i);
    end
end

function x = beeswarm(x,y,varargin)
%function xbee = beeswarm(x,y)
%
% Input arguments:
%   x               column vector of groups (only tested for integer)
%   y               column vector of data
%
% Optional input arguments:
%   sort_style      ('nosort' - default | 'up' | 'down' | 'fan' | 'rand' | 'square' | 'hex')
%   corral_style    ('none' default | 'gutter' | 'omit' | 'rand')
%   dot_size        relative. default=1
%   overlay_style   (false default | 'box' | 'sd' | 'ci' | 'median')
%   use_current_axes (false default | true)
%   colormap        (lines default | 'jet' | 'parula' | 'r' | Nx3 matrix of RGB values]
%
% Output arguments:
%   xbee            optimized layout positions
%
% Known Issues:
%       x locations depend on figure aspect ratio. resizing the figure window and rerunning may give different results
%       setting corral to 'none' still has a gutter when the width is large
%
% Usage example:
% 	x = round(rand(150,1)*5);
%   y = randn(150,1);
%   beeswarm(x,y,3,'sort_style','up','overlay_style','ci')
%
% % Ian Stevenson, CC-BY 2019

p = inputParser;
addRequired(p,'x')
addRequired(p,'y')
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addOptional(p,'sort_style','nosort')
addOptional(p,'corral_style','none')
addOptional(p,'dot_size',11/sqrt(length(x)),validScalarPosNum)
addOptional(p,'overlay_style',false)
addOptional(p,'use_current_axes',false)
addOptional(p,'colormap','lines')
addOptional(p,'MarkerFaceColor','')
addOptional(p,'MarkerFaceAlpha',.5)
addOptional(p,'MarkerEdgeColor','none')
parse(p,x,y,varargin{:});

% extra parameters
rwid = .05; % width of overlay box/dash

dcut=8; % spacing factor
nxloc=512; % resolution for optimization
chanwid = .5; % percent width of channel to use
yl = [min(y) max(y)]; % default y-limits
asp_rat = 1;
keep_hold = false;

% get aspect ratio for a figure window
if isfinite(p.Results.dot_size)
    if ~p.Results.use_current_axes
        % make new axes
        s=scatter(x,y);
        xl=[min(x)-.5 max(x)+.5];
    else
        xl=xlim();
    end
    yl=ylim();
    pasp_rat = get(gca,'PlotBoxAspectRatio');
    dasp_rat = get(gca,'DataAspectRatio');
    asp_rat = pasp_rat(1)/pasp_rat(2);
    
    % pix-scale
    pf = get(gcf,'Position');
    pa = get(gca,'Position');
    as = pf(3:4).*pa(3:4); % width and height of panel in pixels
    dcut = dcut*sqrt(p.Results.dot_size)/as(1)*(range(unique(x))+1);
    if ~ishold
        cla
    else
        keep_hold = true;
    end
end

% sort/round y for different plot styles
yorig=y;
switch lower(p.Results.sort_style)
    case 'up'
        [y,sid]=sort(y);
    case 'fan'
        [~,sid]=sort(abs(y-mean(y)));
        sid=[sid(1:2:end); sid(2:2:end)];
        y=y(sid);
    case 'down'
        [y,sid]=sort(y,'descend');
    case 'rand'
        sid=randperm(length(y));
        y=y(sid);
    case 'square'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/asp_rat));
        [~,e,b]=histcounts(y,edges);
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
    case 'hex'
        nxloc=.9/dcut;
%         [~,e,b]=histcounts(y,ceil((range(x)+1)*chanwid*nxloc/2/sqrt(1-.5.^2)/asp_rat));
        edges = linspace(min(yl),max(yl),ceil((range(x)+1)*chanwid*nxloc/sqrt(1-.5.^2)/asp_rat));
        [n,e,b]=histcounts(y,edges);
        oddmaj=0;
        if sum(mod(n(1:2:end),2)==1)>sum(mod(n(2:2:end),2)==1),
            oddmaj=1;
        end
        y=e(b)'+mean(diff(e))/2;
        [y,sid]=sort(y);
        b=b(sid);
    otherwise
        sid=1:length(y);
end
x=x(sid);
yorig=yorig(sid);
[ux,~,ic] = unique(x);
% rmult=(range(ux)+1)*2;
rmult=5;

% for each group...
for i=1:length(ux)
    fid = find(ic==i);   
    
    % set of possible x locations
    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);

    % rescale y to that things are square visually
    zy=(y(fid)-min(yl))/(max(yl)-min(yl))/asp_rat*(range(ux)+1)*chanwid;
    
    % precalculate y distances so that we only worry about nearby points
    D0=squareform(pdist(zy))<dcut*2;    
    
    if length(fid)>1
        % for each data point in the group sequentially...
        for j=1:length(fid)
            if strcmp(lower(p.Results.sort_style),'hex')
                xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i);
                if mod(b(fid(j)),2)==oddmaj
                    xi = linspace(-chanwid/2*rmult,chanwid/2*rmult,nxloc*rmult+(mod(nxloc*rmult,2)==0))'+ux(i)+mean(diff(xi))/2;
                end
            end
            zid = D0(j,1:j-1);
            e = (xi-ux(i)).^2; % cost function
            if ~strcmp(lower(p.Results.sort_style),'hex') && ~strcmp(lower(p.Results.sort_style),'square')
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D<=dcut)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            else
                if sum(zid)>0
                    D = pdist2([xi ones(length(xi),1)*zy(j)], [x(fid(zid)) zy(zid)]);
                    D(D==0)=Inf;
                    D(D>dcut & isfinite(D))=0;
                    e = e + sum(D,2) + randn(1)*10e-6; % noise to tie-break
                end
            end

            if strcmp(lower(p.Results.sort_style),'one')
                e(xi<ux(i))=Inf;
            end
            [~,mini] = min(e);
            if mini==1 && rand(1)>.5, mini=length(xi); end
            x(fid(j)) = xi(mini);
        end
    end
%     x(fid)=x(fid)-median(x(fid))+ux(i); % center x locations by median
end

if strcmp(lower(p.Results.sort_style),'randn')
    x=ux(ic)+randn(size(ic))/4;
end

% corral any points outside of the channel
out_of_range = abs(x-ux(ic))>chanwid/2;
switch lower(p.Results.corral_style)
    case 'gutter'
        id = (x-ux(ic))>chanwid/2;
        x(id)=chanwid/2+ux(ic(id));
        id = (x-ux(ic))<-chanwid/2;
        x(id)=-chanwid/2+ux(ic(id));
    case 'omit'
        x(out_of_range)=NaN;
    case 'random'
        x(out_of_range)=ux(ic(out_of_range))+rand(sum(out_of_range),1)*chanwid-chanwid/2;
end

% plot groups and add overlay
if isfinite(p.Results.dot_size)
    if isnumeric(p.Results.colormap)
        cmap=p.Results.colormap;
    else
        cmap = feval(p.Results.colormap,length(ux));
    end
    for i=1:length(ux)
        if isempty(p.Results.MarkerFaceColor')
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha',p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',cmap(i,:))
        else
            scatter(x(ic==i),y(ic==i),p.Results.dot_size*36,'filled','MarkerFaceAlpha',p.Results.MarkerFaceAlpha,'MarkerEdgeColor',p.Results.MarkerEdgeColor,'MarkerFaceColor',p.Results.MarkerFaceColor)
        end
        hold on
        iqr = prctile(yorig(ic==i),[25 75]);
        switch lower(p.Results.overlay_style)
            case 'box'
                rectangle('Position',[ux(i)-rwid iqr(1) 2*rwid iqr(2)-iqr(1)],'EdgeColor','k','LineWidth',2)
                line([ux(i)-rwid ux(i)+rwid],[1 1]*median(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'sd'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i)),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'ci'
                line([1 1]*ux(i),mean(yorig(ic==i))+[-1 1]*std(yorig(ic==i))/sqrt(sum(ic==i))*tinv(0.975,sum(ic==i)-1),'Color',cmap(i,:),'LineWidth',2)
                line([ux(i)-2*rwid ux(i)+2*rwid],[1 1]*mean(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
            case 'median'
                line([ux(i)-0.25 ux(i)+0.25],[1 1]*median(yorig(ic==i)),'LineWidth',3,'Color',cmap(i,:))
        end
        
    end
    hold off
    if keep_hold
        hold on
    end
    xlim(xl)
    ylim(yl)
end
end