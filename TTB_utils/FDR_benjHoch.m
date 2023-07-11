% Andrea Insabato 15/01/15
% andrea.insabato@upf.edu
%
% False Discovery Rate - Benjamini-Hochberg method
%
% rejectedH0s = FDR_benjHoch(pvals, q, corr, bootstrapSTD, fig_pvals, fig_var_q)
%
% Controls the False Discovery Rate of a set of tests with p-values pvals to a level q with the
% Benjamini-Hochberg procedure.
% pvals is a vector of p-values for each null hypothesis tested.
% q is the desired threshold.
% corr is an optional string ('positive' or 'arbitrary') indicating if there is correlation between tests (see below).
% bootstrapSTD is a boolean indicating whether to use bootstrap (slower) or
% parametric method (with a binomial distribution) to estimate std of
% number of discoveries.
% fig_pvals is an optional boolean determining whether to plot pvalues and criterias.
% fig_var_q is an optional boolean determining whether to plot the number
% of total discoveries and the number of accepted false discoveries when
% varying the q-value.
% It returns a vector rejectedH0s with the indices of the rejected null hypothesis
% corresponding to the tests in pvals. In case no test is significant an
% empty matrix is returned.
%
%
% Correlation between tests:
% Benjamini and Yakutieli (2001) proved that under positive correlation BH procedure
% controls FDR to a given level q.
% Under arbitrary correlations a modification of the procedure still controls FDR to
% the level q but it is more conservative.
% A resampling method can be used instead.
%
% TODO: implement the resampling method, implement the pFDR (Storey 2003)

function [rejectedH0s, rejectedH0ss, qs] = FDR_benjHoch(pvals, q, correlation, bootstrapSTD, fig_pvals, fig_var_q)

% if pvals is row vector make it column (bootstrp samples from columns)
if isrow(pvals)
    pvals = pvals';
end

if nargin < 3
	c = 1; 
else
	switch correlation
	case 'positive'
		c = 1;
	case 'arbitrary'
		c = cumsum(1./[1:length(pvals)]);
	otherwise
		error('FDR_benjHoch: corr can be either ''positive'' or ''arbitrary''')
	end
end

if nargin < 4
	bootstrapSTD = false;
end

if nargin < 5
	fig_pvals = false;
end

if nargin < 6
	fig_var_q = false;
end

[ordered_pval, pvalIdx] = sort(pvals,'ascend');
criterias = q * [1:length(pvals)]' ./ (length(pvals)*c);

num_rejectedH0s = find(ordered_pval <= criterias,1,'last');  % reject ordered H0s from 1 to this index
rejectedH0s = pvalIdx(1:num_rejectedH0s);

if fig_pvals
	figure
	plot(ordered_pval, '-k', 'linewidth', 3)
	hold on
	plot(criterias, '--r', 'linewidth', 2)
	set(gca, 'xticklabel', pvalIdx, 'fontsize', 16)
	xlabel('indexes of ordered p-values')
	ylabel('ordered p-values')
	legend({'p-values'; 'criterias'})
	set(gca, 'fontsize', 16)
end



count=0;
qs = [0.01, 0.05:0.05:1];
num_rejectedH0ss = zeros(1,length(qs));
for q_alpha = qs
    count=count+1;
    criteriass = q_alpha * [1:length(pvals)]' / (length(pvals)*c);
    reject = find(ordered_pval <= criteriass,1,'last');
    if ~isempty(reject)
        num_rejectedH0ss(count) = reject;
    end
    rejectedH0ss{count} = pvalIdx(1:reject);

    if bootstrapSTD
        % calc standard dev with bootstrap
        rej_var(count) = std(...
            bootstrp(10000,...
            @(bs_pvals) find(bs_pvals <= criteriass,1,'last'),...
            ordered_pval));
    else
        % calc standard dev parametrically (binomial)
        rej_var(count) = sqrt( length(pvals)*num_rejectedH0ss(count)/length(pvals)*(1-num_rejectedH0ss(count)/length(pvals)));
    end
end

if fig_var_q
    figure
    bar(qs, num_rejectedH0ss, 'facecolor',[.9 .9 .9],'edgecolor','k', 'linewidth',3), hold on
    errorbar(qs, num_rejectedH0ss, rej_var,'.','linewidth',3,'color',[.8 .2 0],'markersize',2)
    h1=plot(qs,num_rejectedH0ss.*qs,'r--','linewidth',3);
    h2=plot(qs,num_rejectedH0ss - num_rejectedH0ss.*qs,'k','linewidth',3);
    set(gca, 'fontsize',16)
    xlabel('threshold of FDR q', 'fontsize',16)
    legend([h1, h2], 'false discoveries', 'true discoveries')
end