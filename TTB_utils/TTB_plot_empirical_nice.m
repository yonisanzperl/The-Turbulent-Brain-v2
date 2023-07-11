function TTB_plot_empirical_nice(output,Cfg,metadata)



make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.12 0.015], [0.12 0.01], [0.12 0.01]);
if ~make_it_tight,  clear subplot;  end

LAMBDA = Cfg.iLambda:Cfg.stepsLam:Cfg.fLambda;
LAMBDA = flip(LAMBDA);
[dd ilambda]=min(abs(LAMBDA-Cfg.PlotLambda));

% Turbulence

groups =metadata.group;
stattest = metadata.stattest;

NoComp = Cfg.nBrainStates*(Cfg.nBrainStates-1)*0.5;
pval = zeros(length(LAMBDA),NoComp);
C = cell(1,Cfg.nBrainStates);
f1=figure('Name','Turbulence across scales');
for i=1:length(LAMBDA)
    cont=1;
    for ii=1:Cfg.nBrainStates
        eval(['turbu_' num2str(ii) '=output(ii).Turbulence_sub(i,:);']);
        C{1,ii}= output(ii).Turbulence_sub(i,:);

    end
    subplot(3,ceil(length(LAMBDA)/3),i)
    fprintf('\n \n Stats for Turbu at lambda: %f \n',LAMBDA(i))    ;
    p = swarm(C, groups, tlt=sprintf('lambda %.2f',LAMBDA(i)), overlay_style='boxplot', printPvals=true, stat_test=stattest,name=sprintf('Model-free_turbuStats_lambda %.2f.txt',LAMBDA(i)));
    title(sprintf('lambda %.2f',LAMBDA(i)));
    movefile(sprintf('Model-free_turbuStats_lambda %.2f.txt',LAMBDA(i)), metadata.outdir)


end

f1.Position = [100 100 740 600];
% Information Transfer



K=3; % arbitrary number to create a meaningfull transfer
NoComp = Cfg.nBrainStates*(Cfg.nBrainStates-1)*0.5;
pvalTr = zeros(length(LAMBDA),NoComp);
C = cell(1,Cfg.nBrainStates);
f2=figure('Name','Transfer across scales');
for i=1:length(LAMBDA)
    cont=1;
    for ii=1:Cfg.nBrainStates

        eval(['transfer' num2str(ii) '=output(ii).Turbulence_sub(i,:);']);
        C{1,ii}= K-output(ii).Transfer_sub(i,:);
    end
    subplot(3,ceil(length(LAMBDA)/3),i)
    fprintf('\n \n Stats for Transfer at lambda: %f \n',LAMBDA(i))    ;
    p = swarm(C, groups, tlt=sprintf('lambda %.2f',LAMBDA(i)), overlay_style='boxplot', printPvals=true, stat_test=stattest,name=sprintf('Model-free_transferStats_lambda %.2f.txt',LAMBDA(i)));
    title(sprintf('lambda %.2f',LAMBDA(i)));
        movefile(sprintf('Model-free_transferStats_lambda %.2f.txt',LAMBDA(i)), metadata.outdir)

end

f2.Position = [100 100 740 600];


%Information Flow
col ={'r','b','k','g'};
col_mar ={'r-o','b-o','k-o','g-o'};
%figure('Name','Info Flow');

for ii=1:Cfg.nBrainStates
    clear TL H

    TL = output(ii).TransferLambda_sub;
    InfoCas{1,ii}= nanmean(TL(2:length(LAMBDA),:));
  %  subplot(2,ceil(Cfg.nBrainStates/2),ii)
  %  shadedErrorBar(LAMBDA(2:length(LAMBDA)), mean(TL(2:length(LAMBDA),:),2)',std(TL(2:length(LAMBDA),:),[],2)','lineprops',{col_mar{ii},'markerfacecolor',col{ii}});
   % plot(LAMBDA(2:length(LAMBDA)), mean(TL(2:length(LAMBDA),:),2)')
   % ylabel('Info Flow');xlabel('Lambda Pairs')
end

 

% Information Cascade
figure('Name','Info Cascade');
p = swarm(InfoCas, groups, tlt='', overlay_style='boxplot', printPvals=true, stat_test=stattest,name='ModelfreeInfCasStats.txt');
movefile('ModelfreeInfCasStats.txt', metadata.outdir)


% Turbulence by RSN


make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.1 0.1], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end


NoComp = Cfg.nBrainStates*(Cfg.nBrainStates-1)*0.5;
pval = zeros(length(LAMBDA),NoComp);
C = cell(1,8);

groupsrsn ={'Vis','Som','Att','Sal','Lim','Con','DMN','TP'}

i=ilambda;
figure('Name',sprintf('Turbulence by RSN across at lambda %.2f',LAMBDA(i))');

cont=1;
for ii=1:Cfg.nBrainStates
    for rsn=1:8
        eval(['turbu_rsn_' num2str(ii) '=output(ii).TurbulenceRSN_sub(i,:,:);']);
        C{1,rsn}= squeeze(output(ii).TurbulenceRSN_sub(i,:,rsn));
        %groupsrsn{rsn} = sprintf('RNS %d',rsn);
    end
    
    subplot(Cfg.nBrainStates,2,cont)
    cont = cont+1;
    fprintf('\n \n Stats for Turbu by RSN at lambda: %f \n',LAMBDA(i))    ;
    [p, stat_turRSN] = swarm(C, groupsrsn, tlt=sprintf('lambda %.2f',LAMBDA(i)), overlay_style='boxplot', printPvals=true, stat_test=stattest,name='ModelfreeTurbuStatsRSN.txt',display_on_plot=false);
    title(groups{ii});
    subplot(Cfg.nBrainStates,2,cont)
    cont =cont+1;
    
    stat_turRSN2=zeros(8,8);
    stat_turRSN2(1:7,1:8)=stat_turRSN;
    imagesc(stat_turRSN2'<0.05 & 0<stat_turRSN2');
    set(gca,'fontname','times')  % Set it to times
    set(gca, 'FontSize', 12); hold on;
    xticks([1:1:length(groupsrsn)]);
    xticklabels(groupsrsn);
    yticks([1:1:length(groupsrsn)]);
    yticklabels(groupsrsn);
    axis square
    title(['p_{val}<0.05 at ' groups{ii}]);
    
    %title(sprintf('RSN stats Condition %d',ii)');
    caxis([0, 1])
    myColorMap = jet(256);
    myColorMap(1,:) = 1;
    colormap(myColorMap);
end

figure('Name',sprintf('Turbulence by RSN inter condition across at lambda %.2f',LAMBDA(i))');

for rsn=1:8
    for ii=1:Cfg.nBrainStates
        C2{1,ii}=squeeze(output(ii).TurbulenceRSN_sub(i,:,rsn));
    end
    for kk=1:Cfg.nBrainStates
        for jj=kk+1:Cfg.nBrainStates
            stat_inter(kk,jj)=ranksum(C2{1,kk},C2{1,jj});
            
        end
    end
    subplot(4,2,rsn)
    stat_inter2 =zeros(Cfg.nBrainStates,Cfg.nBrainStates);
    stat_inter2(1:Cfg.nBrainStates-1,1:Cfg.nBrainStates)=stat_inter;
    imagesc(stat_inter2'<0.05 & 0<stat_inter2')
    set(gca,'fontname','times')  % Set it to times
    set(gca, 'FontSize', 12); hold on;
    xticks([1:1:length(groups)]);
    xticklabels(groups);
    yticks([1:1:length(groups)]);
    yticklabels(groups);
    axis square
    title(sprintf('p_{val}<0.05 at %s',groupsrsn{rsn}));
    
    %title(sprintf('RSN stats Condition %d',ii)');
    caxis([0, 1])
    myColorMap = jet(256);
    myColorMap(1,:) = 1;
    colormap(myColorMap);
    
end



movefile('ModelfreeTurbuStatsRSN.txt', metadata.outdir)





% 
% % ---------Brain Render----------
% 

  % lambda to plot renders

for ii=1:Cfg.nBrainStates

    mean_nodeTurb_sub(ii,:,:)=squeeze(mean(output(ii).node_Turbulence_sub,2));
    std_nodeTurb_sub(ii,:,:)=squeeze(std(output(ii).node_Turbulence_sub,[],2));
end
tmp=squeeze(mean_nodeTurb_sub(:,ilambda,:));
lower= min(min(tmp));
upper = max(max(tmp));
for ii=1:Cfg.nBrainStates
    
    rendersurface_schaefer1000(squeeze(mean_nodeTurb_sub(ii,ilambda,:)),lower,upper,0,'Greens9',1,['Mean Node level Turb. for' groups{ii} sprintf(' at lambda %.2f',LAMBDA(ilambda))])
    %rendersurface_schaefer1000(squeeze(std_nodeTurb_sub(ii,ilambda,:)),min(squeeze(std_nodeTurb_sub(ii,ilambda,:))),max(squeeze(std_nodeTurb_sub(ii,ilambda,:))),0,'Greens9',1, ['Std Node level Turb.' sprintf('lambda %.2f',LAMBDA(ilambda))])
end
    
end


