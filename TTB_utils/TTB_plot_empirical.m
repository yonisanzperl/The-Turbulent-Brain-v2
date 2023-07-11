function TTB_plot_empirical(output,Cfg)


LAMBDA = Cfg.iLambda:Cfg.stepsLam:Cfg.fLambda;

% Turbulence

NoComp = Cfg.nBrainStates*(Cfg.nBrainStates-1)*0.5;
pval = zeros(length(LAMBDA),NoComp);
C = cell(1,Cfg.nBrainStates);
figure(1);
for i=1:length(LAMBDA)
    cont=1;
    for ii=1:Cfg.nBrainStates
        for kk=ii+1:Cfg.nBrainStates            
            pval(i,cont) = ranksum(output(ii).Turbulence_sub(i,:),output(kk).Turbulence_sub(i,:));
            comparison{cont}=[ii,kk];
            cont = cont+1;
        end
        eval(['turbu_' num2str(ii) '=output(ii).Turbulence_sub(i,:);']);
        C{1,ii}= output(ii).Turbulence_sub(i,:);
        groups{ii} = sprintf('Cond%d',ii); 
    end
    subplot(4,ceil(length(LAMBDA)/4),i)
    title('Turbulence across scales')
    maxNumEl = max(cellfun(@numel,C));
    Cpad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, C);
    Cmat = cell2mat(Cpad);
    boxplot(Cmat,'Labels',groups)
    H=sigstar(comparison,pval(i,:));

    ylabel(sprintf('lambda %.2f',LAMBDA(i)));
end


% Information Transfer

K=3;
NoComp = Cfg.nBrainStates*(Cfg.nBrainStates-1)*0.5;
pvalTr = zeros(length(LAMBDA),NoComp);
C = cell(1,Cfg.nBrainStates);
figure(2);
for i=1:length(LAMBDA)
    cont=1;
    for ii=1:Cfg.nBrainStates
        for kk=ii+1:Cfg.nBrainStates            
            pvalTr(i,cont) = ranksum(output(ii).Transfer_sub(i,:),output(kk).Transfer_sub(i,:));
            comparison{cont}=[ii,kk];
            cont = cont+1;
        end
        eval(['transfer' num2str(ii) '=output(ii).Turbulence_sub(i,:);']);
        C{1,ii}= K-output(ii).Transfer_sub(i,:);
        groups{ii} = sprintf('Cond%d',ii); 
    end
    subplot(4,ceil(length(LAMBDA)/4),i)
    title('Transfer across scales')
    maxNumEl = max(cellfun(@numel,C));
    Cpad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, C);
    Cmat = cell2mat(Cpad);
    boxplot(Cmat,'Labels',groups)
    H=sigstar(comparison,pvalTr(i,:));
    ylabel(sprintf('lambda %.2f',LAMBDA(i)));
end




%Information Flow
col ={'r','b','k','g'};
figure(3);
for ii=1:Cfg.nBrainStates
    
    TL = output(ii).TransferLambda_sub;
    InfoCas{1,ii}= nanmean(TL(2:length(LAMBDA),:));
    subplot(2,ceil(Cfg.nBrainStates/2),ii)
    shadedErrorBar(LAMBDA(2:length(LAMBDA)), TL(2:length(LAMBDA),:)',{@median,@std},'lineprops',{col{ii},'markerfacecolor',col{ii}});
    hold on
    ylabel('Info Flow');xlabel('Lambda Pairs')
end

% 
figure(4)
cont=1;
maxNumEl = max(cellfun(@numel,InfoCas));
Cpad = cellfun(@(x){padarray(x(:),[maxNumEl-numel(x),0],NaN,'post')}, InfoCas);
Cmat = cell2mat(Cpad); 
pvalICas = [];
for ii=1:Cfg.nBrainStates
    for kk=ii+1:Cfg.nBrainStates
        pvalICas(cont) = ranksum(InfoCas{1,ii},InfoCas{1,kk});
        comparisonIC{cont}=[ii,kk];
        cont = cont+1;
    end
end

boxplot(Cmat,'Labels',groups)
H=sigstar(comparisonIC,pvalICas);


% 
% % ---------Brain Render----------
% 
% 
% 
% turbunodecond1=squeeze(nanmean(con1.node_Turbulence_sub(9,:,:),2));
% % including patient that are out and changed
% turbunodecond2=squeeze(nanmean(con2.node_Turbulence_sub(9,setdiff(1:35,[6 11 35]),:),2));
% 
% con3.node_Turbulence_sub = [con3.node_Turbulence_sub(:,setdiff(1:20,5),:) con2.node_Turbulence_sub(:,[6 11],:)]; 
% 
% con34.node_Turbulence_sub = [con3.node_Turbulence_sub con4.node_Turbulence_sub];
% 
% turbunodecond34=squeeze(nanmean(con34.node_Turbulence_sub(9,:,:),2));
% turbunodecond4=squeeze(nanmean(con4.node_Turbulence_sub(9,:,:),2));
% 
% diff1_new =abs(turbunodecond1-turbunodecond2);
% diff2_new =abs(turbunodecond1-turbunodecond34);
% diff3_new =abs(turbunodecond2-turbunodecond34);
% diff4 =abs(turbunodecond2-turbunodecond3);
% diff5 =abs(turbunodecond2-turbunodecond4);
% diff6 =abs(turbunodecond3-turbunodecond4);
% 
% % rendersurface_schaefer1000(diff1,0,0.02,0,'Greens9',1)
% % rendersurface_schaefer1000(diff2,0,0.02,0,'Greens9',1)
% % rendersurface_schaefer1000(diff3,0,0.02,0,'Greens9',1)
% % rendersurface_schaefer1000(diff4,0,0.02,0,'Greens9',1)
% % rendersurface_schaefer1000(diff5,0,0.02,0,'Greens9',1)
% % rendersurface_schaefer1000(diff6,0,0.02,0,'Greens9',1)