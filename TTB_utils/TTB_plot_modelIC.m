function TTB_plot_modelIC(model_output,Cfg, IC, S, metadata)



make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [0.1 0.1], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

LAMBDA = Cfg.iLambda:Cfg.stepsLam:Cfg.fLambda;
LAMBDA = flip(LAMBDA);
groups=metadata.group;

% Info Cap and Susc

NoComp = Cfg.nBrainStates*(Cfg.nBrainStates-1)*0.5;
pval = zeros(1,NoComp);
C = cell(1,Cfg.nBrainStates);
figure('Name','Turbulence across scales');

cont=1;
for ii=1:Cfg.nBrainStates
    C{1,ii}= IC(ii,:);
    Cs{1,ii}= S(ii,:);
end
subplot(1,2,1)
fprintf('\n \n Stats for Info Enc. Cap. at lambda: %f \n',Cfg.lamICS)    ;

p = swarm(C, groups,tlt='Inf.Encoding Capability', overlay_style='boxplot', printPvals=true, fdr=true,stat_test='ranksum',name='ModelICStats.txt');

subplot(1,2,2)
fprintf('\n \n Stats for Susc. at lambda: %f \n',Cfg.lamICS)    ;
p = swarm(Cs, groups, tlt='Susceptibility', overlay_style='boxplot', printPvals=true, fdr=true,stat_test='ranksum',name='ModelSStats.txt');


movefile('ModelICStats.txt', metadata.outdir)
movefile('ModelSStats.txt', metadata.outdir)

