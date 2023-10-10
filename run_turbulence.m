clear all
% close all
%clearvars -except Answer0
close all
addpath(genpath([pwd '/TTB_utils']));

%% ESTO ES LO ÚNICO QUE SE TOCA (sino armar una copia)

warning('off')

%figure('Position',140 50)
% imshow(imread('Fig1_ttb.png'));
% pause(5)
% close

% if exist('Answer0','var')
%     [Answer,Answer0,Cancelled]=IN_turbulence2(Answer0);
% else
%     [Answer,Answer0,Cancelled]=IN_turbulence2;
% end

a0 = questdlg('What approach do you want to compute?', ...
    'The Turbulent Brain', ...
    'Model-free','Model-based','No thank you');

% Handle response
switch a0
    case 'Model-free'
        
        
        [Answer,Answer0,Cancelled]=IN_turbulence_mf;
        metadata.saveFile = Answer.saveFile;
        metadata.currentDir=pwd;
        
        Cfg.ID = 'TTB computation';
        
        
        metadata.TSCh=Answer.TSCh{1};%Nombre de la variable de la tseries para metadata
        metadata.SCCh=Answer.SCCh{1};%Nombre de la variable de la sc para metadata
        metadata.CoGCh=Answer.GroupCoG{1};
        metadata.RSNCh=Answer.RSNCh{1};

        % load time series, if it is not provided stop the process
        if strcmp(metadata.TSCh, 'timeseries are not included');
            warndlg('Do you need to provide timeseries')

        else
            
            TS=load(Answer.TSDir,Answer.TSCh{1}); %Cell en formato (1,Nsubs) donde cada celda tiene una matriz (nodes x T).
            TS=TS.(Answer.TSCh{1}); % Este es un truco para tomar deshacer el anidado que se genera de hacer el load

         end
        
        Cfg.TR=Answer.TRsec; %Período muestral de las tseries
        
     
        
        
        %load de Cog from file or from data from TTB       
         if strcmp(metadata.CoGCh, 'Centre of mass is not included');
             CoG = load([metadata.currentDir '/TTB_data/schaefercog.mat']);
             CoG = CoG.SchaeferCOG;
         else     
             CoG=load(Answer.CoGDir,Answer.GroupCoG{1}); %Matriz (nNodes x nNodes)
             CoG=CoG.(Answer.GroupCoG{1}); % Este es un truco para deshacer el anidado que se genera de hacer el load
         end
         
         %load de RSN belonging from file or from data from TTB
         if strcmp(metadata.RSNCh, 'RSN belonging is not included');
             RSN = load([metadata.currentDir '/TTB_data/RSN_yeo17_schaefer1000.mat']);
             RSN = RSN.labels(:,2);
         else
             RSN=load(Answer.RSNDir,Answer.RSNCh{1}); %Matriz (nNodes x nNodes)
             RSN = RSN.labels(:,2);
         end
  
         %load de SC from file or from data from TTB
         if strcmp(metadata.SCCh, 'Structural Connectivity is not included');
             SC = load([metadata.currentDir '/TTB_data/sc_schaefer_MK.mat']);
             SC =SC.sc_schaefer;
         else
            SC=load(Answer.SCDir,Answer.SCCh{1}); %Matriz (nNodes x nNodes)
            SC=SC.(Answer.SCCh{1}); % Este es un truco para deshacer el anidado que se genera de hacer el load
  
         end
         if iscell(SC)
             nNodes=length(SC{1});
         else
             nNodes=length(SC);
         end
         
         Cfg.nNodes=nNodes;
        
        
        Cfg.randseed='shuffle'; %Tipo de semilla para los números random
        
        
        Cfg.Tmax=Answer.Tmax; %Esto es el Tmax que va a cortar y que va a tener cada sujeto simulado (la tseries simulada va a ser Cfg.Tmax*nSubs)
        %si es 0 significa que va a tomar el maximo de cada sujeto
        Cfg.filt.bpass =Answer.EnableFilterMode;  %1: filtra la tseries entrante (DFLT=1)
        Cfg.filt.lb= Answer.BPlb; %frecuencia menor del pasabanda (DFLT=0.04)
        Cfg.filt.ub= Answer.BPub; %frecuencia mayor del pasabanda (DFLT=0.07)
        
        
        Cfg.NR = 400; %hardcoded variabled from Turbu paper
        Cfg.NRini = 20;%hardcoded variabled from Turbu paper
        Cfg.NRfin = 80;%hardcoded variabled from Turbu paper
        Cfg.iLambda = Answer.iLambda;
        Cfg.fLambda = Answer.fLambda;
        Cfg.stepsLam = Answer.stepsLam;
        Cfg.PlotEmpYes = Answer.PlotEmpYes;
        Cfg.nBrainStates = str2num(Answer.nBS{1});
        metadata.stattest = Answer.stattest{1}; 
        Cfg.PlotLambda=Answer.PlotLambda;
        
        for kk=1:Cfg.nBrainStates
            group(kk)=inputdlg(sprintf('Enter the name of Condition %d',kk));
        end
                
        
        % compute long range connection from SC
        [Clong lambda_emp numexcSC valexcSC]=TTB_SC_analysis_longrange_EDR_DTI(Cfg,CoG,SC);
        
        % compute turbulence metrics
        for kk=1:Cfg.nBrainStates
            fprintf('Brain state: %d/%d\n',kk,Cfg.nBrainStates)
            [output(kk)]=TTB_turbulence_node_empirical_metaesta_RSN(Cfg, CoG,SC,TS(:,kk),RSN,Clong, lambda_emp, numexcSC, valexcSC);
        end
        
        metadata.group = group;
        metadata.RSNCh=Answer.RSNCh{1}; %
        metadata.currentFolder = pwd;

        out_Dir=[metadata.currentFolder '/TTB_results/' metadata.saveFile '_' date]; % the direction of the bars depends on Win o Lin
        mkdir (out_Dir);
        metadata.outdir=out_Dir;
       
        clear out_Dir
        
        
                % plot emp turbulence
        if Cfg.nBrainStates>1
<<<<<<< HEAD
        if Cfg.PlotEmpYes
            TTB_plot_empirical(output,Cfg,metadata)
=======
            if Cfg.PlotEmpYes
                TTB_plot_empirical(output,Cfg,metadata)
>>>>>>> a6e61b4e41951838721886551af09332f84e78ab
            
            end
        end
        end
        cd(metadata.outdir)
        save('Results_modelfree.mat','output','SC','CoG','RSN','metadata','Cfg');
        cd(metadata.currentDir);

    case 'Model-based'
 
        
        
        a1 = questdlg('Did you already compute model-free appoach?', ...
            'The Turbulent Brain', ...
            'Yes','No','No thank you');
        
        % Handle response
        switch a1
            case 'Yes'
               
                                
                 [Answer,Answer0,Cancelled]=IN_turbulence_mb1;
               
                metadata.TSCh=Answer.TSCh{1};%Nombre de la variable de la tseries para metadata
                metadata.SCCh=Answer.allIt{1};%Nombre de la variable de la sc para metadata
                TS=load(Answer.TSDir,Answer.TSCh{1}); %Cell en formato (1,Nsubs) donde cada celda tiene una matriz (nodes x T).
                TS=TS.(Answer.TSCh{1}); % Este es un truco para tomar deshacer el anidado que se genera de hacer el load
                load(Answer.MFoutdir)
                metadata.saveFile_old = metadata.saveFile;
                metadata.saveFile = Answer.saveFile;
                metadata.currentDir=pwd;
                Cfg.ID = 'TTB model based and free computation';
                % Hopf model section
                Cfg.Tsim=Answer.Tsim;
                Cfg.Glower = Answer.GLb;
                Cfg.Gupper = Answer.GUb;
                Cfg.Gstep = Answer.Gstep;
                Cfg.PlotYes = Answer.PlotYes;
                Cfg.PertuYes=Answer.PertuYes;
                Cfg.NSIM=Answer.NSIM;
                Cfg.lamICS=Answer.lamICS;
                Cfg.nTrials=Answer.nTrials;
                
                
                % to compute the Hopf model end extract Turbu metrics from the modeled
                % signals
                G_range=Cfg.Glower:Cfg.Gstep:Cfg.Gupper;
                [Clong lambda_emp numexcSC valexcSC]=TTB_SC_analysis_longrange_EDR_DTI(Cfg,CoG,SC);
                if Cfg.PlotYes
                    figure('Name','Model fitting');
                end
                for kk=1:Cfg.nBrainStates
                    fprintf('Modeled Brain state: %d/%d\n',kk,Cfg.nBrainStates)
                    [corrfcn,f_diff,fce]=TTB_Empirical_corrfcn_freq(Cfg,SC,CoG,TS(:,kk));
                    [err_hete(kk,:,:) model_output(kk,:)]=TTB_hopf_DTI_Grange_v1(Cfg,TS(:,kk),CoG,SC,corrfcn,f_diff,RSN,Clong, lambda_emp, numexcSC, valexcSC);
                    if Cfg.PlotYes
                       subplot(1,Cfg.nBrainStates,kk)
                        plot(G_range,mean(squeeze(err_hete(kk,:,:)),1))
                    end
                    [min_val(kk) indx_min(kk)]=min(mean(squeeze(err_hete(kk,:,:)),1));
                    if Cfg.PertuYes
                        fprintf('Perturbed Brain state: %d/%d\n',kk,Cfg.nBrainStates)
                        [IC(kk,:) S(kk,:) ]= TTB_infocapacity_hopf(Cfg.nTrials,kk,G_range(indx_min(kk)),Cfg,CoG,SC,f_diff,Clong, lambda_emp, numexcSC, valexcSC);
                    end
                end
                
                

                
                out_Dir=[metadata.currentFolder '/TTB_results/' metadata.saveFile '_model_' date]; % the direction of the bars depends on Win o Lin
                mkdir (out_Dir);
                metadata.outdir=out_Dir;
                
                clear out_Dir
                        cd(metadata.currentDir);

                if Cfg.PertuYes
                    TTB_plot_modelIC(model_output,Cfg, IC, S, metadata)
                    cd(metadata.outdir)
                    save('TTBoutput_modelbased.mat','model_output','IC','S','metadata','Cfg','G_range','err_hete');

                end
                cd(metadata.outdir)

                save('TTBoutput_modelbased.mat','model_output','metadata','Cfg','G_range','err_hete');
                
                        cd(metadata.currentDir);

                
                
            case 'No'
                
                [Answer,Answer0,Cancelled]=IN_turbulence_mb2;
                metadata.saveFile = Answer.saveFile;
                metadata.currentDir=pwd;
                
                Cfg.ID = 'TTB computation';
                
                
                metadata.TSCh=Answer.TSCh{1};%Nombre de la variable de la tseries para metadata
                metadata.SCCh=Answer.SCCh{1};%Nombre de la variable de la sc para metadata
                metadata.RSNCh=Answer.RSNCh{1};%Nombre de la variable de la tseries para metadata

                metadata.CoGCh=Answer.GroupCoG{1};

                
                % load time series, if it is not provided stop the process
                if strcmp(metadata.TSCh, 'timeseries are not included');
                    warndlg('Do you need to provide timeseries')
                    
                else
                    
                    TS=load(Answer.TSDir,Answer.TSCh{1}); %Cell en formato (1,Nsubs) donde cada celda tiene una matriz (nodes x T).
                    TS=TS.(Answer.TSCh{1}); % Este es un truco para tomar deshacer el anidado que se genera de hacer el load
                    
                end
                
                Cfg.TR=Answer.TRsec; %Período muestral de las tseries
                
                
                
                
                %load de Cog from file or from data from TTB
                if strcmp(metadata.CoGCh, 'Centre of mass is not included');
                    CoG = load([metadata.currentDir '/TTB_data/schaefercog.mat']);
                    CoG = CoG.SchaeferCOG;
                else
                    CoG=load(Answer.CoGDir,Answer.GroupCoG{1}); %Matriz (nNodes x nNodes)
                    CoG=CoG.(Answer.GroupCoG{1}); % Este es un truco para deshacer el anidado que se genera de hacer el load
                end
                
                %load de SC from file or from data from TTB
                if strcmp(metadata.SCCh, 'Structural Connectivity is not included');
                    SC = load([metadata.currentDir '/TTB_data/sc_schaefer_MK.mat']);
                    SC =SC.sc_schaefer;
                else
                    SC=load(Answer.SCDir,Answer.SCCh{1}); %Matriz (nNodes x nNodes)
                    SC=SC.(Answer.SCCh{1}); % Este es un truco para deshacer el anidado que se genera de hacer el load
                    
                end
                if iscell(SC)
                    nNodes=length(SC{1});
                else
                    nNodes=length(SC);
                end
                
                Cfg.nNodes=nNodes;
                
                
                
                %load de RSN belonging from file or from data from TTB
                if strcmp(metadata.RSNCh, 'RSN belonging is not included');
                    RSN = load([metadata.currentDir '/TTB_data/RSN_yeo17_schaefer1000.mat']);
                    RSN = RSN.labels(:,2);
                else
                    RSN=load(Answer.RSNDir,Answer.RSNCh{1}); %Matriz (nNodes x nNodes)
                    RSN = RSN.labels(:,2);
                end
                
                
                Cfg.randseed='shuffle'; %Tipo de semilla para los números random
                
                
                %si es 0 significa que va a tomar el maximo de cada sujeto
                Cfg.filt.bpass =Answer.EnableFilterMode;  %1: filtra la tseries entrante (DFLT=1)
                Cfg.filt.lb= Answer.BPlb; %frecuencia menor del pasabanda (DFLT=0.04)
                Cfg.filt.ub= Answer.BPub; %frecuencia mayor del pasabanda (DFLT=0.07)
                
                
                Cfg.NR = 400; %hardcoded variabled from Turbu paper
                Cfg.NRini = 20;%hardcoded variabled from Turbu paper
                Cfg.NRfin = 80;%hardcoded variabled from Turbu paper
                Cfg.iLambda = Answer.iLambda;
                Cfg.fLambda = Answer.fLambda;
                Cfg.stepsLam = Answer.stepsLam;

                Cfg.nBrainStates = str2num(Answer.nBS{1});

                
                for kk=1:Cfg.nBrainStates
                    group(kk)=inputdlg(sprintf('Enter the name of Condition %d',kk));
                end
                
                
                % Hopf model section
                Cfg.Tmax=Answer.Tsim;
                Cfg.Glower = Answer.GLb;
                Cfg.Gupper = Answer.GUb;
                Cfg.Gstep = Answer.Gstep;
                Cfg.PlotYes = Answer.PlotYes;
                Cfg.PertuYes=Answer.PertuYes;
                Cfg.NSIM=Answer.NSIM;
                Cfg.lamICS=Answer.lamICS;
                Cfg.nTrials=Answer.nTrials;
                metadata.group = group;
                
                % to compute the Hopf model end extract Turbu metrics from the modeled
                % signals
                G_range=Cfg.Glower:Cfg.Gstep:Cfg.Gupper;
                
                [Clong lambda_emp numexcSC valexcSC]=TTB_SC_analysis_longrange_EDR_DTI(Cfg,CoG,SC);
                
                if Cfg.PlotYes
                    figure('Name','Model fitting');
                end
                for kk=1:Cfg.nBrainStates
                    fprintf('Modeled Brain state: %d/%d\n',kk,Cfg.nBrainStates)
                    [corrfcn,f_diff,fce]=TTB_Empirical_corrfcn_freq(Cfg,SC,CoG,TS(:,kk));
                    [err_hete(kk,:,:) model_output(kk,:)]=TTB_hopf_DTI_Grange_v1(Cfg,TS(:,kk),CoG,SC,corrfcn,f_diff,RSN,Clong, lambda_emp, numexcSC, valexcSC);
                    if Cfg.PlotYes
                        subplot(1,Cfg.nBrainStates,kk)
                        plot(G_range,mean(squeeze(err_hete(kk,:,:)),1))
                    end
                    [min_val(kk) indx_min(kk)]=min(mean(squeeze(err_hete(kk,:,:)),1));
                    if Cfg.PertuYes
                        fprintf('Perturbed Brain state: %d/%d\n',kk,Cfg.nBrainStates)
                        [IC(kk,:) S(kk,:) ]= TTB_infocapacity_hopf(Cfg.nTrials,kk,G_range(indx_min(kk)),Cfg,CoG,SC,f_diff,Clong, lambda_emp, numexcSC, valexcSC);
                    end
                end
                
                
                metadata.currentFolder = pwd;
                
                
                
                out_Dir=[metadata.currentFolder '/TTB_results/' metadata.saveFile '_' date]; % the direction of the bars depends on Win o Lin
                mkdir (out_Dir);
                metadata.outdir=out_Dir;
                
                clear out_Dir
                cd(metadata.currentDir);
                
                if Cfg.PertuYes
                    TTB_plot_modelIC(model_output,Cfg, IC, S, metadata)
                    cd(metadata.outdir)          
                    save('TTBoutput_modelbased.mat','model_output','IC','S','metadata','Cfg','G_range','err_hete');

                end
                cd(metadata.outdir)
                
                save('TTBoutput_modelbased.mat','model_output','IC','S','metadata','Cfg','G_range','err_hete');
                
                cd(metadata.currentDir);
                
                
        end
        
end



cd(metadata.currentDir);
%clear all;


