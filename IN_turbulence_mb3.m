
% Si se buguea la GUI
% set(0,'ShowHiddenHandles','on')
% delete(get(0,'Children'))

function [Answer,Answer0,Cancelled]=IN_turbulence_mb3(varargin)
% % 
% if nargin==1
%     Answer0=varargin{1};
% 
%         
%     
%     if ~isfile(Answer0.TSDir)
%     [pathStr, name, ext] = fileparts(Answer0.TSDir)
%     Answer0.TSDir=strcat(pwd,'\Medidas\',name,ext);
%     fprintf('TRATASTE DE CARGAR UN DEFAULT QUE NO EXISTE, usando path local /medidas')
%     %Answer0.TSCh{1}='TS_W';
%     end
%     
%     if ~isfile(Answer0.SCDir)
%     [pathStr, name, ext] = fileparts(Answer0.SCDir)
%     Answer0.SCDir=strcat(pwd,'\Medidas\',name,ext);
%     fprintf('TRATASTE DE CARGAR UN DEFAULT QUE NO EXISTE, usando path local /medidas')
%     %Answer0.SCCh{1}='SC';
%     end
%     
%     if ~isfile(Answer0.GroupDir)
%     [pathStr, name, ext] = fileparts(Answer0.GroupDir)
%     Answer0.GroupDir=strcat(pwd,'\',name,ext);
%     fprintf('TRATASTE DE CARGAR UN DEFAULT QUE NO EXISTE, usando path local ')
%     %Answer0.GroupCh{1}=Rsn(original);
%      end
%     
%     clear pathStr name ext
% end
Title = 'The Turbulent Brain: model-based';

%%%% SETTING DIALOG OPTIONS
% Options.WindowStyle = 'modal';
Options.Resize = 'on';
Options.Interpreter = 'tex';
Options.CancelButton = 'on';
% Options.ApplyButton = 'on';
Options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration
Option.Dim = 6; % Horizontal dimension in fields
Options.AlignControls='on';
Options.FontSize=9;

size_control=[150 50];
size_control_short=[150 50];

Prompt = {};
Formats = {};
DefAns = struct([]);





Prompt(end+1,:) = {'Folder name', 'saveFile',[]};
Formats(1,3).type = 'edit';
Formats(1,3).format = 'text';
Formats(1,3).size = 100; % automatically assign the height
DefAns(1).saveFile = 'TTB measures ';
% 
% Prompt(end+1,:) = {['Cargar la tseries'],[],[]};
% Formats(2,1).type = 'text';
% Formats(2,1).size = [-1 0];
% Formats(2,1).span = [1 1]; % item is 1 field x 4 fields




% 
% Prompt(end+1,:) = {['Cargar la SC'],[],[]};
% Formats(3,1).labelloc = 'leftmiddle';
% Formats(3,1).type = 'text';
% Formats(3,1).size = [-1 0];
% Formats(3,1).span = [1 1]; % item is 1 field x 4 fields

Prompt(end+1,:) = {'Load Structural Connectivity (mat file)','',''};
Formats(3,2).type = 'button';
Formats(3,2).size = [150 40]; % Tama�o
%Formats(5,1).callback = @(~,~,handles,k)msgbox(sprintf('You just pressed %s button',get(handles(k),'String')),'modal');
Formats(3,2).callback = @(~,~,h,k)dataload(h,k);



Prompt(end+1,:) = {'Choose SC','SCCh',[]};
Formats(3,3).labelloc = 'leftmiddle';
Formats(3,3).type = 'list';
Formats(3,3).style = 'list';
Formats(3,3).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(3,3).items = {'Structural Connectivity is not included'};
Formats(3,3).limits = [1 1]; % multi-select
Formats(3,3).size =  size_control;
DefAns.SCCh = {'Structural Connectivity is not included'};


Prompt(end+1,:) = {'SC Folder','SCDir',[]};
Formats(3,4).labelloc = 'topcenter';
Formats(3,4).type = 'edit';
Formats(3,4).format = 'text';
Formats(3,4).size = [-1 0];
Formats(3,4).span = [1 1];  % item is 1 field x 3 fields
DefAns.SCDir = pwd;



% 
% % 
% 
Prompt(end+1,:) = {'Which connection?','connectivity',[]};
Formats(3,5).labelloc = 'topcenter';
Formats(3,5).type = 'list';
Formats(3,5).style = 'list';
Formats(3,5).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(3,5).items = {'EDR','EDR+LR','SC'};
Formats(3,5).limits = [1 1]; % multi-select
Formats(3,5).size =  size_control;
DefAns.connectivity = {'EDR+LR'};
% 

Prompt(end+1,:) = {'Enable bandpass filter' 'EnableFilterMode',[]};
Formats(3,6).labelloc = 'topcenter';
Formats(3,6).type = 'check';
DefAns.EnableFilterMode = true;



Prompt(end+1,:) = {'Load RSN (mat file)','',''};
Formats(4,2).type = 'button';
Formats(4,2).size = [150 40]; % Tama�o
%Formats(2,1).callback = @(~,~,handles,k)msgbox(sprintf('You just pressed %s button',get(handles(k),'String')),'modal');
Formats(4,2).callback = @(~,~,h,k)dataload(h,k);


Prompt(end+1,:) = {'Choose RSN','RSNCh',[]};
Formats(4,3).labelloc = 'leftmiddle';
Formats(4,3).type = 'list';
Formats(4,3).style = 'list';
Formats(4,3).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(4,3).items = {'RSN belonging is not included'};
Formats(4,3).limits = [1 1]; % multi-select
Formats(4,3).size =  size_control;
DefAns.RSNCh = {'RSN belonging is not included'};
% 
Prompt(end+1,:) = {'RSN folder','RSNDir',[]};
Formats(4,4).labelloc = 'topcenter';
Formats(4,4).type = 'edit';
Formats(4,4).format = 'text';
Formats(4,4).size = [-1 0];
Formats(4,4).span = [1 1];  % item is 1 field x 3 fields
DefAns.RSNDir = pwd;



Prompt(end+1,:) = {'Lower Freq. bandpass', 'BPlb',[]};
Formats(4,5).labelloc = 'topcenter';
Formats(4,5).type = 'edit';
Formats(4,5).format = 'float';
Formats(4,5).size = 80;
Formats(4,5).limits = [0 inf] ;% non-negative decimal number
DefAns.BPlb = 0.04;


Prompt(end+1,:) = {'Upper Freq. bandpass', 'BPub',[]};
Formats(4,6).labelloc = 'topcenter';
Formats(4,6).type = 'edit';
Formats(4,6).format = 'float';
Formats(4,6).size = 80;
Formats(4,6).limits = [0 inf] ;% non-negative decimal number
DefAns.BPub = 0.07;




Prompt(end+1,:) = {'Load Centre of mass (mat file)','',''};
Formats(5,2).type = 'button';
Formats(5,2).size = [150 40]; % Tama�o
%Formats(2,1).callback = @(~,~,handles,k)msgbox(sprintf('You just pressed %s button',get(handles(k),'String')),'modal');
Formats(5,2).callback = @(~,~,h,k)dataload(h,k);

% 
Prompt(end+1,:) = {'Choose CoG ','GroupCoG',[]};
Formats(5,3).labelloc = 'leftmiddle';
Formats(5,3).type = 'list';
Formats(5,3).style = 'list';
Formats(5,3).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(5,3).items = {'Centre of mass is not included'};
Formats(5,3).limits = [1 1]; % multi-select
Formats(5,3).size =  size_control;
DefAns.GroupCoG = {'Centre of mass is not included'};

Prompt(end+1,:) = {'CoG folder','CoGDir',[]};
Formats(5,4).labelloc = 'topcenter';
Formats(5,4).type = 'edit';
Formats(5,4).format = 'text';
Formats(5,4).size = [-1 0];
Formats(5,4).span = [1 1];  % item is 1 field x 3 fields
DefAns.CoGDir = pwd;


Prompt(end+1,:) = {'Model ts up to...', 'Tsim','(0: max)'};
Formats(5,5).labelloc = 'topcenter';
Formats(5,5).type = 'edit';
Formats(5,5).format = 'integer';
Formats(5,5).limits = [0 999999999]; % 9-digits (positive #)
Formats(5,5).size = 80;
Formats(5,5).unitsloc = 'bottomleft';
DefAns.Tsim = 200;


Prompt(end+1,:) = {'Sample Rate (TR)', 'TRsec',[]};
Formats(5,6).labelloc = 'topcenter';
Formats(5,6).type = 'edit';
Formats(5,6).format = 'float';
Formats(5,6).size = 80;
Formats(5,6).limits = [0 inf] ;% non-negative decimal number
DefAns.TRsec = 2.6;


Prompt(end+1,:) = {'How many brain states','nBS',[]};
Formats(6,3).labelloc = 'topcenter';
Formats(6,3).type = 'list';
Formats(6,3).style = 'list';
Formats(6,3).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(6,3).items = {'1','2','3','4'};
Formats(6,3).limits = [1 1]; % multi-select
Formats(6,3).size =  size_control;
DefAns.nBS = {'2'};

Prompt(end+1,:) = {'Initial Lambda', 'iLambda',[]};
Formats(6,4).labelloc = 'topcenter';
Formats(6,4).type = 'edit';
Formats(6,4).format = 'float';
Formats(6,4).limits = [0 999999999]; % 9-digits (positive #)
Formats(6,4).size = 80;
Formats(6,4).unitsloc = 'bottomleft';
DefAns.iLambda = 0.01;

Prompt(end+1,:) = {'final Lambda', 'fLambda',[]};
Formats(6,5).labelloc = 'topcenter';
Formats(6,5).type = 'edit';
Formats(6,5).format = 'float';
Formats(6,5).size = 80;
Formats(6,5).limits = [0 inf] ;% non-negative decimal number
DefAns.fLambda = 0.3;


Prompt(end+1,:) = {'Lambda steps', 'stepsLam',[]};
Formats(6,6).labelloc = 'topcenter';
Formats(6,6).type = 'edit';
Formats(6,6).format = 'float';
Formats(6,6).size = 80;
Formats(6,6).limits = [0 inf] ;% non-negative decimal number
DefAns.stepsLam = 0.03;


Prompt(end+1,:) = {'G range lower limit', 'GLb',[]};
Formats(7,2).labelloc = 'topcenter';
Formats(7,2).type = 'edit';
Formats(7,2).format = 'float';
Formats(7,2).size = 80;
Formats(7,2).limits = [-inf inf]; % non-negative decimal number
DefAns.GLb = 0;
% 

Prompt(end+1,:) = {'G range upper limit', 'GUb',[]};
Formats(7,3).labelloc = 'topcenter';
Formats(7,3).type = 'edit';
Formats(7,3).format = 'float';
Formats(7,3).size = 80;
Formats(7,3).limits = [-inf inf] ;% non-negative decimal number
DefAns.GUb = 3;


Prompt(end+1,:) = {'G range step', 'Gstep',[]};
Formats(7,4).labelloc = 'topcenter';
Formats(7,4).type = 'edit';
Formats(7,4).format = 'float';
Formats(7,4).size = 80;
Formats(7,4).limits = [-inf inf] ;% non-negative decimal number
DefAns.Gstep = 0.25;

Prompt(end+1,:) = {'NSUBSIM (if =0,= NSUB)', 'NSIM',[]};
Formats(7,5).labelloc = 'topcenter';
Formats(7,5).type = 'edit';
Formats(7,5).format = 'integer';
Formats(7,5).limits = [0 999999999]; % 9-digits (positive #)
Formats(7,5).size = 80;
Formats(7,5).unitsloc = 'bottomleft';
DefAns.NSIM = 5;
% 

Prompt(end+1,:) = {'Model Plots?' 'PlotYes',[]};
Formats(7,6).labelloc = 'topcenter';
Formats(7,6).type = 'check';
DefAns.PlotYes = true;

% perturbations

Prompt(end+1,:) = {'Do you want model perturbation?' 'PertuYes',[]};
Formats(8,2).labelloc = 'topcenter';
Formats(8,2).type = 'check';
DefAns.PertuYes = true;
Formats(8,2).callback = @(~,~,h,k)disablePerts(h,k);


Prompt(end+1,:) = {'Statistical Test','stattest',[]};
Formats(8,3).labelloc = 'topcenter';
Formats(8,3).type = 'list';
Formats(8,3).style = 'list';
Formats(8,3).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(8,3).items = {'ttest','ttest2','signrank','ranksum','permutation'};
Formats(8,3).limits = [1 1]; % multi-select
Formats(8,3).size =  size_control;
DefAns.stattest = {'ranksum'};



Prompt(end+1,:) = {'Lambda perturbed', 'lamICS',[]};
Formats(8,5).labelloc = 'topcenter';
Formats(8,5).type = 'edit';
Formats(8,5).format = 'float';
Formats(8,5).size = 80;
Formats(8,5).limits = [-inf inf] ;% non-negative decimal number
DefAns.lamICS=0.18;


Prompt(end+1,:) = {'Trials', 'nTrials',[]};
Formats(8,6).labelloc = 'topcenter';
Formats(8,6).type = 'edit';
Formats(8,6).format = 'float';
Formats(8,6).size = 80;
Formats(8,6).limits = [-inf inf] ;% non-negative decimal number
DefAns.nTrials=10;

[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);

Answer0=Answer;
% 
% %{'FC prom en Fisher';'FC por suj';'FCD concat';'FCD por suj';'Ventana de FCD';'FMI'};
% switch Answer.ObsCh{1}
%     case 'FC prom con Fisher'
%         Answer.fx_emp=@(tseries) observable_FC(tseries);
%         Answer.fx_TS=@(tseries) obtain_TS_nochange(tseries);
%         Answer.fx_obs=@(tseries) observable_FC(tseries);
%         
%     case 'FC por suj'
%         Answer.fx_emp=@(tseries) observable_FCsev(tseries);
%         Answer.fx_TS=@(tseries) obtain_TS_nochange(tseries);
%         Answer.fx_obs=@(tseries) observable_FC(tseries);
%     case 'FCD concat'
%         Answer.fx_emp=@(tseries) observable_FCD(tseries,60,40,Answer.TRsec);
%         Answer.fx_TS=@(tseries) obtain_TS_nochange(tseries);
%         Answer.fx_obs=@(tseries) observable_FCD(tseries,60,40,Answer.TRsec);
%     case 'FCD por suj'
%         Answer.fx_emp=@(tseries) observable_FCDsev(tseries,60,40,Answer.TRsec);
%         Answer.fx_TS=@(tseries) obtain_TS_nochange(tseries);
%         Answer.fx_obs=@(tseries) observable_FCD(tseries,60,40,Answer.TRsec);
%     case 'FCD mean'
%         Answer.fx_emp=@(tseries) observable_FCDmean(tseries,60,40,Answer.TRsec);
%         Answer.fx_TS=@(tseries) obtain_TS_nochange(tseries);
%         Answer.fx_obs=@(tseries) observable_FCD(tseries,60,40,Answer.TRsec);
%     case 'Ventana de FCD'
%         Answer.fx_emp=@(tseries) observable_FCDvent(tseries,60,40,Answer.TRsec);
%         Answer.fx_TS=@(tseries) obtain_TS_fct(tseries,60,40,Answer.TRsec);
%         Answer.fx_obs=@(tseries) observable_FC(tseries);
%     case 'FMI'
%         Answer.fx_emp=@(tseries) observable_FMI(tseries);
%         Answer.fx_TS=@(tseries) obtain_TS_nochange(tseries);
%         Answer.fx_obs=@(tseries) observable_FMI(tseries);
%     case 'GAMMA'
%         Answer.fx_emp=@(tseries) observable_FCDmean(tseries,60,40,Answer.TRsec);
%         Answer.fx_TS=@(tseries) obtain_TS_nochange(tseries);
%         Answer.fx_obs=@(tseries) observable_fctaFCD(tseries);
% end
% 
% switch Answer.MetCh{1}
%     
%     case '1-SSIM'
%         Answer.fx_metric = @(FSim,FEmp) 1-ssim(FSim,FEmp);
%     case 'Frobenius'
%         Answer.fx_metric = @(FSim,FEmp) norm(tril(FEmp,-1)-tril(FSim,-1),'fro');
%     case 'Kolmogorov-Smirnov'
%         Answer.fx_metric = @(FSim,FEmp) metrica_KS(FSim,FEmp);
% end
% 
% 
% switch Answer.ParcellCh{1}
%     
%     case 'aes por grupo, G=G_Ini'
%         Answer.Parcell = 1;
%     case 'aes por grupo, G homog'
%         Answer.Parcell = 2;
%     case 'a por grupo, G por grupo'
%         Answer.Parcell = 3;
%     case 'a homog, G por grupo'
%         Answer.Parcell = 4;
% end
% 
% switch Answer.aIni{1}
%     
%     case 'zeros(nNodes,1)'
%         Answer.fx_aIni=@(nNodes)Answer.aesIniAmp*zeros(nNodes,1);
%     case 'rand(nNodes,1)'
%         Answer.fx_aIni=@(nNodes)Answer.aesIniAmp*rand(nNodes,1);
%     case 'ones(nNodes,1)'
%         Answer.fx_aIni=@(nNodes)Answer.aesIniAmp*ones(nNodes,1);
%         
% end
% 
% 
% switch Answer.GIni{1}
%     
%     case 'zeros(nNodes,1)'
%         Answer.fx_GIni=@(nNodes)Answer.GIniAmp*zeros(nNodes,1);
%     case 'rand(nNodes,1)'
%         Answer.fx_GIni=@(nNodes)Answer.GIniAmp*rand(nNodes,1);
%     case 'ones(nNodes,1)'
%         Answer.fx_GIni=@(nNodes)Answer.GIniAmp*ones(nNodes,1);
%         
% end
%         
%         
%     
% 
