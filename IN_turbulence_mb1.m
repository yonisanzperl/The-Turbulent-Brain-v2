
% Si se buguea la GUI
% set(0,'ShowHiddenHandles','on')
% delete(get(0,'Children'))

function [Answer,Answer0,Cancelled]=IN_turbulence_mb1(varargin)
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




% 

Prompt(end+1,:) = {'Folder name', 'saveFile',[]};
Formats(1,3).type = 'edit';
Formats(1,3).format = 'text';
Formats(1,3).size = 100; % automatically assign the height
DefAns(1).saveFile = 'TTB measures';
% 
% Prompt(end+1,:) = {['Cargar la tseries'],[],[]};
% Formats(2,1).type = 'text';
% Formats(2,1).size = [-1 0];
% Formats(2,1).span = [1 1]; % item is 1 field x 4 fields



Prompt(end+1,:) = {'Load timeseires (mat file)','',''};
Formats(2,2).type = 'button';
Formats(2,2).size = [150 40]; % Tamaño
%Formats(2,1).callback = @(~,~,handles,k)msgbox(sprintf('You just pressed %s button',get(handles(k),'String')),'modal');
Formats(2,2).callback = @(~,~,h,k)dataload(h,k);
 

Prompt(end+1,:) = {'Choose timeseries:','TSCh',[]};
Formats(2,3).labelloc = 'leftmiddle';
Formats(2,3).type = 'list';
Formats(2,3).style = 'list';
Formats(2,3).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(2,3).items = {'TS_W' ; 'there is not ts'};
Formats(2,3).limits = [1 1]; % multi-select
Formats(2,3).size = size_control;
DefAns.TSCh = {'there is not ts'};


Prompt(end+1,:) = {'Timeseries Folder','TSDir',[]};
Formats(2,4).labelloc = 'topcenter';
Formats(2,4).type = 'edit';
Formats(2,4).format = 'text';
Formats(2,4).size = [-1 0];
Formats(2,4).span = [1 1];  % item is 1 field x 3 fields
DefAns.TSDir = pwd;



Prompt(end+1,:) = {'Load output from model-based approach (mat file)','',''};
Formats(3,2).type = 'button';
Formats(3,2).size = [150 40]; % Tamaño
%Formats(5,1).callback = @(~,~,handles,k)msgbox(sprintf('You just pressed %s button',get(handles(k),'String')),'modal');
Formats(3,2).callback = @(~,~,h,k)dataload(h,k);



Prompt(end+1,:) = {'View variables (will load all of them)','allIt',[]};
Formats(3,3).labelloc = 'leftmiddle';
Formats(3,3).type = 'list';
Formats(3,3).style = 'list';
Formats(3,3).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(3,3).items = {'SC','there is not SC'};
Formats(3,3).limits = [1 1]; % multi-select
Formats(3,3).size =  size_control;
DefAns.allIt = {'there is not SC'};


Prompt(end+1,:) = {'SC Folder','MFoutdir',[]};
Formats(3,4).labelloc = 'topcenter';
Formats(3,4).type = 'edit';
Formats(3,4).format = 'text';
Formats(3,4).size = [-1 0];
Formats(3,4).span = [1 1];  % item is 1 field x 3 fields
DefAns.MFoutdir = pwd;



Prompt(end+1,:) = {'Model ts up to...', 'Tsim','(0: max)'};
Formats(3,5).labelloc = 'topcenter';
Formats(3,5).type = 'edit';
Formats(3,5).format = 'integer';
Formats(3,5).limits = [0 999999999]; % 9-digits (positive #)
Formats(3,5).size = 80;
Formats(3,5).unitsloc = 'bottomleft';
DefAns.Tsim = 200;



Prompt(end+1,:) = {'How many brain states','nBS',[]};
Formats(4,3).labelloc = 'topcenter';
Formats(4,3).type = 'list';
Formats(4,3).style = 'list';
Formats(4,3).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(4,3).items = {'1','2','3','4','5','6'};
Formats(4,3).limits = [1 1]; % multi-select
Formats(4,3).size =  size_control;
DefAns.nBS = {'2'};
%Formats(7,2).callback = @(~,~,h,k)disableG(h,k);

Prompt(end+1,:) = {'Initial Lambda', 'iLambda',[]};
Formats(4,4).labelloc = 'topcenter';
Formats(4,4).type = 'edit';
Formats(4,4).format = 'float';
Formats(4,4).limits = [0 999999999]; % 9-digits (positive #)
Formats(4,4).size = 80;
Formats(4,4).unitsloc = 'bottomleft';
DefAns.iLambda = 0.01;

Prompt(end+1,:) = {'final Lambda', 'fLambda',[]};
Formats(4,5).labelloc = 'topcenter';
Formats(4,5).type = 'edit';
Formats(4,5).format = 'float';
Formats(4,5).size = 80;
Formats(4,5).limits = [0 inf] ;% non-negative decimal number
DefAns.fLambda = 0.3;


Prompt(end+1,:) = {'Lambda steps', 'stepsLam',[]};
Formats(4,6).labelloc = 'topcenter';
Formats(4,6).type = 'edit';
Formats(4,6).format = 'float';
Formats(4,6).size = 80;
Formats(4,6).limits = [0 inf] ;% non-negative decimal number
DefAns.stepsLam = 0.03;



%% modelado





% 
% 
% 
Prompt(end+1,:) = {'G range lower limit', 'GLb',[]};
Formats(5,2).labelloc = 'topcenter';
Formats(5,2).type = 'edit';
Formats(5,2).format = 'float';
Formats(5,2).size = 80;
Formats(5,2).limits = [-inf inf]; % non-negative decimal number
DefAns.GLb = 0;
% 

Prompt(end+1,:) = {'G range upper limit', 'GUb',[]};
Formats(5,3).labelloc = 'topcenter';
Formats(5,3).type = 'edit';
Formats(5,3).format = 'float';
Formats(5,3).size = 80;
Formats(5,3).limits = [-inf inf] ;% non-negative decimal number
DefAns.GUb = 3;


Prompt(end+1,:) = {'G range step', 'Gstep',[]};
Formats(5,4).labelloc = 'topcenter';
Formats(5,4).type = 'edit';
Formats(5,4).format = 'float';
Formats(5,4).size = 80;
Formats(5,4).limits = [-inf inf] ;% non-negative decimal number
DefAns.Gstep = 0.25;

Prompt(end+1,:) = {'NSUBSIM (if =0,= NSUB)', 'NSIM',[]};
Formats(5,5).labelloc = 'topcenter';
Formats(5,5).type = 'edit';
Formats(5,5).format = 'integer';
Formats(5,5).limits = [0 999999999]; % 9-digits (positive #)
Formats(5,5).size = 80;
Formats(5,5).unitsloc = 'bottomleft';
DefAns.NSIM = 5;
% 

Prompt(end+1,:) = {'Model Plots?' 'PlotYes',[]};
Formats(5,6).labelloc = 'topcenter';
Formats(5,6).type = 'check';
DefAns.PlotYes = true;

% perturbations

Prompt(end+1,:) = {'Do you want model perturbation?' 'PertuYes',[]};
Formats(7,3).labelloc = 'topcenter';
Formats(7,3).type = 'check';
DefAns.PertuYes = true;
Formats(7,3).callback = @(~,~,h,k)disablePerts(h,k);



Prompt(end+1,:) = {'Lambda perturbed', 'lamICS',[]};
Formats(7,4).labelloc = 'topcenter';
Formats(7,4).type = 'edit';
Formats(7,4).format = 'float';
Formats(7,4).size = 80;
Formats(7,4).limits = [-inf inf] ;% non-negative decimal number
DefAns.lamICS=0.18;


Prompt(end+1,:) = {'Trials', 'nTrials',[]};
Formats(7,5).labelloc = 'topcenter';
Formats(7,5).type = 'edit';
Formats(7,5).format = 'float';
Formats(7,5).size = 80;
Formats(7,5).limits = [-inf inf] ;% non-negative decimal number
DefAns.nTrials=10;

% Prompt(end+1,:) = {'Selección de a inicial:','aIni',[]};
% Formats(8,2).labelloc = 'topcenter';
% Formats(8,2).type = 'list';
% Formats(8,2).style = 'list';
% Formats(8,2).format = 'text'; % Answer will give value shown in items, disable to get integer
% Formats(8,2).items = {'zeros(nNodes,1)','rand(nNodes,1)','ones(nNodes,1)'};
% Formats(8,2).limits = [1 1]; % multi-select
% Formats(8,2).size =  size_control_short;
% DefAns.aIni = {'rand(nNodes,1)'};
% Formats(8,2).callback = @(~,~,h,k)disableaIni(h,k);
% 
% 
% Prompt(end+1,:) = {'Amplitud aIni', 'aesIniAmp',[]};
% Formats(8,3).labelloc = 'topcenter';
% Formats(8,3).type = 'edit';
% Formats(8,3).format = 'float';
% Formats(8,3).size = 80;
% Formats(8,3).limits = [-inf inf] ;% non-negative decimal number
% DefAns.aesIniAmp = 0;
% 
% 
% Prompt(end+1,:) = {'Selección de G inicial:','GIni',[]};
% Formats(8,4).labelloc = 'topcenter';
% Formats(8,4).type = 'list';
% Formats(8,4).style = 'list';
% Formats(8,4).format = 'text'; % Answer will give value shown in items, disable to get integer
% Formats(8,4).items = {'zeros(nNodes,1)','rand(nNodes,1)','ones(nNodes,1)'};
% Formats(8,4).limits = [1 1]; % multi-select
% Formats(8,4).size =  size_control_short;
% DefAns.GIni = {'rand(nNodes,1)'};
% Formats(8,4).callback = @(~,~,h,k)disableaIni(h,k);
% 
% 
% Prompt(end+1,:) = {'Amplitud G inicial', 'GIniAmp',[]};
% Formats(8,5).labelloc = 'topcenter';
% Formats(8,5).type = 'edit';
% Formats(8,5).format = 'float';
% Formats(8,5).size = 80;
% Formats(8,5).limits = [-inf inf] ;% non-negative decimal number
% DefAns.GIniAmp = 0.5;
% 
% 
% 
% Prompt(end+1,:) = {'Numero de runs por fitness', 'Prom',''};
% Formats(9,2).labelloc = 'topcenter';
% Formats(9,2).type = 'edit';
% Formats(9,2).format = 'integer';
% Formats(9,2).limits = [1 999999999]; % 9-digits (positive #)
% Formats(9,2).size = 80;
% Formats(9,2).unitsloc = 'bottomleft';
% DefAns.Prom = 1;
% 
% 
% Prompt(end+1,:) = {'Paso dt', 'dt',[]};
% Formats(9,3).labelloc = 'topcenter';
% Formats(9,3).type = 'edit';
% Formats(9,3).format = 'float';
% Formats(9,3).size = 80;
% Formats(9,3).limits = [-inf inf] ;% non-negative decimal number
% DefAns.dt = 0.1;
% 
% 
% Prompt(end+1,:) = {'Enable permuta' 'EnablePermutaMode',[]};
% Formats(9,4).labelloc = 'topcenter';
% Formats(9,4).type = 'check';
% DefAns.EnablePermutaMode = false;
% 
% 
% Prompt(end+1,:) = {'Enable SC.*sign(FEmp)' 'EnableAnticorrelMode',[]};
% Formats(9,5).labelloc = 'topcenter';
% Formats(9,5).type = 'check';
% DefAns.EnablePermutaMode = false;
% 


%% Fin cosas particulares del integrador

% exist('Answer0')
% if exist('Answer0')
%     
%     DefAns=Answer0;
%     %Necesito agregarle las listas a los items
%     Formats(2,3).items = who( matfile( DefAns.TSDir ) );
%     Formats(3,3).items = who( matfile( DefAns.SCDir ) );    
% end
% 
% save('errorinput','DefAns','Formats')


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
