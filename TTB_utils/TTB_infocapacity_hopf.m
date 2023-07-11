function [infocapacity,susceptibility]=  infocapacity_hopf(ss,kk,G,Cfg,CoG,SC,f_diff,Clong, lambda_emp, numexcSC, valexcSC)



NPARCELLS=Cfg.nNodes;
NR=Cfg.NR;
NRini=Cfg.NRini;
NRfin=Cfg.NRfin;
NSUBSIM=Cfg.NSIM;

SchaeferCOG=CoG;



%lambda=round(lambda_emp,1);
lambda = Cfg.lamICS;

rr=zeros(NPARCELLS,NPARCELLS);
for i=1:NPARCELLS
    for j=1:NPARCELLS
        rr(i,j)=norm(SchaeferCOG(i,:)-SchaeferCOG(j,:));
    end
end
range=max(max(rr));
delta=range/NR;

for i=1:NR
    xrange(i)=delta/2+delta*(i-1);
end

Isubdiag = find(tril(ones(NPARCELLS),-1));

C=zeros(NPARCELLS,NPARCELLS);

LAMBDA = Cfg.iLambda:Cfg.stepsLam:Cfg.fLambda;
LAMBDA = flip(LAMBDA);




[aux indsca]=min(abs(LAMBDA-lambda));
LAMBDA = LAMBDA(indsca);
NLAMBDA=length(LAMBDA);
C1=zeros(NLAMBDA,NPARCELLS,NPARCELLS);

ilam=1;
for lambda2=LAMBDA
    for i=1:NPARCELLS
        for j=1:NPARCELLS
            C1(ilam,i,j)=exp(-lambda2*rr(i,j));
        end
    end
    ilam=ilam+1;
end

%%%
% Parameters of the data
TR=Cfg.TR;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELLS),-1));

% Parameters HOPF
Tmax=Cfg.Tmax;
omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
dt=0.1*TR/2;
sig=0.01;
dsig = sqrt(dt)*sig;

%%

lam_mean_spatime_enstrophy=zeros(NLAMBDA,NPARCELLS,Tmax);
ensspasub=zeros(NLAMBDA,NSUBSIM,NPARCELLS);
ensspasub1=zeros(NLAMBDA,NSUBSIM,NPARCELLS);


IClong=find(Clong>0);
for i=1:NPARCELLS
    for j=1:NPARCELLS
        C(i,j)=exp(-lambda*rr(i,j));
    end
    C(i,i)=0;
end
C(IClong)=Clong(IClong);

factor=max(max(C));
C=C/factor*0.2;

for s=1:ss
   fprintf('Subject Number: %d/ %d \n',s,ss)    ;
for sub=1:NSUBSIM
   % sub    
    wC = G*C;
    sumC = repmat(sum(wC,2),1,2);
    
    %% Hopf Simulation
    a=-0.02*ones(NPARCELLS,2);
    xs=zeros(Tmax,NPARCELLS);
    %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
    z = 0.1*ones(NPARCELLS,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    % discard first 2000 time steps
    for t=0:dt:2000
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
    end
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    
    for i=1:NPARCELLS
        %%% enstrophy
        ilam=1;
        for lam=LAMBDA
            enstrophy=nansum(repmat(squeeze(C1(ilam,i,:)),1,Tmax).*complex(cos(Phases),sin(Phases)))/sum(C1(ilam,i,:));
            lam_mean_spatime_enstrophy(ilam,i,:)=abs(enstrophy);
            Rspatime=squeeze(lam_mean_spatime_enstrophy(ilam,:,:));
            ensspasub(ilam,sub,:)=(nanmean(Rspatime,2))';
            ilam=ilam+1;
        end
    end

    
    %%% Perturbation
    
    a=-0.02+0.02*repmat(rand(NPARCELLS,1),1,2).*ones(NPARCELLS,2);
    nn=0;
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(NPARCELLS,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    ts=xs';
    Rspatime1=zeros(NPARCELLS,Tmax);
    
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end

    
        
    for i=1:NPARCELLS
        %%% enstrophy
        ilam=1;
        for lam=LAMBDA
            enstrophy=nansum(repmat(squeeze(C1(ilam,i,:)),1,Tmax).*complex(cos(Phases),sin(Phases)))/sum(C1(ilam,i,:));
            lam_mean_spatime_enstrophy(ilam,i,:)=abs(enstrophy);
            Rspatime=squeeze(lam_mean_spatime_enstrophy(ilam,:,:));
            ensspasub1(ilam,sub,:)=(nanmean(Rspatime,2))';
            ilam=ilam+1;
        end
    end

    
    
    
end
% for ii=1:NLAMBDA
%     infocapacity(s,ii,:)=nanmean(nanstd(squeeze(ensspasub1(ii,:,:))-ones(NSUBSIM,1)*nanmean(squeeze(ensspasub(ii,:,:)))));
%     susceptibility(s,ii,:)=nanmean(nanmean(squeeze(ensspasub1(ii,:,:))-ones(NSUBSIM,1)*nanmean(squeeze(ensspasub(ii,:,:)))));
% end
    ii=indsca;
    ii=1;

    infocapacity(s,:)=nanmean(nanstd(squeeze(ensspasub1(ii,:,:))-ones(NSUBSIM,1)*nanmean(squeeze(ensspasub(ii,:,:)))));
    susceptibility(s,:)=nanmean(nanmean(squeeze(ensspasub1(ii,:,:))-ones(NSUBSIM,1)*nanmean(squeeze(ensspasub(ii,:,:)))));

end
%save(sprintf('Wtrials_full_%d.mat',kk),'infocapacity','susceptibility');