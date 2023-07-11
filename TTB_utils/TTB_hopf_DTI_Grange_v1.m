function [err_hete model_output]=TTB_hopf_DTI_Grange_v1(Cfg,TS,CoG,SC,corrfcn,f_diff,RSN,Clong, lambda_emp, numexcSC, valexcSC);


xs=TS; 

NPARCELLS=Cfg.nNodes;
NR=Cfg.NR;
NRini=Cfg.NRini;
NRfin=Cfg.NRfin;
NSUBSIM=Cfg.NSIM;


if NSUBSIM==0
    NSUBSIM=NSUB;
end
lambda=round(lambda_emp,2);
Tmax=Cfg.Tmax;

G_range=Cfg.Glower:Cfg.Gstep:Cfg.Gupper;
empcorrfcn=corrfcn;

rr=zeros(NPARCELLS,NPARCELLS);
for i=1:NPARCELLS
    for j=1:NPARCELLS
        rr(i,j)=norm(CoG(i,:)-CoG(j,:));
    end
end
range=max(max(rr));
delta=range/NR;

for i=1:NR
    xrange(i)=delta/2+delta*(i-1);
end

Isubdiag = find(tril(ones(NPARCELLS),-1));

C=zeros(NPARCELLS,NPARCELLS);


%%%
% Parameters of the data
TR=Cfg.TR;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = Cfg.filt.lb;                    % lowpass frequency of filter (Hz)
fhi = Cfg.filt.ub;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Isubdiag = find(tril(ones(NPARCELLS),-1));

% Parameters HOPF
omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
dt=0.1*TR/2;
sig=0.01;
dsig = sqrt(dt)*sig;

%%

corrfcn=zeros(NPARCELLS,NR);
err_hete=zeros(NSUBSIM,size(G_range,2));


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

for nG=1:size(G_range,2)
    G = G_range(nG);
    fprintf('\n \n Coupling Strength G: %f \n',G)    ;
    for sub=1:NSUBSIM
        
        
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
        ts_out{sub}=xs';
        ts = xs';
        for seed=1:NPARCELLS
            ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
            Xanalytic = hilbert(demean(signal_filt(seed,:)));
            Phases(seed,:) = angle(Xanalytic);
        end
        fcsimul=corrcoef(signal_filt');
        
        for i=1:NPARCELLS
            numind=zeros(1,NR);
            corrfcn_1=zeros(1,NR);
            for j=1:NPARCELLS
                r=rr(i,j);
                index=floor(r/delta)+1;
                if index==NR+1
                    index=NR;
                end
                mcc=fcsimul(i,j);
                if ~isnan(mcc)
                    corrfcn_1(index)=corrfcn_1(index)+mcc;
                    numind(index)=numind(index)+1;
                end
            end
            corrfcn(i,:)=corrfcn_1./numind;
            
            
            for i=1:NPARCELLS
                for k=NRini:NRfin
                    err11(k)=(corrfcn(i,k)-empcorrfcn(i,k))^2;
                end
                err1(i)=(nanmean(err11(NRini:NRfin)));
            end
        
        err_hete(sub,nG)=sqrt(nanmean(err1));
        end
    end
        
   [model_output(nG)] = TTB_turbulence_node_empirical_metaesta_RSN(Cfg, CoG,SC,ts_out,RSN,Clong, lambda_emp, numexcSC, valexcSC);               
            
    end
    
