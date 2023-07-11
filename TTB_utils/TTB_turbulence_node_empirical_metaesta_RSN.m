function [output]=TTB_turbulence_node_empirical_metaesta_RSN(Cfg, CoG,SC,TS,RSN,Clong, lambda_emp, numexcSC, valexcSC);

% Parameters of the data
xs=TS(:);                         %
NSUB=size(find(~cellfun(@isempty,xs)),1);   
TR=Cfg.TR;                                  % Repetition Time (seconds)
NPARCELLS=Cfg.nNodes;                       % Atlas schaefer 1000 nodes
NR=Cfg.NR;
NRini=Cfg.NRini;
NRfin=Cfg.NRfin;

labels=RSN;

LAMBDA = Cfg.iLambda:Cfg.stepsLam:Cfg.fLambda;
LAMBDA = flip(LAMBDA);
NLAMBDA=length(LAMBDA);

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = Cfg.filt.lb;                  % lowpass frequency of filter (Hz)
fhi = Cfg.filt.ub;                   % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter




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


C1=zeros(NLAMBDA,NPARCELLS,NPARCELLS);
[aux indsca]=min(abs(LAMBDA-lambda_emp));
ilam=1;
for lambda=LAMBDA
    for i=1:NPARCELLS
        for j=1:NPARCELLS
            C1(ilam,i,j)=exp(-lambda*rr(i,j));
        end
    end
    ilam=ilam+1;
end



TransferLambda_sub=zeros(NLAMBDA,NSUB);
Turbulence_sub=zeros(NLAMBDA,NSUB);
node_Turbulence_sub=zeros(NLAMBDA,NSUB,NPARCELLS);
InformationCascade_sub=zeros(1,NSUB);
Transfer_sub=zeros(1,NSUB);
fcr=zeros(NLAMBDA,NSUB);
fclam=zeros(NLAMBDA,NPARCELLS,NPARCELLS);



for sub=1:NSUB
    fprintf('Subject Number: %d/ %d \n',sub,NSUB)    ;
    ts=TS{sub};
    if Cfg.Tmax>0
        Tmax=Cfg.Tmax;                                   % Timepoints
    else
        Tmax= size(ts{1},2);
    end
    
    ts = ts(:,1:Tmax);
    clear Phases Xanalytic signal_filt
    enstrophy=zeros(NLAMBDA,NPARCELLS,Tmax);
    signal_filt=zeros(NPARCELLS,Tmax);
    Phases=zeros(NPARCELLS,Tmax);

   
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        
        if sum(isnan(ts(seed,:)))<1
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        else
            ts(seed,:)=ts(seed,:);
            signal_filt(seed,:) = ts(seed,:);
        end
        
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end

    for i=1:NPARCELLS
       
        for ilam=1:NLAMBDA
            C1lam=squeeze(C1(ilam,:,:));
            sumphases=nansum(repmat(C1lam(i,:)',1,Tmax).*complex(cos(Phases),sin(Phases)))/nansum(C1lam(i,:));
            enstrophy1(ilam,i,:)=abs(sumphases);
        end
    end
    TransferLambda_sub(1,:)=NaN;

    for ii=1:max(labels)
        indx_rsn=find(labels(:)==ii);
%        NPARCELLS=size(indx_rsn,1);
        clear enstrophy
        for ilam=1:NLAMBDA
            enstrophy(ilam,:,:) = enstrophy1(ilam,indx_rsn,:);
            TurbulenceRSN_sub(ilam,sub,ii)=nanstd(squeeze(enstrophy(ilam,:)));
        end
    end     
    % compute Turbe and node_tubr
    for ilam=1:NLAMBDA
        Turbulence_sub(ilam,sub)=nanstd(squeeze(enstrophy1(ilam,:)));
        for j=1:NPARCELLS
            node_Turbulence_sub(ilam,sub,j) = nanstd(squeeze(enstrophy1(ilam,j,:)));
        end
    end
    
    %metaestability
    gKoP(sub)=nanmean(abs(nansum(complex(cos(Phases(:,:)),sin(Phases(:,:))),1))/NPARCELLS);
    Meta(sub) = nanstd(abs(nansum(complex(cos(Phases(:,:)),sin(Phases(:,:))),1))/NPARCELLS);


    for ilam=1:NLAMBDA-1
        [cc pp]=corr(squeeze(enstrophy(ilam+1,:,2:end))',squeeze(enstrophy(ilam,:,1:end-1))');
        TransferLambda_sub(ilam+1,sub)=nanmean(abs(cc(find(pp(:)<0.05))));
    end
    
    InformationCascade_sub(sub)=nanmean(TransferLambda_sub(2:NLAMBDA,sub),1);
    
    %%% Transfer across space
    clear fclam
    for ilam=1:NLAMBDA
        fclam(ilam,:,:)=corrcoef(squeeze(enstrophy1(ilam,:,:))');
    end
    
    
    for lam=1:NLAMBDA
        numind=zeros(1,NR);
        fcra=zeros(1,NR);
        for i=1:NPARCELLS
            for j=1:NPARCELLS
                r=rr(i,j);
                index=floor(r/delta)+1;
                if index==NR+1
                    index=NR;
                end
                mcc=fclam(lam,i,j);
                if ~isnan(mcc)
                    fcra(index)=fcra(index)+mcc;
                    numind(index)=numind(index)+1;
                end
            end
        end
        %%% Powerlaw a
        
        grandcorrfcn=fcra./numind;
        clear xcoor;
        clear ycoor;
        nn=1;
        indxx = find(~isnan(grandcorrfcn));
        indices = indxx(indxx>NRini & indxx<NRfin);
        for k=1:length(indices)
            if grandcorrfcn(k)>0
                xcoor(nn)=log(xrange(k));
                ycoor(nn)=log(grandcorrfcn(k)/grandcorrfcn(indices(1)));
                nn=nn+1;
            end
            
        end
        linfunc = @(A, x)(A(1)*x+A(2));
        options=optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');
        A0=[-1 1];
        [Afit Residual]= lsqcurvefit(linfunc,A0,xcoor,ycoor,[-4 -10],[4 10],options);
        Transfer_sub(lam,sub)=abs(Afit(1));       
    end
end
InformationCascade_sub=nanmean(TransferLambda_sub,1);


output.LAMBDA= LAMBDA;
output.TurbulenceRSN_sub=TurbulenceRSN_sub;
output.Turbulence_sub=Turbulence_sub;
output.node_Turbulence_sub=node_Turbulence_sub;
output.TurbulenceRSN_sub=TurbulenceRSN_sub;
output.Transfer_sub=Transfer_sub;
output.InformationCascade_sub=InformationCascade_sub;
output.TransferLambda_sub=TransferLambda_sub;
output.gKoP=gKoP;
output.Meta=Meta;
output.xcoor=xcoor;
output.ycoor=ycoor;

%save (sprintf('turbu_measurements%d_RSN_TEST.mat',cond), 'LAMBDA','TurbulenceRSN_sub','node_Turbulence_sub', 'TurbulenceRSN_sub','Transfer_sub', 'InformationCascade_sub', 'TransferLambda_sub','gKoP','Meta');
