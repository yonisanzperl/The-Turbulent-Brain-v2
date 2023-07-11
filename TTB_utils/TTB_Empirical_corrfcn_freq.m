function [corrfcn,f_diff,fce]=TTB_Empirical_corrfcn_freq(Cfg,SC,CoG,TSk);


NPARCELLS=Cfg.nNodes;
NR=Cfg.NR;

TR=Cfg.TR;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = Cfg.filt.lb;                    % lowpass frequency of filter (Hz)
fhi = Cfg.filt.ub;                    % highpass
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


xs=TSk;
NSUB=size(find(~cellfun(@isempty,xs)),1);
ensspasub=zeros(NSUB,NPARCELLS);
Rsub=zeros(1,NSUB);
DTsub=zeros(1,NSUB);
corrfcnsub=zeros(NSUB,NPARCELLS,NR);
for sub=1:NSUB
    
    ts=xs{sub,1};
    clear signal_filt 
    for seed=1:NPARCELLS
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        
        if sum(isnan(ts(seed,:)))<1
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        else
            ts(seed,:)=ts(seed,:);
            signal_filt(seed,:) = ts(seed,:);
        end
    end
    
    fce=corrcoef(signal_filt');
    
    for i=1:NPARCELLS
        numind=zeros(1,NR);
        corrfcn_1=zeros(1,NR);
        for j=1:NPARCELLS
            r=rr(i,j);
            index=floor(r/delta)+1;
            if index==NR+1
                index=NR;
            end
            mcc=fce(i,j);
            if ~isnan(mcc)
                corrfcn_1(index)=corrfcn_1(index)+mcc;
                numind(index)=numind(index)+1;
            end
        end
        corrfcnsub(sub,i,:)=corrfcn_1./numind;
    end
end

corrfcn=squeeze(nanmean(corrfcnsub));

ts=xs{sub,1};
[Ns, Tmax]=size(ts);
TT=Tmax;
Ts = TT*TR;
freq = (0:TT/2-1)/Ts;
nfreqs=length(freq);

for seed=1:NPARCELLS
    
    if sum(isnan(ts(seed,:)))<1
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        tss(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
    else
        tss(seed,:)=ts(seed,:);
    end
    
    
    pw = abs(fft(tss(seed,:)));
    PowSpect(:,seed,sub) = pw(1:floor(TT/2)).^2/(TT/TR);
end


Power_Areas=squeeze(mean(PowSpect,3));
for seed=1:NPARCELLS
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);
f_diff(find(f_diff==0))=mean(f_diff(find(f_diff~=0)));



end