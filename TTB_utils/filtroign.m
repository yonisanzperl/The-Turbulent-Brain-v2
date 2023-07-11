
function [new_tseries]=filtroign(old_tseries,TRsec,lb,ub)
    
    delt = TRsec;                               % sampling interval
    fnq = 1/(2*delt);                           % Nyquist frequency
    k = 4;                                      % 2nd order butterworth filter

    if fnq<ub
        error('La cagaste, la frec max es > que frec de Nyquist');
    end
    
    if (iscell(old_tseries))
      
        nSubs=length(old_tseries);
        nNodes=size(old_tseries{1},1);
        for i=1:nSubs
            
            
            each_Tmax(i)=size(old_tseries{i},2);
            
        end
    else
        nSubs=size(old_tseries,3);
        nNodes=size(old_tseries,1);
        for i=1:nSubs
                       
            each_Tmax(i)=size(old_tseries(:,:,i),2);
            
        end
        
    end
    %siempre chequea si es una celda, porque por ahi es solo un array
    %asume todos los sujetos tienen el mismo TRsec

    

for i=1:nSubs
    
    Tmax=each_Tmax(i);
    
    TT=Tmax;
    Ts = TT*TRsec;
    freq = (0:TT/2-1)/Ts;

    %NARROW LOW BANDPASS
    flp = lb;                                  % lowpass frequency of filter
    fhi = ub;                                  % highpass
    Wn=[flp/fnq fhi/fnq];                       % butterworth bandpass non-dimensional frequency
    [bfilt_narrow,afilt_narrow] = butter(k,Wn); % construct the filter

        if (iscell(old_tseries))
            xtemp=old_tseries{i};
        else
            xtemp=old_tseries;
        end

        
        xOut_narrow=zeros(nNodes,each_Tmax(i));
        xOut_wide=zeros(nNodes,each_Tmax(i));

        for j=1:nNodes

            Xpre=xtemp(j,:);


            %[f,power1]=fatfourier(Xpif,fs);

            %subplot(4,1,1)
            %plot(Xpif)

            %subplot(4,1,2)
            %plot(f,power1)

            %xlabel('Frequency')
            %ylabel('Power')
            %ESTA TODO HECHO EN VECTOR PORQUE NO ESTA CONFIRMADO SI zscore
            %transforma sobre columnas o vectores de una matriz
            X=detrend(demean(Xpre));
            
            dataOut_narrow =zscore(filtfilt(bfilt_narrow,afilt_narrow,X));
            %dataOut_narrow =zscore(bandpass(X,[lb ub],1/TRsec));
            
                     

            %[f,power2]=fatfourier(dataOut,fs);

            %subplot(4,1,3)
            %plot(f,power2)

            %subplot(4,1,4)
            %plot(dataOut)

            xOut_narrow(j,:)=dataOut_narrow;
            
            
        end
        
        if (iscell(old_tseries))
        new_tseries{i}=xOut_narrow;
        
        else
        new_tseries(:,:,i)=xOut_narrow;
        
        end

    end

end

