function [Clong lambda numexcSC valexcSC]=TTB_SC_analysis_longrange_EDR_DTI(Cfg,CoG,SC);


linfunc = @(A, x)(A(1)*x+A(2));
options=optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');

NPARCELLS=Cfg.nNodes;
NR=Cfg.NR;
NRini=Cfg.NRini;
NRfin=Cfg.NRfin;

% so far harcoded
NSTD=1;
DistRange=0;

Isubdiag = find(tril(ones(NPARCELLS),-1));


C=SC;
C=C/max(max(C));


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

%% SC
numind3=zeros(1,NR);
sc_3=zeros(1,NR);
sc_density=cell(1,NR);
sc_density_i=cell(1,NR);
sc_density_j=cell(1,NR);
for i=1:NPARCELLS
    for j=1:NPARCELLS
        r=rr(i,j);
        index=floor(r/delta)+1;
        if index==NR+1
            index=NR;
        end
        sc_density{index}=[sc_density{index} C(i,j)];
        sc_density_i{index}=[sc_density_i{index} i];
        sc_density_j{index}=[sc_density_j{index} j];
        sc_3(index)=sc_3(index)+C(i,j);
        numind3(index)=numind3(index)+1;
    end
end
sctotal=sc_3./numind3;

%%%%

for k=1:NR
    xcoor(k)=xrange(k);
    ycoor(k)=log(sctotal(k));
    ycoor2(k)=sctotal(k);
end
A0=[-1 1];
Afit = lsqcurvefit(linfunc,A0,xcoor(20:50),ycoor(20:50),[-4 -10],[4 10],options);
lambda=Afit(1);
yl=Afit(1)*xcoor+Afit(2);

expfunc = @(A, x)(A(1)*exp(-A(2)*x));
options=optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');
A0=[0.15 0.18];
Afit = lsqcurvefit(expfunc,A0,xcoor(25:end),ycoor2(25:end),[-100 -100],[100 100],options);
yl=Afit(1)*exp(-Afit(2)*xcoor);
lambda=Afit(2);


%%%%

Clong=zeros(NPARCELLS,NPARCELLS);
numexcSC=zeros(NR,1);
valexcSC=zeros(NR,1);

for i=NRini:NRfin
    mv=nanmean(sc_density{i});
    st=nanstd(sc_density{i});
    numexcSC(i)=length(find(sc_density{i}>mv+NSTD*st))/length(sc_density{i});
    valexcSC(i)=mean(sc_density{i}(find(sc_density{i}>mv+NSTD*st)))/mean(sc_density{i});
    for ind=find(sc_density{i}>mv+NSTD*st)
        if rr(sc_density_i{i}(ind),sc_density_j{i}(ind))>DistRange
            Clong(sc_density_i{i}(ind),sc_density_j{i}(ind))=sc_density{i}(ind);
        end
    end
end

Clong=Clong/Afit(1);
SC=C/Afit(1);


%save SClongrange.mat SC Clong lambda numexcSC valexcSC;
