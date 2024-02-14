clear all
path2=[ '../../Nonequilibrium/'];
addpath(genpath(path2));
path3=[ '../../Fluctuations_FDT/Classify_FDT'];
addpath(genpath(path3));

N=80;

indexN=1:N; %[1:31 50:80];  %% Cortical areas
Tau=2;
sigma=0.01;

epsFC=0.0004;
epsFCtau=0.0001;
maxC=0.2;

load SC_dbs80HARDIFULL.mat;

Isubdiag = find(tril(ones(N),-1));

% Parameters of the data
TR=1.25;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

C = SC_dbs80HARDI(indexN,indexN);  %% anatomie..DTI tractography
C = C/max(max(C))*maxC;
%%C(C<maxC*0.005)=0;

load psilodep2_extrainfo.mat;

NSUB_P=length(find(drug==2));
NSUB_E=length(find(drug==1));
idx_P=find(drug==2);
idx_E=find(drug==1);

%%
%% Psilo Before
load empirical_psilodep2_before.mat;
load psilodep2_before.mat;
f_diff_B=f_diff;

%% Group
for nsub=1:NSUB_P
    ts=psilodep2_before{idx_P(nsub),1};  % fMRI PB
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCPB(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    COVtauPB(nsub,:,:)=COVtauemp;
end
FCemp=squeeze(mean(FCPB));
COVtauemp=squeeze(mean(COVtauPB));
Cnew=C;
olderror=100000;
for iter=1:5000
    % Linear Hopf FC
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_B,sigma);
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    errorFC(iter)=mean(mean((FCemp-FCsim).^2));
    errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

    if mod(iter,100)<0.1
        errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
        if  (olderror-errornow)/errornow<0.001
            break;
        end
        if  olderror<errornow
            break;
        end
        olderror=errornow;
    end

    for i=1:N  %% learning
        for j=1:N
            if (C(i,j)>0 || j==N-i+1)
                Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                    +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC;
end
CeffgroupPB=Cnew;

%% Individual PB
for nsub=1:NSUB_P
    nsub
    ts=psilodep2_before{idx_P(nsub),1};  % fMRI PB
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCPB(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    Cnew=CeffgroupPB;
    olderror=100000;
    for iter=1:5000
        % Linear Hopf FC
        [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_B,sigma);
        COVtausim=expm((Tau*TR)*A)*COVsimtotal;
        COVtausim=COVtausim(1:N,1:N);
        for i=1:N
            for j=1:N
                sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
            end
        end
        COVtausim=COVtausim.*sigratiosim;
        errorFC(iter)=mean(mean((FCemp-FCsim).^2));
        errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

        if mod(iter,100)<0.1
            errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
            if  (olderror-errornow)/errornow<0.001
                break;
            end
            if  olderror<errornow
                break;
            end
            olderror=errornow;
        end

        for i=1:N  %% learning
            for j=1:N
                if (C(i,j)>0 || j==N-i+1)
                    Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                        +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                end
            end
        end
        Cnew = Cnew/max(max(Cnew))*maxC;
    end
    Ceff=Cnew;
    CeffPB(nsub,:,:)=Ceff;
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Ceff,f_diff_B,sigma);
    fittFC_PB(nsub)=corr2(FCemp(Isubdiag),FCsim(Isubdiag));
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    fittCVtau_PB(nsub)=corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag));
end


%%%%%%%%  
%% Psilo After
%%
load empirical_psilodep2_after.mat;
load psilodep2_after.mat;
f_diff_A=f_diff;

%% Group
for nsub=1:NSUB_P
    ts=psilodep2_after{idx_P(nsub),1};  % fMRI PA
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCPA(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    COVtauPA(nsub,:,:)=COVtauemp;
end
FCemp=squeeze(mean(FCPA));
COVtauemp=squeeze(mean(COVtauPA));
Cnew=C;
olderror=100000;
for iter=1:5000
    % Linear Hopf FC
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_A,sigma);
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    errorFC(iter)=mean(mean((FCemp-FCsim).^2));
    errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

    if mod(iter,100)<0.1
        errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
        if  (olderror-errornow)/errornow<0.001
            break;
        end
        if  olderror<errornow
            break;
        end
        olderror=errornow;
    end

    for i=1:N  %% learning
        for j=1:N
            if (C(i,j)>0 || j==N-i+1)
                Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                    +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC;
end
CeffgroupPA=Cnew;

%% Individual PA
for nsub=1:NSUB_P
    nsub
    ts=psilodep2_after{idx_P(nsub),1};  % fMRI PA
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCPA(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    Cnew=CeffgroupPA;
    olderror=100000;
    for iter=1:5000
        % Linear Hopf FC
        [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_A,sigma);
        COVtausim=expm((Tau*TR)*A)*COVsimtotal;
        COVtausim=COVtausim(1:N,1:N);
        for i=1:N
            for j=1:N
                sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
            end
        end
        COVtausim=COVtausim.*sigratiosim;
        errorFC(iter)=mean(mean((FCemp-FCsim).^2));
        errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

        if mod(iter,100)<0.1
            errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
            if  (olderror-errornow)/errornow<0.001
                break;
            end
            if  olderror<errornow
                break;
            end
            olderror=errornow;
        end

        for i=1:N  %% learning
            for j=1:N
                if (C(i,j)>0 || j==N-i+1)
                    Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                        +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                end
            end
        end
        Cnew = Cnew/max(max(Cnew))*maxC;
    end
    Ceff=Cnew;
    CeffPA(nsub,:,:)=Ceff;
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Ceff,f_diff_A,sigma);
    fittFC_PA(nsub)=corr2(FCemp(Isubdiag),FCsim(Isubdiag));
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    fittCVtau_PA(nsub)=corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag));
end


%%%%%%%%
%% Escitalopram before

%% Group
for nsub=1:NSUB_E
    ts=psilodep2_before{idx_E(nsub),1};  % fMRI EB
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCEB(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    COVtauEB(nsub,:,:)=COVtauemp;
end
FCemp=squeeze(mean(FCEB));
COVtauemp=squeeze(mean(COVtauEB));
Cnew=C;
olderror=100000;
for iter=1:5000
    % Linear Hopf FC
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_B,sigma);
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    errorFC(iter)=mean(mean((FCemp-FCsim).^2));
    errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

    if mod(iter,100)<0.1
        errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
        if  (olderror-errornow)/errornow<0.001
            break;
        end
        if  olderror<errornow
            break;
        end
        olderror=errornow;
    end

    for i=1:N  %% learning
        for j=1:N
            if (C(i,j)>0 || j==N-i+1)
                Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                    +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC;
end
CeffgroupEB=Cnew;

%% Individual EB
for nsub=1:NSUB_E
    nsub
    ts=psilodep2_before{idx_E(nsub),1};  % fMRI EB
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCEB(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    Cnew=CeffgroupEB;
    olderror=100000;
    for iter=1:5000
        % Linear Hopf FC
        [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_B,sigma);
        COVtausim=expm((Tau*TR)*A)*COVsimtotal;
        COVtausim=COVtausim(1:N,1:N);
        for i=1:N
            for j=1:N
                sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
            end
        end
        COVtausim=COVtausim.*sigratiosim;
        errorFC(iter)=mean(mean((FCemp-FCsim).^2));
        errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

        if mod(iter,100)<0.1
            errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
            if  (olderror-errornow)/errornow<0.001
                break;
            end
            if  olderror<errornow
                break;
            end
            olderror=errornow;
        end

        for i=1:N  %% learning
            for j=1:N
                if (C(i,j)>0 || j==N-i+1)
                    Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                        +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                end
            end
        end
        Cnew = Cnew/max(max(Cnew))*maxC;
    end
    Ceff=Cnew;
    CeffEB(nsub,:,:)=Ceff;
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Ceff,f_diff_B,sigma);
    fittFC_EB(nsub)=corr2(FCemp(Isubdiag),FCsim(Isubdiag));
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    fittCVtau_EB(nsub)=corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag));
end


%%%%%%%%  
%% Escitalopram After
%%

%% Group
for nsub=1:NSUB_E
    ts=psilodep2_after{idx_E(nsub),1};  % fMRI EA
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCEA(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    COVtauEA(nsub,:,:)=COVtauemp;
end
FCemp=squeeze(mean(FCEA));
COVtauemp=squeeze(mean(COVtauEA));
Cnew=C;
olderror=100000;
for iter=1:5000
    % Linear Hopf FC
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_A,sigma);
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    errorFC(iter)=mean(mean((FCemp-FCsim).^2));
    errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

    if mod(iter,100)<0.1
        errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
        if  (olderror-errornow)/errornow<0.001
            break;
        end
        if  olderror<errornow
            break;
        end
        olderror=errornow;
    end

    for i=1:N  %% learning
        for j=1:N
            if (C(i,j)>0 || j==N-i+1)
                Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                    +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC;
end
CeffgroupEA=Cnew;

%% Individual EA
for nsub=1:NSUB_E
    nsub
    ts=psilodep2_after{idx_E(nsub),1};  % fMRI EA
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCEA(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    Cnew=CeffgroupEA;
    olderror=100000;
    for iter=1:5000
        % Linear Hopf FC
        [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_A,sigma);
        COVtausim=expm((Tau*TR)*A)*COVsimtotal;
        COVtausim=COVtausim(1:N,1:N);
        for i=1:N
            for j=1:N
                sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
            end
        end
        COVtausim=COVtausim.*sigratiosim;
        errorFC(iter)=mean(mean((FCemp-FCsim).^2));
        errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

        if mod(iter,100)<0.1
            errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
            if  (olderror-errornow)/errornow<0.001
                break;
            end
            if  olderror<errornow
                break;
            end
            olderror=errornow;
        end

        for i=1:N  %% learning
            for j=1:N
                if (C(i,j)>0 || j==N-i+1)
                    Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                        +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                end
            end
        end
        Cnew = Cnew/max(max(Cnew))*maxC;
    end
    Ceff=Cnew;
    CeffEA(nsub,:,:)=Ceff;
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Ceff,f_diff_A,sigma);
    fittFC_EA(nsub)=corr2(FCemp(Isubdiag),FCsim(Isubdiag));
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    fittCVtau_EA(nsub)=corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Hierarchy Ceff

%% PB
for nsub=1:NSUB_P
    Ceff=squeeze(CeffPB(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevelsPB(nsub,:)=gamma';
    H=(meshgrid(gamma)-meshgrid(gamma)'-1).^2;
    F0=sum(sum((A.*H)))/sum(sum(A));
    trophiccoherencePB(nsub)=1-F0;
end
%% PA
for nsub=1:NSUB_P
    Ceff=squeeze(CeffPA(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevelsPA(nsub,:)=gamma';
    H=(meshgrid(gamma)-meshgrid(gamma)'-1).^2;
    F0=sum(sum((A.*H)))/sum(sum(A));
    trophiccoherencePA(nsub)=1-F0;
end

%% EB
for nsub=1:NSUB_E
    Ceff=squeeze(CeffEB(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevelsEB(nsub,:)=gamma';
    H=(meshgrid(gamma)-meshgrid(gamma)'-1).^2;
    F0=sum(sum((A.*H)))/sum(sum(A));
    trophiccoherenceEB(nsub)=1-F0;
end

%% EA
for nsub=1:NSUB_E
    Ceff=squeeze(CeffEA(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevelsEA(nsub,:)=gamma';
    H=(meshgrid(gamma)-meshgrid(gamma)'-1).^2;
    F0=sum(sum((A.*H)))/sum(sum(A));
    trophiccoherenceEA(nsub)=1-F0;
end

figure(1);
subplot(1,3,1)
boxplot([trophiccoherencePB' trophiccoherencePA']);
a=trophiccoherencePB;
b=trophiccoherencePA;
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'signrank');
min(stats.pvals)

subplot(1,3,2)
boxplot([trophiccoherenceEB' trophiccoherenceEA']);
a=trophiccoherenceEB;
b=trophiccoherenceEA;
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'signrank');
min(stats.pvals)

subplot(1,3,3)
boxplot([trophiccoherencePA trophiccoherenceEA],[ones(length(trophiccoherencePA),1)', 2*ones(length(trophiccoherenceEA),1)']);
a=trophiccoherencePA;
b=trophiccoherenceEA;
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
min(stats.pvals)

a=trophiccoherencePB;
b=trophiccoherenceEB;
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
min(stats.pvals)

figure(2);
subplot(2,2,1);
boxplot([fittFC_PB' fittFC_PA']);
subplot(2,2,2);
boxplot([fittFC_EB' fittFC_EA']);
subplot(2,2,3);
boxplot([fittCVtau_PB' fittCVtau_PA']);
subplot(2,2,4);
boxplot([fittCVtau_EB' fittCVtau_EA']);

save results_Ceff_psilodep2.mat CeffPB CeffPA CeffEB CeffEA ...
     trophiccoherenceEB trophiccoherenceEA trophiccoherencePB trophiccoherencePA ...
     hierarchicallevelsPB hierarchicallevelsPA hierarchicallevelsEB hierarchicallevelsEA;
