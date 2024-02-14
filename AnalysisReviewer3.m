clear all
path2=[ '../../Nonequilibrium/'];
addpath(genpath(path2));
path3=[ '../../Fluctuations_FDT/Classify_FDT'];
addpath(genpath(path3));

load results_Ceff_psilodep2.mat;
load psilodep2_extrainfo.mat;
load BDI_baseline.mat;
load psilodep2_gender.mat;

idx_P=find(drug==2);
idx_E=find(drug==1);

idx_P_Nresp=find(BDIresponse(idx_P)==0);
idx_P_resp=find(BDIresponse(idx_P)==1);
idx_E_Nresp=find(BDIresponse(idx_E)==0);
idx_E_resp=find(BDIresponse(idx_E)==1);

idx_P_male=find(gender(idx_P)==1);
idx_P_female=find(gender(idx_P)==0);
idx_E_male=find(gender(idx_E)==1);
idx_E_female=find(gender(idx_E)==0);

%% Gender
figure(1)
subplot(2,2,1)
a=trophiccoherencePB(idx_P_female);
b=trophiccoherencePA(idx_P_female);
boxplot([a' b']);
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ttest2');
min(stats.pvals)
subplot(2,2,2)
a=trophiccoherencePB(idx_P_male);
b=trophiccoherencePA(idx_P_male);
boxplot([a' b']);
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ttest2');
min(stats.pvals)
subplot(2,2,3)
a=trophiccoherenceEB(idx_E_female);
b=trophiccoherenceEA(idx_E_female);
boxplot([a' b']);
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ttest2');
min(stats.pvals)
subplot(2,2,4)
a=trophiccoherenceEB(idx_E_male);
b=trophiccoherenceEA(idx_E_male);
boxplot([a' b']);
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ttest2');
min(stats.pvals)

%%

aa=BDI_baseline>23;
bb=BDIscore<-10;
cc=aa+bb;
dd=cc>1;
idx_P_similarBDI=find(dd(idx_P)==1);
idx_E_similarBDI=find(dd(idx_E)==1);

%%% similar Baseline and scorechange

figure(1);

a=trophiccoherencePB(idx_P_similarBDI);
b=trophiccoherencePA(idx_P_similarBDI);
boxplot([a' b']);
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ttest2');
min(stats.pvals)

figure(2);

a=trophiccoherenceEB(idx_E_similarBDI);
b=trophiccoherenceEA(idx_E_similarBDI);
boxplot([a' b']);
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ttest2');
min(stats.pvals)
%%

figure(1);

a=trophiccoherencePB(idx_P_resp);
b=trophiccoherencePA(idx_P_resp);
boxplot([a' b']);
signrank(a,b)


figure(2);

a=trophiccoherenceEB(idx_E_resp);
b=trophiccoherenceEA(idx_E_resp);
boxplot([a' b']);
signrank(a,b)

figure(3);

a2=trophiccoherenceEB(idx_E_Nresp);
b2=trophiccoherenceEA(idx_E_Nresp);
boxplot([a2' b2']);
signrank(a2,b2)

figure(4)
boxplot([(b-a)' ;(b2-a2)'],[ones(length((b-a)'),1)', 2*ones(length((b2-a2)'),1)']);
ranksum(b-a,b2-a2)


%%% FCtau

clear all
path2=[ '../../Nonequilibrium/'];
addpath(genpath(path2));
path3=[ '../../Fluctuations_FDT/Classify_FDT'];
addpath(genpath(path3));

N=80;

indexN=1:N; %[1:31 50:80];  %% Cortical areas

% Parameters of the data
TR=1.25;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

load psilodep2_extrainfo.mat;

NSUB_P=length(find(drug==2));
NSUB_E=length(find(drug==1));
idx_P=find(drug==2);
idx_E=find(drug==1);

%%
%% Psilo Before
load psilodep2_before.mat;

Tau=2;
NSUB=42;

%% Group
for nsub=1:NSUB
    ts=psilodep2_before{nsub,1};  % fMRI PB
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCe=corrcoef(ts2');
    FCemp=corr(ts2(:,1:Tm-Tau)',ts2(:,1+Tau:Tm)');
    FCemp=FCemp-eye(N).*diag(FCemp);
    vec=FCemp-FCemp';
    FCTauB(nsub)=mean(mean(abs(vec)));
    FCB(nsub)=mean(mean(abs(FCe)));
end

load psilodep2_after.mat;
for nsub=1:NSUB
    ts=psilodep2_after{nsub,1};  % fMRI PB
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCe=corrcoef(ts2');
    FCemp=corr(ts2(:,1:Tm-Tau)',ts2(:,1+Tau:Tm)');
    FCemp=FCemp-eye(N).*diag(FCemp);
    vec=FCemp-FCemp';
    FCTauA(nsub)=mean(mean(abs(vec)));
    FCA(nsub)=mean(mean(abs(FCe)));
end

figure(1);
subplot(2,2,1)
boxplot([FCB(idx_P)' FCA(idx_P)']);
a=FCB(idx_P);
b=FCA(idx_P);
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'signrank');
min(stats.pvals)

subplot(2,2,2)
boxplot([FCB(idx_E)' FCA(idx_E)']);
a=FCB(idx_E);
b=FCA(idx_E);
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'signrank');
min(stats.pvals)

subplot(2,2,3)
boxplot([FCTauB(idx_P)' FCTauA(idx_P)']);
a=FCTauB(idx_P);
b=FCTauA(idx_P);
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'signrank');
min(stats.pvals)

subplot(2,2,4)
boxplot([FCTauB(idx_E)' FCTauA(idx_E)']);
a=FCTauB(idx_E);
b=FCTauA(idx_E);
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'signrank');
min(stats.pvals)


%%%  Analysis per Node

hierarchicallevelsPBm=squeeze(mean(hierarchicallevelsPB));
hierarchicallevelsPAm=squeeze(mean(hierarchicallevelsPA));
hierarchicallevelsEBm=squeeze(mean(hierarchicallevelsEB));
hierarchicallevelsEAm=squeeze(mean(hierarchicallevelsEA));

diffP=hierarchicallevelsPBm-hierarchicallevelsPAm;
diffE=hierarchicallevelsEBm-hierarchicallevelsEAm;

boxplot([diffP' diffE']);
a=diffP;
b=diffE;
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
min(stats.pvals)

for i=1:N
    i
    dif1=hierarchicallevelsPB-hierarchicallevelsPA;
    dif2=hierarchicallevelsEB-hierarchicallevelsEA;
    a=dif1(:,i)';
    b=dif2(:,i)';
    stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
    pp(i)=min(stats.pvals);
end
nodos=FDR_benjHoch(pp,0.05)
