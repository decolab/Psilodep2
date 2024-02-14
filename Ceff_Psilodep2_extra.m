%% correlation change in hierarchy vs BDIscore (p value, all or resp? siehe paper Nat Med.)
%% Select top hierachicallevels change per condition SSND, mean (Ceffin-Ceffout), class...

clear all
path2=[ '../../Nonequilibrium/'];
addpath(genpath(path2));
path3=[ '../../Fluctuations_FDT/Classify_FDT'];
addpath(genpath(path3));

load results_Ceff_psilodep2.mat;
load psilodep2_extrainfo.mat;

%% group change of Hierarchy
N=80;

%%

idx_P=find(drug==2);
idx_E=find(drug==1);

idx_P_Nresp=find(BDIresponse(idx_P)==0);
idx_P_resp=find(BDIresponse(idx_P)==1);
idx_E_Nresp=find(BDIresponse(idx_E)==0);
idx_E_resp=find(BDIresponse(idx_E)==1);

%%  P A-B changes Nr vs R   (CASE=1)

CASE=1;

Ceff1=CeffPB;
Ceff2=CeffPA;
NSUB1=1:length(idx_P);
NSUB2=1:length(idx_P);

sub=1;
for nsub=idx_P_Nresp'
    Ceff=squeeze(Ceff1(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub1_NR_P(sub,:)=gamma';
    sub=sub+1;
end
sub=1;
for nsub=idx_P_Nresp'
    Ceff=squeeze(Ceff2(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub2_NR_P(sub,:)=gamma';
    sub=sub+1;
end
 changedhierarchyNR=abs(hierarchicallevels_sub2_NR_P-hierarchicallevels_sub1_NR_P);

 sub=1;
for nsub=idx_P_resp'
    Ceff=squeeze(Ceff1(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub1_R_P(sub,:)=gamma';
    sub=sub+1;
end
sub=1;
for nsub=idx_P_resp'
    Ceff=squeeze(Ceff2(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub2_R_P(sub,:)=gamma';
    sub=sub+1;
end
 changedhierarchyR=abs(hierarchicallevels_sub2_R_P-hierarchicallevels_sub1_R_P);

for n=1:N
    SSND(n)=abs(mean(changedhierarchyR(:,n))-mean(changedhierarchyNR(:,n))) ...
        /sqrt(var(changedhierarchyR(:,n))+var(changedhierarchyNR(:,n)));
end
[soSSND index]=sort(SSND,'descend');

for n=1:10
    clear trophiccoherence1NR_P trophiccoherence2NR_P trophiccoherence1R_P trophiccoherence2R_P;
    sub=1;
    for nsub=idx_P_Nresp'
        gamma=hierarchicallevels_sub1_NR_P(sub,index(1:n));
        trophiccoherence1NR_P(sub,:)=gamma;
        sub=sub+1;
    end
    sub=1;
    for nsub=idx_P_Nresp'
        gamma=hierarchicallevels_sub2_NR_P(sub,index(1:n));
        trophiccoherence2NR_P(sub,:)=gamma;
        sub=sub+1;
    end
    sub=1;
    for nsub=idx_P_resp'
        gamma=hierarchicallevels_sub1_R_P(sub,index(1:n));
        trophiccoherence1R_P(sub,:)=gamma;
        sub=sub+1;
    end
    sub=1;
    for nsub=idx_P_resp'
        gamma=hierarchicallevels_sub2_R_P(sub,index(1:n));
        trophiccoherence2R_P(sub,:)=gamma;
        sub=sub+1;
    end
    a=mean(abs(trophiccoherence2NR_P-trophiccoherence1NR_P),2)';
    b=mean(abs(trophiccoherence2R_P-trophiccoherence1R_P),2)';
    stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
    pp(n)=min(stats.pvals);
end
[aux minn]=min(pp);
n=minn;
clear trophiccoherence1NR_P trophiccoherence2NR_P trophiccoherence1R_P trophiccoherence2R_P;
sub=1;
for nsub=idx_P_Nresp'
    gamma=hierarchicallevels_sub1_NR_P(sub,index(1:n));
    trophiccoherence1NR_P(sub,:)=gamma;
    sub=sub+1;
end
sub=1;
for nsub=idx_P_Nresp'
    gamma=hierarchicallevels_sub2_NR_P(sub,index(1:n));
    trophiccoherence2NR_P(sub,:)=gamma;
    sub=sub+1;
end
sub=1;
for nsub=idx_P_resp'
    gamma=hierarchicallevels_sub1_R_P(sub,index(1:n));
    trophiccoherence1R_P(sub,:)=gamma;
    sub=sub+1;
end
sub=1;
for nsub=idx_P_resp'
    gamma=hierarchicallevels_sub2_R_P(sub,index(1:n));
    trophiccoherence2R_P(sub,:)=gamma;
    sub=sub+1;
end
a=mean(abs(trophiccoherence2NR_P-trophiccoherence1NR_P),2)';
b=mean(abs(trophiccoherence2R_P-trophiccoherence1R_P),2)';
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
pvalue(CASE)=min(stats.pvals)

figure(CASE)
boxplot([a b],[zeros(1,length(a)) ones(1,length(b))]);

changetrP(idx_P_Nresp')=a;
changetrP(idx_P_resp')=b;
[cc ppv]=corrcoef(changetrP',BDIscore(idx_P));
cc(1,2)
ppv(1,2)
index(1:n)

%%  E A-B changes Nr vs R   (CASE=1)

CASE=2;

Ceff1=CeffEB;
Ceff2=CeffEA;
NSUB1=1:length(idx_E);
NSUB2=1:length(idx_E);

sub=1;
for nsub=idx_E_Nresp'
    Ceff=squeeze(Ceff1(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub1_NR_E(sub,:)=gamma';
    sub=sub+1;
end
sub=1;
for nsub=idx_E_Nresp'
    Ceff=squeeze(Ceff2(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub2_NR_E(sub,:)=gamma';
    sub=sub+1;
end
 changedhierarchyNR=abs(hierarchicallevels_sub2_NR_E-hierarchicallevels_sub1_NR_E);

 sub=1;
for nsub=idx_E_resp'
    Ceff=squeeze(Ceff1(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub1_R_E(sub,:)=gamma';
    sub=sub+1;
end
sub=1;
for nsub=idx_E_resp'
    Ceff=squeeze(Ceff2(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub2_R_E(sub,:)=gamma';
    sub=sub+1;
end
 changedhierarchyR=abs(hierarchicallevels_sub2_R_E-hierarchicallevels_sub1_R_E);

for n=1:N
    SSND(n)=abs(mean(changedhierarchyR(:,n))-mean(changedhierarchyNR(:,n))) ...
        /sqrt(var(changedhierarchyR(:,n))+var(changedhierarchyNR(:,n)));
end
[soSSND index]=sort(SSND,'descend');

for n=1:10
    clear trophiccoherence1NR_E trophiccoherence2NR_E trophiccoherence1R_E trophiccoherence2R_E;
    sub=1;
    for nsub=idx_E_Nresp'
        gamma=hierarchicallevels_sub1_NR_E(sub,index(1:n));
        trophiccoherence1NR_E(sub,:)=gamma;
        sub=sub+1;
    end
    sub=1;
    for nsub=idx_E_Nresp'
        gamma=hierarchicallevels_sub2_NR_E(sub,index(1:n));
        trophiccoherence2NR_E(sub,:)=gamma;
        sub=sub+1;
    end
    sub=1;
    for nsub=idx_E_resp'
        gamma=hierarchicallevels_sub1_R_E(sub,index(1:n));
        trophiccoherence1R_E(sub,:)=gamma;
        sub=sub+1;
    end
    sub=1;
    for nsub=idx_E_resp'
        gamma=hierarchicallevels_sub2_R_E(sub,index(1:n));
        trophiccoherence2R_E(sub,:)=gamma;
        sub=sub+1;
    end
    a=mean(abs(trophiccoherence2NR_E-trophiccoherence1NR_E),2)';
    b=mean(abs(trophiccoherence2R_E-trophiccoherence1R_E),2)';
    stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
    pp(n)=min(stats.pvals);
end
[aux minn]=min(pp);
n=minn;
clear trophiccoherence1NR_E trophiccoherence2NR_E trophiccoherence1R_E trophiccoherence2R_E;
sub=1;
for nsub=idx_E_Nresp'
    gamma=hierarchicallevels_sub1_NR_E(sub,index(1:n));
    trophiccoherence1NR_E(sub,:)=gamma;
    sub=sub+1;
end
sub=1;
for nsub=idx_E_Nresp'
    gamma=hierarchicallevels_sub2_NR_E(sub,index(1:n));
    trophiccoherence2NR_E(sub,:)=gamma;
    sub=sub+1;
end
sub=1;
for nsub=idx_E_resp'
    gamma=hierarchicallevels_sub1_R_E(sub,index(1:n));
    trophiccoherence1R_E(sub,:)=gamma;
    sub=sub+1;
end
sub=1;
for nsub=idx_E_resp'
    gamma=hierarchicallevels_sub2_R_E(sub,index(1:n));
    trophiccoherence2R_E(sub,:)=gamma;
    sub=sub+1;
end
a=mean(abs(trophiccoherence2NR_E-trophiccoherence1NR_E),2)';
b=mean(abs(trophiccoherence2R_E-trophiccoherence1R_E),2)';
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
pvalue(CASE)=min(stats.pvals)

figure(CASE)
boxplot([a b],[zeros(1,length(a)) ones(1,length(b))]);

changetrE(idx_E_Nresp')=a;
changetrE(idx_E_resp')=b;
[cc ppv]=corrcoef(changetrE',BDIscore(idx_E));
cc(1,2)
ppv(1,2)
index(1:n)

%%  P Before R vs NR

CASE=3;
Ceff1=CeffPB;
Ceff2=CeffPB;

sub=1;
for nsub=idx_P_Nresp'
    Ceff=squeeze(Ceff1(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub_NR_P(sub,:)=gamma';
    sub=sub+1;
end

sub=1;
for nsub=idx_P_resp'
    Ceff=squeeze(Ceff2(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub_R_P(sub,:)=gamma';
    sub=sub+1;
end

for n=1:N
    SSND(n)=abs(mean(hierarchicallevels_sub_R_P(:,n))-mean(hierarchicallevels_sub_NR_P(:,n))) ...
        /sqrt(var(hierarchicallevels_sub_R_P(:,n))+var(hierarchicallevels_sub_NR_P(:,n)));
end
[soSSND index]=sort(SSND,'descend');

for n=1:10
    clear trophiccoherenceNR_P trophiccoherenceR_P;
    sub=1;
    for nsub=idx_P_Nresp'
        gamma=hierarchicallevels_sub_NR_P(sub,index(1:n));
        trophiccoherenceNR_P(sub,:)=gamma;
        sub=sub+1;
    end
    sub=1;
    for nsub=idx_P_resp'
        gamma=hierarchicallevels_sub_R_P(sub,index(1:n));
        trophiccoherenceR_P(sub,:)=gamma;
        sub=sub+1;
    end
    a=mean(abs(trophiccoherenceNR_P),2)';
    b=mean(abs(trophiccoherenceR_P),2)';
    stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
    pp(n)=min(stats.pvals);
end
[aux minn]=min(pp);
n=minn;
clear trophiccoherenceNR_P trophiccoherenceR_P;
sub=1;
for nsub=idx_P_Nresp'
    gamma=hierarchicallevels_sub_NR_P(sub,index(1:n));
    trophiccoherenceNR_P(sub,:)=gamma;
    sub=sub+1;
end
sub=1;
for nsub=idx_P_resp'
    gamma=hierarchicallevels_sub_R_P(sub,index(1:n));
    trophiccoherenceR_P(sub,:)=gamma;
    sub=sub+1;
end
a=mean(abs(trophiccoherenceNR_P),2)';
b=mean(abs(trophiccoherenceR_P),2)';
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
pvalue(CASE)=min(stats.pvals)
figure(CASE)
boxplot([a b],[zeros(1,length(a)) ones(1,length(b))]);

baseP(idx_P_Nresp')=a;
baseP(idx_P_resp')=b;
[cc ppv]=corrcoef(baseP',BDIscore(idx_P));
cc(1,2)
ppv(1,2)
index(1:n)

%%  E Before R vs NR

CASE=4;
Ceff1=CeffEB;
Ceff2=CeffEB;

sub=1;
for nsub=idx_E_Nresp'
    Ceff=squeeze(Ceff1(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub_NR_E(sub,:)=gamma';
    sub=sub+1;
end

sub=1;
for nsub=idx_E_resp'
    Ceff=squeeze(Ceff2(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub_R_E(sub,:)=gamma';
    sub=sub+1;
end

for n=1:N
    SSND(n)=abs(mean(hierarchicallevels_sub_R_E(:,n))-mean(hierarchicallevels_sub_NR_E(:,n))) ...
        /sqrt(var(hierarchicallevels_sub_R_E(:,n))+var(hierarchicallevels_sub_NR_E(:,n)));
end
[soSSND index]=sort(SSND,'descend');

for n=1:10
    clear trophiccoherenceNR_E trophiccoherenceR_E;
    sub=1;
    for nsub=idx_E_Nresp'
        gamma=hierarchicallevels_sub_NR_E(sub,index(1:n));
        trophiccoherenceNR_E(sub,:)=gamma;
        sub=sub+1;
    end
    sub=1;
    for nsub=idx_E_resp'
        gamma=hierarchicallevels_sub_R_E(sub,index(1:n));
        trophiccoherenceR_E(sub,:)=gamma;
        sub=sub+1;
    end
    a=mean(abs(trophiccoherenceNR_E),2)';
    b=mean(abs(trophiccoherenceR_E),2)';
    stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
    pp(n)=min(stats.pvals);
end
[aux minn]=min(pp);
n=minn;
clear trophiccoherenceNR_E trophiccoherenceR_E;
sub=1;
for nsub=idx_E_Nresp'
    gamma=hierarchicallevels_sub_NR_E(sub,index(1:n));
    trophiccoherenceNR_E(sub,:)=gamma;
    sub=sub+1;
end
sub=1;
for nsub=idx_E_resp'
    gamma=hierarchicallevels_sub_R_E(sub,index(1:n));
    trophiccoherenceR_E(sub,:)=gamma;
    sub=sub+1;
end
a=mean(abs(trophiccoherenceNR_E),2)';
b=mean(abs(trophiccoherenceR_E),2)';
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
pvalue(CASE)=min(stats.pvals)
figure(CASE)
boxplot([a b],[zeros(1,length(a)) ones(1,length(b))]);

baseE(idx_E_Nresp')=a;
baseE(idx_E_resp')=b;
[cc ppv]=corrcoef(baseE',BDIscore(idx_E));
cc(1,2)
ppv(1,2)
index(1:n)
