clear all


load results_Ceff_psilodep2.mat;
load psilodep2_extrainfo.mat;

kfold=1000;


%% group change of Hierarchy
N=80;

%%  P A-B changes Nr vs R   (CASE=1)
%%  E A-B changes Nr vs R   (CASE=2)
%%  P B - Nr vs R           (CASE=3)
%%  E B - Nr vs R           (CASE=4)
%%  P A - Nr vs R           (CASE=5)
%%  E A - Nr vs R           (CASE=6)
%%  A R - E vs P            (CASE=7)
%%  P  -  A vs B            (CASE=8)
%%  E  -  A vs B            (CASE=9)

idx_P=find(drug==2);
idx_E=find(drug==1);

idx_P_Nresp=find(BDIresponse(idx_P)==0);
idx_P_resp=find(BDIresponse(idx_P)==1);
idx_E_Nresp=find(BDIresponse(idx_E)==0);
idx_E_resp=find(BDIresponse(idx_E)==1);


for CASE=[1 2 4 7 8 9]
    if CASE==1
        Ceff1=CeffPB;
        Ceff2=CeffPA;
        NSUB1=length(idx_P_Nresp);
        NSUB2=length(idx_P_resp);

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

        renderingdiff(CASE,:)=mean(hierarchicallevels_sub2_R_P-hierarchicallevels_sub1_R_P)-mean(hierarchicallevels_sub2_NR_P-hierarchicallevels_sub1_NR_P);

        for n=1:N
            clear trophiccoherence1NR_P trophiccoherence2NR_P trophiccoherence1R_P trophiccoherence2R_P;
            sub=1;
            for nsub=idx_P_Nresp'
                gamma=hierarchicallevels_sub1_NR_P(sub,n);
                trophiccoherence1NR_P(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_Nresp'
                gamma=hierarchicallevels_sub2_NR_P(sub,n);
                trophiccoherence2NR_P(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_resp'
                gamma=hierarchicallevels_sub1_R_P(sub,n);
                trophiccoherence1R_P(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_resp'
                gamma=hierarchicallevels_sub2_R_P(sub,n);
                trophiccoherence2R_P(sub,:)=gamma;
                sub=sub+1;
            end
            a=mean((trophiccoherence2NR_P-trophiccoherence1NR_P).^2,2)';
            b=mean((trophiccoherence2R_P-trophiccoherence1R_P).^2,2)';
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ranksum');
            pp(n)=min(stats.pvals);
        end
        [aa1 indx] = sort(pp,'ascend');
        nTot = 6;
        nsig=indx(1:nTot);

        if length(nsig)>0
            clear trophiccoherence1NR_P trophiccoherence2NR_P trophiccoherence1R_P trophiccoherence2R_P;
            sub=1;
            for nsub=idx_P_Nresp'
                gamma=hierarchicallevels_sub1_NR_P(sub,nsig);
                trophiccoherence1NR_P(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_Nresp'
                gamma=hierarchicallevels_sub2_NR_P(sub,nsig);
                trophiccoherence2NR_P(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_resp'
                gamma=hierarchicallevels_sub1_R_P(sub,nsig);
                trophiccoherence1R_P(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_resp'
                gamma=hierarchicallevels_sub2_R_P(sub,nsig);
                trophiccoherence2R_P(sub,:)=gamma;
                sub=sub+1;
            end
            a=trophiccoherence2NR_P-trophiccoherence1NR_P;
            b=trophiccoherence2R_P-trophiccoherence1R_P;
            xxdata=zscore([a;b]);


            %%
        end
    elseif CASE==2
        Ceff1=CeffEB;
        Ceff2=CeffEA;
        NSUB1=length(idx_E_Nresp);
        NSUB2=length(idx_E_resp);

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

        renderingdiff(CASE,:)=mean(hierarchicallevels_sub2_R_E-hierarchicallevels_sub1_R_E)-mean(hierarchicallevels_sub2_NR_E-hierarchicallevels_sub1_NR_E);

        for n=1:N
            clear trophiccoherence1NR_E trophiccoherence2NR_E trophiccoherence1R_E trophiccoherence2R_E;
            sub=1;
            for nsub=idx_E_Nresp'
                gamma=hierarchicallevels_sub1_NR_E(sub,n);
                trophiccoherence1NR_E(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_E_Nresp'
                gamma=hierarchicallevels_sub2_NR_E(sub,n);
                trophiccoherence2NR_E(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_E_resp'
                gamma=hierarchicallevels_sub1_R_E(sub,n);
                trophiccoherence1R_E(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_E_resp'
                gamma=hierarchicallevels_sub2_R_E(sub,n);
                trophiccoherence2R_E(sub,:)=gamma;
                sub=sub+1;
            end
            a=mean((trophiccoherence2NR_E-trophiccoherence1NR_E).^2,2)';
            b=mean((trophiccoherence2R_E-trophiccoherence1R_E).^2,2)';
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ranksum');
            pp(n)=min(stats.pvals);
        end
        [aa1 indx] = sort(pp,'ascend');
        nTot = 5;
        nsig=indx(1:nTot);
        
        if length(nsig)>0
            clear trophiccoherence1NR_E trophiccoherence2NR_E trophiccoherence1R_E trophiccoherence2R_E;
            sub=1;
            for nsub=idx_E_Nresp'
                gamma=hierarchicallevels_sub1_NR_E(sub,nsig);
                trophiccoherence1NR_E(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_E_Nresp'
                gamma=hierarchicallevels_sub2_NR_E(sub,nsig);
                trophiccoherence2NR_E(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_E_resp'
                gamma=hierarchicallevels_sub1_R_E(sub,nsig);
                trophiccoherence1R_E(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_E_resp'
                gamma=hierarchicallevels_sub2_R_E(sub,nsig);
                trophiccoherence2R_E(sub,:)=gamma;
                sub=sub+1;
            end
            a=trophiccoherence2NR_E-trophiccoherence1NR_E;
            b=trophiccoherence2R_E-trophiccoherence1R_E;
            xxdata=zscore([a;b]);


        end
    elseif CASE==3
        Ceff1=CeffPB;
        Ceff2=CeffPB;
        NSUB1=length(idx_P_Nresp);
        NSUB2=length(idx_P_resp);
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

        renderingdiff(CASE,:)=mean(hierarchicallevels_sub_R_P)-mean(hierarchicallevels_sub_NR_P);

        for n=1:N
            clear trophiccoherenceNR_P trophiccoherenceR_P;
            sub=1;
            for nsub=idx_P_Nresp'
                gamma=hierarchicallevels_sub_NR_P(sub,n);
                trophiccoherenceNR_P(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_resp'
                gamma=hierarchicallevels_sub_R_P(sub,n);
                trophiccoherenceR_P(sub,:)=gamma;
                sub=sub+1;
            end
            a=mean((trophiccoherenceNR_P).^2,2)';
            b=mean((trophiccoherenceR_P).^2,2)';
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ranksum');
            pp(n)=min(stats.pvals);
        end

        nsig=FDR_benjHoch(pp,0.05);
        if length(nsig)>0

            clear trophiccoherenceNR_P trophiccoherenceR_P;
            sub=1;
            for nsub=idx_P_Nresp'
                gamma=hierarchicallevels_sub_NR_P(sub,nsig);
                trophiccoherenceNR_P(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_resp'
                gamma=hierarchicallevels_sub_R_P(sub,nsig);
                trophiccoherenceR_P(sub,:)=gamma;
                sub=sub+1;
            end
            a=trophiccoherenceNR_P;
            b=trophiccoherenceR_P;
            xxdata=zscore([a;b]);
        end
    elseif CASE==4
        Ceff1=CeffEB;
        Ceff2=CeffEB;
        NSUB1=length(idx_E_Nresp);
        NSUB2=length(idx_E_resp);
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

        renderingdiff(CASE,:)=mean(hierarchicallevels_sub_R_E)-mean(hierarchicallevels_sub_NR_E);

        for n=1:N
            clear trophiccoherenceNR_E trophiccoherenceR_E;
            sub=1;
            for nsub=idx_E_Nresp'
                gamma=hierarchicallevels_sub_NR_E(sub,n);
                trophiccoherenceNR_E(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_E_resp'
                gamma=hierarchicallevels_sub_R_E(sub,n);
                trophiccoherenceR_E(sub,:)=gamma;
                sub=sub+1;
            end
            a=mean((trophiccoherenceNR_E).^2,2)';
            b=mean((trophiccoherenceR_E).^2,2)';
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ranksum');
            pp(n)=min(stats.pvals);
        end
        [aa1 indx] = sort(pp,'ascend');
        nTot = 6;
        nsig=indx(1:nTot);
        
        if length(nsig)>0

            clear trophiccoherenceNR_E trophiccoherenceR_E;
            sub=1;
            for nsub=idx_E_Nresp'
                gamma=hierarchicallevels_sub_NR_E(sub,nsig);
                trophiccoherenceNR_E(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_E_resp'
                gamma=hierarchicallevels_sub_R_E(sub,nsig);
                trophiccoherenceR_E(sub,:)=gamma;
                sub=sub+1;
            end
            a=trophiccoherenceNR_E;
            b=trophiccoherenceR_E;
            xxdata=zscore([a;b]);
            %%
            baseE(idx_E_Nresp')=mean(a,2);
            baseE(idx_E_resp')=mean(b,2);
            [cc ppv]=corrcoef(baseE',BDIscore(idx_E));
            cc(1,2)
            ppv(1,2)
            %%
        end
    elseif CASE==5
        Ceff1=CeffPA;
        Ceff2=CeffPA;
        NSUB1=length(idx_P_Nresp);
        NSUB2=length(idx_P_resp);
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
            hierarchicallevels_sub_NR_PA(sub,:)=gamma';
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
            hierarchicallevels_sub_R_PA(sub,:)=gamma';
            sub=sub+1;
        end

        renderingdiff(CASE,:)=mean(hierarchicallevels_sub_R_PA)-mean(hierarchicallevels_sub_NR_PA);


        for n=1:N
            clear trophiccoherenceNR_PA trophiccoherenceR_PA;
            sub=1;
            for nsub=idx_P_Nresp'
                gamma=hierarchicallevels_sub_NR_PA(sub,n);
                trophiccoherenceNR_PA(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_resp'
                gamma=hierarchicallevels_sub_R_PA(sub,n);
                trophiccoherenceR_PA(sub,:)=gamma;
                sub=sub+1;
            end
            a=mean((trophiccoherenceNR_PA).^2,2)';
            b=mean((trophiccoherenceR_PA).^2,2)';
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ranksum');
            pp(n)=min(stats.pvals);
        end

        nsig=FDR_benjHoch(pp,0.05);
        if length(nsig)>0

            clear trophiccoherenceNR_PA trophiccoherenceR_PA;
            sub=1;
            for nsub=idx_P_Nresp'
                gamma=hierarchicallevels_sub_NR_PA(sub,nsig);
                trophiccoherenceNR_PA(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_resp'
                gamma=hierarchicallevels_sub_R_PA(sub,nsig);
                trophiccoherenceR_PA(sub,:)=gamma;
                sub=sub+1;
            end
            a=trophiccoherenceNR_PA;
            b=trophiccoherenceR_PA;
            xxdata=zscore([a;b]);
        end
    elseif CASE==6
        Ceff1=CeffEA;
        Ceff2=CeffEA;
        NSUB1=length(idx_E_Nresp);
        NSUB2=length(idx_E_resp);
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
            hierarchicallevels_sub_NR_EA(sub,:)=gamma';
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
            hierarchicallevels_sub_R_EA(sub,:)=gamma';
            sub=sub+1;
        end

        renderingdiff(CASE,:)=mean(hierarchicallevels_sub_R_EA)-mean(hierarchicallevels_sub_NR_EA);


        for n=1:N
            clear trophiccoherenceNR_EA trophiccoherenceR_EA;
            sub=1;
            for nsub=idx_E_Nresp'
                gamma=hierarchicallevels_sub_NR_EA(sub,n);
                trophiccoherenceNR_EA(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_E_resp'
                gamma=hierarchicallevels_sub_R_EA(sub,n);
                trophiccoherenceR_EA(sub,:)=gamma;
                sub=sub+1;
            end
            a=mean((trophiccoherenceNR_EA).^2,2)';
            b=mean((trophiccoherenceR_EA).^2,2)';
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ranksum');
            pp(n)=min(stats.pvals);
        end

        nsig=FDR_benjHoch(pp,0.05);
        if length(nsig)>0

            clear trophiccoherenceNR_EA trophiccoherenceR_EA;
            sub=1;
            for nsub=idx_E_Nresp'
                gamma=hierarchicallevels_sub_NR_EA(sub,nsig);
                trophiccoherenceNR_EA(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_E_resp'
                gamma=hierarchicallevels_sub_R_EA(sub,nsig);
                trophiccoherenceR_EA(sub,:)=gamma;
                sub=sub+1;
            end
            a=trophiccoherenceNR_EA;
            b=trophiccoherenceR_EA;
            xxdata=zscore([a;b]);
        end
    elseif CASE==7
        Ceff1=CeffEA;
        Ceff2=CeffPA;
        NSUB1=length(idx_E_resp);
        NSUB2=length(idx_P_resp);
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
            hierarchicallevels_sub_EA(sub,:)=gamma';
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
            hierarchicallevels_sub_PA(sub,:)=gamma';
            sub=sub+1;
        end

        renderingdiff(CASE,:)=mean(hierarchicallevels_sub_PA)-mean(hierarchicallevels_sub_EA);

        for n=1:N
            clear trophiccoherenceEA trophiccoherencePA;
            sub=1;
            for nsub=idx_E_resp'
                gamma=hierarchicallevels_sub_EA(sub,n);
                trophiccoherenceEA(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_resp'
                gamma=hierarchicallevels_sub_PA(sub,n);
                trophiccoherencePA(sub,:)=gamma;
                sub=sub+1;
            end
            a=mean((trophiccoherenceEA).^2,2)';
            b=mean((trophiccoherencePA).^2,2)';
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ranksum');
            pp(n)=min(stats.pvals);
        end

        [aa1 indx] = sort(pp,'ascend');
        nTot = 33;
        nsig=indx(1:nTot);
        
        if length(nsig)>0

            clear trophiccoherenceEA trophiccoherencePA;
            sub=1;
            for nsub=idx_E_resp'
                gamma=hierarchicallevels_sub_EA(sub,nsig);
                trophiccoherenceEA(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=idx_P_resp'
                gamma=hierarchicallevels_sub_PA(sub,nsig);
                trophiccoherencePA(sub,:)=gamma;
                sub=sub+1;
            end
            a=trophiccoherenceEA;
            b=trophiccoherencePA;
            xxdata=zscore([a;b]);
        end
    elseif CASE==8
        Ceff1=CeffPB;
        Ceff2=CeffPA;
        NSUB1=length(idx_P);
        NSUB2=length(idx_P);
        sub=1;
        for nsub=1:NSUB1
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
            hierarchicallevels_sub_PBall(sub,:)=gamma';
            sub=sub+1;
        end

        sub=1;
        for nsub=1:NSUB1
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
            hierarchicallevels_sub_PAall(sub,:)=gamma';
            sub=sub+1;
        end

        renderingdiff(CASE,:)=mean(hierarchicallevels_sub_PAall)-mean(hierarchicallevels_sub_PBall);

        for n=1:N
            clear trophiccoherencePAall trophiccoherencePBall;
            sub=1;
            for nsub=1:NSUB1
                gamma=hierarchicallevels_sub_PAall(sub,n);
                trophiccoherencePAall(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=1:NSUB1
                gamma=hierarchicallevels_sub_PBall(sub,n);
                trophiccoherencePBall(sub,:)=gamma;
                sub=sub+1;
            end
            a=mean((trophiccoherencePAall).^2,2)';
            b=mean((trophiccoherencePBall).^2,2)';
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ranksum');
            pp(n)=min(stats.pvals);
        end

        [aa1 indx] = sort(pp,'ascend');
        nTot = 69;
        nsig=indx(1:nTot);
        
        if length(nsig)>0

            clear trophiccoherencePAall trophiccoherencePBall;
            sub=1;
            for nsub=1:NSUB1
                gamma=hierarchicallevels_sub_PAall(sub,nsig);
                trophiccoherencePAall(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=1:NSUB1
                gamma=hierarchicallevels_sub_PBall(sub,nsig);
                trophiccoherencePBall(sub,:)=gamma;
                sub=sub+1;
            end
            a=trophiccoherencePAall;
            b=trophiccoherencePBall;
            xxdata=zscore([a;b]);
        end
    elseif CASE==9
        Ceff1=CeffEB;
        Ceff2=CeffEA;
        NSUB1=length(idx_E);
        NSUB2=length(idx_E);
        sub=1;
        for nsub=1:NSUB1
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
            hierarchicallevels_sub_EBall(sub,:)=gamma';
            sub=sub+1;
        end

        sub=1;
        for nsub=1:NSUB1
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
            hierarchicallevels_sub_EAall(sub,:)=gamma';
            sub=sub+1;
        end

        renderingdiff(CASE,:)=mean(hierarchicallevels_sub_EAall)-mean(hierarchicallevels_sub_EBall);

        for n=1:N
            clear trophiccoherenceEAall trophiccoherenceEBall;
            sub=1;
            for nsub=1:NSUB1
                gamma=hierarchicallevels_sub_EAall(sub,n);
                trophiccoherenceEAall(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=1:NSUB1
                gamma=hierarchicallevels_sub_EBall(sub,n);
                trophiccoherenceEBall(sub,:)=gamma;
                sub=sub+1;
            end
            a=mean((trophiccoherenceEAall).^2,2)';
            b=mean((trophiccoherenceEBall).^2,2)';
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ranksum');
            pp(n)=min(stats.pvals);
        end

        [aa1 indx] = sort(pp,'ascend');
        nTot = 6;
        nsig=indx(1:nTot); 
        
        if length(nsig)>0

            clear trophiccoherenceEAall trophiccoherenceEBall;
            sub=1;
            for nsub=1:NSUB1
                gamma=hierarchicallevels_sub_EAall(sub,nsig);
                trophiccoherenceEAall(sub,:)=gamma;
                sub=sub+1;
            end
            sub=1;
            for nsub=1:NSUB1
                gamma=hierarchicallevels_sub_EBall(sub,nsig);
                trophiccoherenceEBall(sub,:)=gamma;
                sub=sub+1;
            end
            a=trophiccoherenceEAall;
            b=trophiccoherenceEBall;
            xxdata=zscore([a;b]);
        end
    end

    if length(nsig)>0
        vecnsig{CASE}=nsig;
        %%
        NSUB=min(NSUB1,NSUB2);
        DataAll1=xxdata(1:NSUB1,:);
        DataAll2=xxdata(NSUB1+1:NSUB1+NSUB2,:);
        %%

        cl=1:2;
        pc=zeros(2,2);
        pc2=zeros(2,2);
        for nfold=1:kfold
            shuffling=randperm(NSUB1);
            Data=DataAll1;
            Data=Data(shuffling(1:NSUB),:);
            TrainData1=Data(1:NSUB-1,:);
            XValidation1=Data(NSUB,:);
            Responses1=categorical(ones(NSUB-1,1),cl);
            YValidation1=categorical(ones(1,1),cl);

            shuffling=randperm(NSUB2);
            Data=DataAll2;
            Data=Data(shuffling(1:NSUB),:);
            TrainData2=Data(1:NSUB-1,:);
            XValidation2=Data(NSUB,:);
            Responses2=categorical(2*ones(NSUB-1,1),cl);
            YValidation2=categorical(2*ones(1,1),cl);

            TrainData=vertcat(TrainData1,TrainData2);
            XValidation=vertcat(XValidation1,XValidation2);
            Responses=vertcat(Responses1,Responses2);
            YValidation=vertcat(YValidation1,YValidation2);

            %% RBF
            t = templateSVM('KernelFunction','gaussian');
            svmmodel=fitcecoc(TrainData,Responses,'Learners',t);

            %% compute

            con=zeros(2,2);
            test1=predict(svmmodel,XValidation1);
            winclass=test1;
            con(1,winclass)=con(1,winclass)+1;
            test2=predict(svmmodel,XValidation2);
            winclass=test2;
            con(2,winclass)=con(2,winclass)+1;
            accdist(CASE,nfold,:)=sum(diag(con))/2;
            pc=pc+con;

            %% Linear
            t = templateSVM('KernelFunction','linear');
            svmmodel2=fitcecoc(TrainData,Responses,'Learners',t);

            %% compute

            con=zeros(2,2);
            test1=predict(svmmodel2,XValidation1);
            winclass=test1;
            con(1,winclass)=con(1,winclass)+1;
            test2=predict(svmmodel2,XValidation2);
            winclass=test2;
            con(2,winclass)=con(2,winclass)+1;
            accdist2(CASE,nfold,:)=sum(diag(con))/2;
            pc2=pc2+con;
        end
        pc=pc/kfold;
        pcmat(CASE,:,:)=pc;
        acc(CASE)=sum(diag(pc))/2

        pc2=pc2/kfold;
        pcmat2(CASE,:,:)=pc2;
        acc2(CASE)=sum(diag(pc2))/2
    end
end
acc
acc2
%%

save results_class_psilodep2_sorted.mat vecnsig renderingdiff acc pcmat accdist acc2 pcmat2 accdist2;
