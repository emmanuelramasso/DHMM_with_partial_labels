% This file allows to reproduce the results of figure 2
% It makes use of ECM (evidential Cmeans, download from T. Denoeux Homepage)
% It takes several hours to run, due mainly to ECM.m, did try to optimise it

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings of the paper: PHMM or PGMM (with continuous data)
% or DPHMM with discrete data
%
MODELETYPE = 'HMM';% Hidden Markov Models
%MODELETYPE = 'GMM';% Gaussian Mixture Models
% OBSTYPE = 'discrete'; % discrete observations
OBSTYPE = 'continuous';% continuous ones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Parameters of the model

K=3; % states %% MUST REMAIN 3...
d=2; % components in the mixture %% MUST REMAIN 2...

MU=2*[2 0;0 2;1 1];
SIG=[7 7 ;  7 7 ; 2 2];% 3rd line ~ doubt
Pi=ones(K,1)/K; % prior of HMM

% in case of HMM => transition
if strcmp(MODELETYPE,'HMM')
    A=[0.6 0.1 0.3
        0.1 0.6 0.3
        0.15 0.15 0.7];
end

V = 4; % symbols;
B = rand(K,V);
B = mk_stochastic(B);

ntests = 50;
ARIECMalone = zeros(ntests,1);% ECM

ARIBBAE2M = zeros(ntests,1);% BBA
ARIPLE2M = zeros(ntests,1);% PL
ARIBETPE2M = zeros(ntests,1);% BETP
ARIMODELEINIT = zeros(ntests,1);% unsupervised
ARIBBAE2MNOISY = zeros(ntests,1);% BBA
ARIPLE2MNOISY = zeros(ntests,1);% PL
ARIBETPE2MNOISY = zeros(ntests,1);% BETP

ARIBBAE2M_GMM = zeros(ntests,1);% BBA
ARIPLE2M_GMM = zeros(ntests,1);% PL
ARIBETPE2M_GMM = zeros(ntests,1);% BETP
ARIMODELEINIT_GMM = zeros(ntests,1);% unsupervised
ARIBBAE2MNOISY_GMM = zeros(ntests,1);% BBA
ARIPLE2MNOISY_GMM = zeros(ntests,1);% PL
ARIBETPE2MNOISY_GMM = zeros(ntests,1);% BETP

ARIBBAE2M_DHMM = zeros(ntests,1);% BBA
ARIPLE2M_DHMM = zeros(ntests,1);% PL
ARIBETPE2M_DHMM = zeros(ntests,1);% BETP
ARIMODELEINIT_DHMM = zeros(ntests,1);% unsupervised
ARIBBAE2MNOISY_DHMM = zeros(ntests,1);% BBA
ARIPLE2MNOISY_DHMM = zeros(ntests,1);% PL
ARIBETPE2MNOISY_DHMM = zeros(ntests,1);% BETP

% loop over tests
for uu=1:ntests
    
    % sample a sequence
    T=1000;
    y=zeros(T,1);
    y(1)=find(mnrnd(1,Pi));
    x=zeros(T,d);
    x(1,:)=mvnrnd(MU(y(1),:),SIG(y(1),:));
    for t=2:T;
        if strcmp(MODELETYPE,'HMM')
            y(t)= find(mnrnd(1,A(y(t-1),:)));
        else
            y(t)= find(mnrnd(1,Pi));
        end
        x(t,:)=mvnrnd(MU(y(t),:),SIG(y(t),:));
    end;
    I=eye(K);
    pl0=I(y,:);% truth
    plvide=ones(size(pl0));
    
    yy=y; save yy yy
    
    %plotmatrix_mine(x,y)
    
    % run ECM to generate the BBA
    clear m g F pl BetP histJ N LL
    disp('%%%% OPTIM ECM %%%%')
    for c=1:5
        [m{c},g{c},F{c},pl{c},BetP{c},histJ{c},N{c}] = ECM(x,K,0,1,2,100,0);
        % perf based on BBA => MIN L 
        [a b]=max(m{c}(:,[2 3 5]),[],2);
        % pl => MAX L
        %[a b]=max(pl{c},[],2); 
        LL(c,1) = 1-valid_RandIndex(b,y);
    end

    % select the best model using ground truth
    [a b]=min(LL(:,1));% min si BBA, max si PL
    m = m{b}; g = g{b}; F=F{b}; pl=pl{b}; BetP=BetP{b}; histJ=histJ{b}; N=N{b};
    disp(sprintf('Best model ECM alone => perf %f',1-a))
    ARIECMalone(uu) = 1-a;
    
    %[a b]=max(pl,[],2); valid_RandIndex(b,y), plotmatrix_mine(x,b)
    %[a b]=max(m(:,[2 3 5]),[],2); valid_RandIndex(b,y), plotmatrix_mine(x,b)   
    
    
%    IL FAUT APPLIQUER GMM ET HMM ET DHMM SUR MEME DONNEES INITIALES!!!
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% CONTINUOUS
    %%%
    % HMM

    parametersAlgorithm = setHMMDefaultParameters;
    parametersAlgorithm.hmmOrgmm = 'HMM';
       
        % INIT a HMM => unsupervised
        pr = ones(size(x,1),K); 
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIMODELEINIT(uu) = valid_RandIndex(b,y);% perf
                
        % use the same init for all following test for fair comparison
        parametersAlgorithm.phmmInit.mu = parametersHMM.muf;
        parametersAlgorithm.phmmInit.sig = parametersHMM.Sigf;
        parametersAlgorithm.phmmInit.mix = parametersHMM.mixmatf;
        parametersAlgorithm.phmmInit.Pi = parametersHMM.Pif;
        parametersAlgorithm.phmmInit.A = parametersHMM.Af;
        parametersAlgorithm.phmmInit.gamma = outputsInference.gamma;
        parametersAlgorithm.phmmInit.gamma2 = outputsInference.gamma2;

        % from here use prior
        parametersAlgorithm.init = true;

        % BBA in E2M
        pr = m(:,[2 3 5]); % m
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIBBAE2M(uu) = valid_RandIndex(b,y);
        
        % PL in E2M
        pr = pl; % pl
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIPLE2M(uu) = valid_RandIndex(b,y);
        
        % BETP in E2M
        pr = BetP; % BetP
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIBETPE2M(uu) = valid_RandIndex(b,y);
                
        % BBA in E2M NOISY
        pr = m(:,[2 3 5]); % m
        [a b] = max(pr,[],2);
        pr = zeros(T,K); for t=1:T, pr(t,b(t)) = 1; end
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIBBAE2MNOISY(uu) = valid_RandIndex(b,y);
        
        % PL in E2M NOISY
        pr = pl; % pl
        [a b] = max(pr,[],2);
        pr = zeros(T,K); for t=1:T, pr(t,b(t)) = 1; end
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIPLE2MNOISY(uu) = valid_RandIndex(b,y);
        
        % BETP in E2M NOISY
        pr = BetP; % BetP
        [a b] = max(pr,[],2);
        pr = zeros(T,K); for t=1:T, pr(t,b(t)) = 1; end
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIBETPE2MNOISY(uu) = valid_RandIndex(b,y);
        
        
        %%%
    % GMM

    parametersAlgorithm = setHMMDefaultParameters;
    parametersAlgorithm.hmmOrgmm = 'GMM';
       
        % INIT a HMM => unsupervised
        pr = ones(size(x,1),K); 
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIMODELEINIT_GMM(uu) = valid_RandIndex(b,y);% perf
                
        % use the same init for all following test for fair comparison
        parametersAlgorithm.phmmInit.mu = parametersHMM.muf;
        parametersAlgorithm.phmmInit.sig = parametersHMM.Sigf;
        parametersAlgorithm.phmmInit.mix = parametersHMM.mixmatf;
        parametersAlgorithm.phmmInit.Pi = parametersHMM.Pif;
        parametersAlgorithm.phmmInit.A = parametersHMM.Af;
        parametersAlgorithm.phmmInit.gamma = outputsInference.gamma;
        parametersAlgorithm.phmmInit.gamma2 = outputsInference.gamma2;

        % from here use prior
        parametersAlgorithm.init = true;

        % BBA in E2M
        pr = m(:,[2 3 5]); % m
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIBBAE2M_GMM(uu) = valid_RandIndex(b,y);
        
        % PL in E2M
        pr = pl; % pl
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIPLE2M_GMM(uu) = valid_RandIndex(b,y);
        
        % BETP in E2M
        pr = BetP; % BetP
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIBETPE2M_GMM(uu) = valid_RandIndex(b,y);
                
        % BBA in E2M NOISY
        pr = m(:,[2 3 5]); % m
        [a b] = max(pr,[],2);
        pr = zeros(T,K); for t=1:T, pr(t,b(t)) = 1; end
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIBBAE2MNOISY_GMM(uu) = valid_RandIndex(b,y);
        
        % PL in E2M NOISY
        pr = pl; % pl
        [a b] = max(pr,[],2);
        pr = zeros(T,K); for t=1:T, pr(t,b(t)) = 1; end
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIPLE2MNOISY_GMM(uu) = valid_RandIndex(b,y);
        
        % BETP in E2M NOISY
        pr = BetP; % BetP
        [a b] = max(pr,[],2);
        pr = zeros(T,K); for t=1:T, pr(t,b(t)) = 1; end
        [parametersHMM, outputsInference] = ...
            phmm_gauss_mix_learn(x, pr, K, 1, parametersAlgorithm);
        
        [a b]=max(outputsInference.gamma,[],2);% MAP
        ARIBETPE2MNOISY_GMM(uu) = valid_RandIndex(b,y);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% DISCRETE HMM
        
        % quantization
        disp('In kmeans...')
        [a b] = kmeans(x,V,'replicates',100);
        obsquant = a';
        
        % DHMM VACUOUS
        disp('In dhmm init first model...')
        clear LL ll prior2 transmat2 obsmat2
        for c = 1:10
	      [ll prior2{c}, transmat2{c}, obsmat2{c}, ~] = dhmm_em_partialLabels({obsquant}, normalise(rand(K,1)), ...
                mk_stochastic(rand(K,K)), B, {ones(size(obsquant,2),K)'}, 'max_iter', 50,'verbose',0);
            LL(c) = ll(end);
        end
        [a b]=max(LL);
        prior2 = prior2{b}; transmat2 = transmat2{b}; obsmat2 = obsmat2{b};
        obslik = multinomial_prob(obsquant, obsmat2);
        
        [~, ~, gamma, ~, ~] = fwdback(prior2, transmat2, obslik);
        [~,states]=max(gamma);
        ARIMODELEINIT_DHMM(uu)=valid_RandIndex(states,y);
        
        % BBA in E2M
        disp('In dhmm BBA...')
        pr = m(:,[2 3 5]); % m
        [~, prior3, transmat3, obsmat3, ~] = dhmm_em_partialLabels({obsquant}, prior2, ...
            transmat2, obsmat2, {pr'}, 'max_iter', 500,'verbose',0);
        
        obslik = multinomial_prob(obsquant, obsmat3);
        [~, ~, gamma, ~, ~] = fwdback(prior3, transmat3, obslik);
        [a b]=max(gamma);
        ARIBBAE2M_DHMM(uu) = valid_RandIndex(b,y);
        
        % PL in E2M
        disp('In dhmm PL...')
        pr = pl;
        [~, prior3, transmat3, obsmat3, ~] = dhmm_em_partialLabels({obsquant}, prior2, ...
            transmat2, obsmat2, {pr'}, 'max_iter', 500,'verbose',0);
        
        obslik = multinomial_prob(obsquant, obsmat3);
        [~, ~, gamma, ~, ~] = fwdback(prior3, transmat3, obslik);
        [a b]=max(gamma);
        ARIPLE2M_DHMM(uu) = valid_RandIndex(b,y);
        
        % BETP in E2M
        disp('In dhmm BETP...')
        pr = BetP;
        [~, prior3, transmat3, obsmat3, ~] = dhmm_em_partialLabels({obsquant}, prior2, ...
            transmat2, obsmat2, {pr'}, 'max_iter', 500,'verbose',0);
        
        obslik = multinomial_prob(obsquant, obsmat3);
        [~, ~, gamma, ~, ~] = fwdback(prior3, transmat3, obslik);
        [a b]=max(gamma);
        ARIBETPE2M_DHMM(uu) = valid_RandIndex(b,y);
        
        
        % BBA in E2M NOISY
        disp('In dhmm BBA NOISY...')
        pr = m(:,[2 3 5]); % m
        [~, prior3, transmat3, obsmat3, ~] = dhmm_em_partialLabels({obsquant}, prior2, ...
            transmat2, obsmat2, {pr'}, 'max_iter', 500,'verbose',0);
        
        obslik = multinomial_prob(obsquant, obsmat3);
        [~, ~, gamma, ~, ~] = fwdback(prior3, transmat3, obslik);
        [a b]=max(gamma);
        ARIBBAE2MNOISY_DHMM(uu) = valid_RandIndex(b,y);
        
        % PL in E2M
        disp('In dhmm PL NOISY...')
        pr = pl;
        [a b] = max(pr,[],2);
        pr = zeros(T,K); for t=1:T, pr(t,b(t)) = 1; end
        [~, prior3, transmat3, obsmat3, ~] = dhmm_em_partialLabels({obsquant}, prior2, ...
            transmat2, obsmat2, {pr'}, 'max_iter', 500,'verbose',0);
        
        obslik = multinomial_prob(obsquant, obsmat3);
        [~, ~, gamma, ~, ~] = fwdback(prior3, transmat3, obslik);
        [a b]=max(gamma);
        ARIPLE2MNOISY_DHMM(uu) = valid_RandIndex(b,y);
        
        % BETP in E2M
        disp('In dhmm BETP NOISY...')
        pr = BetP;
        [a b] = max(pr,[],2);
        pr = zeros(T,K); for t=1:T, pr(t,b(t)) = 1; end
        [~, prior3, transmat3, obsmat3, ~] = dhmm_em_partialLabels({obsquant}, prior2, ...
            transmat2, obsmat2, {pr'}, 'max_iter', 500,'verbose',0);
        
        obslik = multinomial_prob(obsquant, obsmat3);
        [~, ~, gamma, ~, ~] = fwdback(prior3, transmat3, obslik);
        [a b]=max(gamma);
        ARIBETPE2MNOISY_DHMM(uu) = valid_RandIndex(b,y);
        
            
    % PLOT
    if uu>1
        ff=figure; clf(ff); bh=boxplot([ARIECMalone(1:uu) ARIBBAE2M(1:uu)  ARIPLE2M(1:uu)  ARIBETPE2M(1:uu) ...
            ARIBBAE2MNOISY(1:uu)  ARIPLE2MNOISY(1:uu)  ARIBETPE2MNOISY(1:uu) ARIMODELEINIT(1:uu)]);
        xtix = {'ECM alone','E2M+BBA','E2M+pl','E2M+BetP','E2M+BBA noisy','E2M+pl noisy','E2M+BetP noisy','Unsupervised'};
        xtixloc = [1:length(xtix)];
        set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
        ylabel(sprintf('ARI on %d runs',uu))
        improve_figure
        pause(3), close(ff), pause(1)
        
        ff=figure; clf(ff); bh=boxplot([ARIECMalone(1:uu) ARIBBAE2M_GMM(1:uu)  ARIPLE2M_GMM(1:uu)  ARIBETPE2M_GMM(1:uu) ...
            ARIBBAE2MNOISY_GMM(1:uu)  ARIPLE2MNOISY_GMM(1:uu)  ARIBETPE2MNOISY_GMM(1:uu) ARIMODELEINIT_GMM(1:uu)]);
        xtix = {'ECM alone','E2M+BBA','E2M+pl','E2M+BetP','E2M+BBA noisy','E2M+pl noisy','E2M+BetP noisy','Unsupervised'};
        xtixloc = [1:length(xtix)];
        set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
        ylabel(sprintf('ARI on %d runs',uu))
        improve_figure
        pause(3), close(ff), pause(1)
        
        ff=figure; clf(ff); bh=boxplot([ARIECMalone(1:uu) ARIBBAE2M_DHMM(1:uu)  ARIPLE2M_DHMM(1:uu)  ARIBETPE2M_DHMM(1:uu) ...
            ARIBBAE2MNOISY_DHMM(1:uu)  ARIPLE2MNOISY_DHMM(1:uu)  ARIBETPE2MNOISY_DHMM(1:uu) ARIMODELEINIT_DHMM(1:uu)]);
        xtix = {'ECM alone','E2M+BBA','E2M+pl','E2M+BetP','E2M+BBA noisy','E2M+pl noisy','E2M+BetP noisy','Unsupervised'};
        xtixloc = [1:length(xtix)];
        set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
        ylabel(sprintf('ARI on %d runs',uu))
        improve_figure
        pause(3), close(ff), pause(1)
        
    end
end

% for saving
% r=round(rand*10000);
% save(sprintf('boxp_%s_doute_%d',MODELETYPE,r));
% final plot

data = [ARIECMalone ARIBBAE2M  ARIPLE2M  ARIBETPE2M ...
    ARIBBAE2MNOISY  ARIPLE2MNOISY  ARIBETPE2MNOISY ARIMODELEINIT];
xtix = {'ECM alone','E2M+BBA','E2M+pl','E2M+BetP','E2M+BBA noisy','E2M+pl noisy','E2M+BetP noisy','Unsupervised'};
boxplot_change_labels(data, xtix, 36)
ylabel(sprintf('ARI of PHMM on %d runs',ntests))
improve_figure
figure_pdf_cropped(gcf,sprintf('boxpfig_%s_doute','PHMM'))

data = [ARIECMalone ARIBBAE2M_GMM  ARIPLE2M_GMM  ARIBETPE2M_GMM ...
    ARIBBAE2MNOISY_GMM  ARIPLE2MNOISY_GMM  ARIBETPE2MNOISY_GMM ARIMODELEINIT_GMM];
xtix = {'ECM alone','E2M+BBA','E2M+pl','E2M+BetP','E2M+BBA noisy','E2M+pl noisy','E2M+BetP noisy','Unsupervised'};
boxplot_change_labels(data, xtix, 36)
ylabel(sprintf('ARI of PGMM on %d runs',ntests))
improve_figure
figure_pdf_cropped(gcf,sprintf('boxpfig_%s_doute','PGMM'))

data = [ARIECMalone ARIBBAE2M_DHMM  ARIPLE2M_DHMM  ARIBETPE2M_DHMM ...
    ARIBBAE2MNOISY_DHMM  ARIPLE2MNOISY_DHMM  ARIBETPE2MNOISY_DHMM ARIMODELEINIT_DHMM];
xtix = {'ECM alone','E2M+BBA','E2M+pl','E2M+BetP','E2M+BBA noisy','E2M+pl noisy','E2M+BetP noisy','Unsupervised'};
boxplot_change_labels(data, xtix, 36)
ylabel(sprintf('ARI of DPHMM on %d runs',ntests))
improve_figure
figure_pdf_cropped(gcf,sprintf('boxpfig_%s_doute','DPHMM'))

data = [ARIBBAE2M  ARIPLE2M  ARIBETPE2M ARIMODELEINIT];
xtix = {'E2M+BBA','E2M+pl','E2M+BetP','Unsupervised'};
boxplot_change_labels(data, xtix, 36)
ylabel(sprintf('ARI on %d runs',ntests))
improve_figure
figure_pdf_cropped(gcf,sprintf('boxpfig_%s_%s_doute',MODELETYPE,OBSTYPE))

%  save E2M_PHMM
%   save E2M_DHMM
%   save E2M_PGMM

% save figure...
% figure_pdf_cropped(gcf,sprintf('boxpfig_%s_doute_%d',MODELETYPE,r))
