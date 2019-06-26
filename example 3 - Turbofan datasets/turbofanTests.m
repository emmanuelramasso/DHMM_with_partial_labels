% This file allows to reproduce results of DPHMM on turbofan engines (CMAPSS)

clear all
LADATA = input('Which dataset (1 to 4) ?: ')

if LADATA == 1
    load dataset1_hi %=> charge data_train data_train_modele_morceaux
elseif LADATA==2
    load dataset2_hi
elseif LADATA == 3
    load dataset3_hi
elseif LADATA==4
    load dataset4_hi
end
disp('OC = operting conditions, useless for dataset 1 or 3')
TAKEOC = input('Take OC in modeling (1 or 0): ');

if TAKEOC
    data_train  = data_train_modele_morceaux;
end

%%%%%%%%%%% SEGMENTATION INTO STATES %%%%%%%%%%%%%
% use the same code as in "segmentation_cmapss_classification.m" available
% on my Matlab central account but also generate some plausibilities
% See the doc segmentation_cmapss_classification.m for more details
dataHMM=cell(1,length(data_train)); 
pl=cell(1,length(data_train)); 
etats=cell(1,length(data_train));
for i=1:length(data_train)
    
    x= data_train{i};% non buite => data_train_modele_morceaux
    
    r=round(20*length(x)/100);
    deb=x(1:r);
    bdeb = robustfit((1:r)'./r,deb);
    fin=x(end-r+1:end);
    bfin = robustfit((length(x)-r+1:length(x))'./length(x),fin);
    zdeb = ((1:length(x))'./length(x)) .*repmat(bdeb(2),length(x),1) + bdeb(1);
    
    % figure, plot(x)
    % hold on, plot(1:length(x),zdeb,'k')
    
    zfin = ((1:length(x))'./length(x)) .*repmat(bfin(2),length(x),1) + bfin(1);
    
    % hold on, plot(1:length(x),zfin,'k')
    
    % residue
    r=abs(zdeb - x);
    r(1:10)=-inf;
    r = smooth(r,25);
    
    % figure,plot(r)
    
    f = find(r>=10/100);
    
    etats{i} = 2+zeros(size(data_train{i}));
    etats{i}(1:f(1))=1;
    dd = f(1);
    r=abs(zfin - x);
    r = smooth(r,25);
    r(end-9:end)=+inf;
    % figure,plot(r)
    
    f = find(r<=10/100);
    f(find(f<=dd+10/100*length(r)))=[];
    if isempty(f), f = find(r<=15/100); end
    etats{i}(f(1):end)=3;
    
    % PLAUSIBILITIES
    QQ=6;
    pl{i} = zeros(QQ, length(x));
    
    pl{i}(1,find(etats{i}==1))=1;
    pl{i}(2,find(etats{i}==2))=1;
    pl{i}(3,find(etats{i}==3))=1;
    
    d=diff(etats{i}); f=find(d);
    pl{i}(4,round(f(1)/2):f(1))=1;
    pl{i}(1,round(f(1)/2):f(1))=0;
    
    pl{i}(5,round(f(2)-(f(2)-f(1))/2):f(2))=1;
    pl{i}(2,round(f(2)-(f(2)-f(1))/2):f(2))=0;
    
    pl{i}(6,round(f(2)+(length(x)-f(2))/2:length(x)))=1;
    pl{i}(3,round(f(2)+(length(x)-f(2))/2:length(x)))=0;
    
    % Data for learning
    dataHMM{i} = data_train{i}';
    
    % figure,plot(etats{i})
    
end

% figure, hold on, for i=1:100, plot(data_train{i}), end
% figure, hold on, for i=1:100, plot(etats{i}), end
% improve_figure
% xlabel('Time unit')
% ylabel('Health indicator')
% title('Dataset #1')
%%%%%%%%%%% END OF SEGMENTATION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Here begins the learning and inference phases of DPHMM
for N=[3 5 15]

    % Observation symbols, equation 23
    cl=cell(1,length(dataHMM));
    mx=-1;
    for i=1:length(dataHMM)
        cl{i} = fix(dataHMM{i}*N);
        f = min(cl{i});
        if f<=0, cl{i}=cl{i}+abs(f)+1; end
        mx=max(mx,max(cl{i}));
    end
    
    % Perf estimation
    Q=QQ;
    S=mx;clear mx
    ARI_par_nu = []; nit_par_nu = []; ll_par_nu = [];
    tps = [];
    for j=[0:0.1:1] % rho
        
        tic;

        %%%%%%%%%%% NOISY/UNCERTAIN LABELS %%%%%%%%%%%%%%
        disp('###################')
        disp(sprintf('Param rho=%f',j))
        
        disp('Initialise model...')
        if 0 
            % NON SPECIFICITY        
            pl2=pl;
            for k=1:length(pl)
                pl2{k}(find(pl2{k}==0)) = j; % non specificity
            end
            
        else
            % UNCERTAIN AND NOISY LABELS           
            for k=1:length(pl)
                if j==0,
                    perr=zeros(size(pl{k},2),1);
                elseif j==1,
                    perr=ones(size(pl{k},2),1);
                else,
                    [a,b]=param_beta(j,(0.2).^2);
                    perr=betarnd(a,b,size(pl{k},2),1);
                end;
                %[a l] = max(pl{k});% verite
                l = etats{k}; Q=3;
                [pluncertain,labelsnoisy,plnoisy] = add_noise1(l(:),perr,Q);
                if j==1
                    pluncertain = ones(size(pluncertain));
                    plnoisy = ones(size(plnoisy));
                end
                pluncertain = pluncertain./repmat(max(pluncertain,[],2),1,Q);
                Pluncertain{k}=pluncertain'; clear pluncertain
                Plnoisy{k}=plnoisy'; clear plnoisy
            end
        end
        
        %pl2 = Pluncertain;
        pl2 = Plnoisy;% better with noisy labels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Init HMM
        disp('prior...')
        % INIT
        prior0 = normalise(rand(Q,1));
        for i=1:length(pl2), prior0 = prior0 + pl2{i}(:,1);
        end
        prior0 = normalise(prior0);
        
        disp('transmat...')
        transmat0 = mk_stochastic(rand(Q,Q));
        for i=1:length(pl2), for t=2:size(pl2{i},2)
                transmat0 = transmat0 + mk_stochastic(pl2{i}(:,t-1)*pl2{i}(:,t)');
            end, end
        transmat0 = mk_stochastic(transmat0);
        
        disp('Obsmat...')
        obsmat0 = mk_stochastic(rand(Q,S));
        for i=1:length(pl2), for t=1:size(pl2{i},2)
                c = sampleDiscrete(pl2{i}(:,t));
                obsmat0(c,cl{i}(t)) = obsmat0(c,cl{i}(t))+1;
            end, end
        obsmat0 = mk_stochastic(obsmat0);
        
        disp('Learn model...')
        [LL, prior1_2, transmat1_2, obsmat2,nit] = ...
            dhmm_em_partialLabels(cl, prior0, transmat0, obsmat0, pl2, 'max_iter', 100);
        tps = [tps ; toc];
        
        disp('Infer & eval...')
        ARI=[];
        for i=1:length(dataHMM)
            obslik1_2 = multinomial_prob(cl{i}, obsmat2);
            if 0 % use posterior for decision
                [a b c] = fwdback(prior1_2, transmat1_2, obslik1_2);
                [a p1_2] = max(c); clear a b c
            else% viterbi algo.
                p1_2 = viterbi_path(prior1_2, transmat1_2, obslik1_2);
            end
            [AR,RI,MI,HI]=valid_RandIndex(p1_2,etats{i});
            ARI = [ARI ; AR];
        end
        
        % Perf
        %   figure,plot(ARI)
        %figure,boxplot(ARI),title(sprintf('Pour \rho=%f',j'))
        ARI_par_nu = [ARI_par_nu ARI];
        nit_par_nu = [nit_par_nu ; nit];
        ll_par_nu = [ll_par_nu ; LL(end)];
    end
    
    figure,bh=boxplot(ARI_par_nu);
    for i=1:size(bh,1) % <- # graphics handles/x
        %for j=1:size(bh,2)
        set(bh(i,:),'linewidth',3);
        disp(sprintf('working on component: %3d = %s',i,get(bh(i,1),'tag')));
        %set(gca,'FontSize',28)
    end
    xlabel('\rho')
    ylabel('ARI')
    xtix = {'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};   % Your labels
    xtixloc = [1:11];      % Your label locations
    set(gca,'XTickMode','auto','XTickLabel',xtix,'XTick',xtixloc);
    improve_figure% YOU SHOULD MODIFY THE XLABEL BEFORE SAVING IF ON UBUNTU
    ylim([0 1])
    
    %figure_pdf_cropped(gcf,'ari_nu_data1')
    figure_pdf_cropped(gcf,sprintf('rev3_ari_dhmm_turbofan_dataset_%d_OC%d_rho_N%d_Q%d_PlnoisyOnly',LADATA,TAKEOC,N,Q))
    save(sprintf('rev3_dhmm_turbofan_dataset%d_OC%d_rho_N%d_Q%d_PlnoisyOnly',LADATA,TAKEOC,N,Q));
end
