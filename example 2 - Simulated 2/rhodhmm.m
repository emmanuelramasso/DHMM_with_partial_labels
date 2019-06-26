function res=rhodhmm(T,K,y,yt,x,obs,obst,prior1,transmat1,obsmat1,ind,res)
row=1;
y = num2cell(y', 2);
pl0={ones(K,1000)};
pl=cell(length(y),1);
pl1=cell(length(y),1);
for rho=0:0.1:1
    
    for ex=1:length(y)
        
        st=y{ex};
        [plex, y1, pl1ex] = bruite_labels(rho, st', K) ;
        pl{ex}=plex';
        pl1{ex}=pl1ex';
    end
    
    tic
    % uncertain
    [LL, prior2, transmat2, obsmat2,~] = dhmm_em_partialLabels(obs', prior1, transmat1, obsmat1, pl, 'max_iter', 5,'verbose',0);
    
    % noisy
    [LL1, prior21, transmat21, obsmat21,~] = dhmm_em_partialLabels(obs', prior1, transmat1, obsmat1, pl1, 'max_iter', 5,'verbose',0);
    
    % unsupervised
    [LL0, prior0, transmat0, obsmat0,~] = dhmm_em_partialLabels(obst', prior1, transmat1, obsmat1, pl0, 'max_iter', 5,'verbose',0);
    
    time=toc;
    
    levN = 0.01;
    obsmat20 = mk_stochastic(obsmat0+levN*rand(size(obsmat0)));
    obslik0 = multinomial_prob(obst, obsmat20);
    
    obsmat2 = mk_stochastic(obsmat2+levN*rand(size(obsmat2)));
    obslik = multinomial_prob(obst, obsmat2);
    
    obsmat21 = mk_stochastic(obsmat21+levN*rand(size(obsmat21)));
    obslik1 = multinomial_prob(obst, obsmat21);
    
    [~, ~, gamma, ~, ~] = fwdback(prior0, transmat0, obslik0);
    [~,states0]=max(gamma); % unsup
    
    [~, ~, gamma, ~, ~] = fwdback(prior2, transmat2, obslik);
    [~,states]=max(gamma);% uncertain
    
    [~, ~, gamma, ~, ~] = fwdback(prior21, transmat21, obslik1);
    [~,states1]=max(gamma);% noisy
    
    % not stable for this task
    % states0 = viterbi_path(prior0, transmat0, obslik0)';
    % states = viterbi_path(prior2, transmat2, obslik)';
    % states1 = viterbi_path(prior21, transmat21, obslik1)';
    
    %hold on
    %plot(y,'r');
    
    ARI0=valid_RandIndex(states0,yt);
    ARI=valid_RandIndex(states,yt);
    ARI1=valid_RandIndex(states1,yt);
    
    
    res(row,ind).rho=rho;
    res(row,ind).data=obs';
    res(row,ind).signal=x';
    res(row,ind).ARI0=ARI0;
    res(row,ind).ARI=ARI;
    res(row,ind).ARI1=ARI1;
    res(row,ind).timelapse=time;
    res(row,ind).calcstates=states;
    res(row,ind).realstates=yt';
    res(row,ind).prior1=prior1;
    res(row,ind).transmat1=transmat1;
    res(row,ind).obsmat1=obsmat1;
    row=row+1;
    
end
