function res=nudhmm(T,K,y,yt,x,obs,obst,prior1,transmat1,obsmat1,ind,res)
row=1;

y = num2cell(y', 2);
plnu=cell(length(y),1);
for nu=0:0.1:1
    
    for ex=1:length(y)
        plnuvect=zeros(T,K)+nu;
        st=y{ex};
        for i=1:T
            plnuvect(i,st(i))=1;
        end
        plnu{ex}=plnuvect';
    end
    
    tic
    
    [~, prior2, transmat2, obsmat2, ~] = dhmm_em_partialLabels(obs', prior1, transmat1, obsmat1, plnu', 'max_iter', 50,'verbose',0);
    % obsmat2 = mk_stochastic(obsmat2+0.01*rand(size(obsmat2)));
    % transmat2 = mk_stochastic(transmat2+0.01*rand(size(transmat2)));
    obslik = multinomial_prob(obst, obsmat2);
    %states = viterbi_path(prior2, transmat2, obslik)';
    [~, ~, gamma, ~, ~] = fwdback(prior1, transmat1, obslik);
    [~,states]=max(gamma);
    loglik = dhmm_logprob(obs', prior2, transmat2, obsmat2);
    ARI=valid_RandIndex(states,yt);
    time=toc;
    
    res(row,ind).nu=nu;
    res(row,ind).data=obs';
    res(row,ind).signal=x';
    res(row,ind).ARI=ARI;
    res(row,ind).timelapse=time;
    res(row,ind).calcstates=states;
    res(row,ind).realstates=yt;
    res(row,ind).prior1=prior1;
    res(row,ind).transmat1=transmat1;
    res(row,ind).obsmat1=obsmat1;
    row=row+1;
    
end
