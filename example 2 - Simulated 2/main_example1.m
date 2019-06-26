function [resnu,resrho]=main_example1(CL,T)
% Reproduce results of paper
%  for nbClusters=[3 5 10 15 20 25 30]
%     [resnu,resrho]=main_example1(nbClusters,100);
%     save(sprintf('results_dphmm_%dClusters_T%d',nbClusters,T))
% end
% The same for T=300
%
% Authors: Pablo Juesas, Emmanuel Ramasso
% 2015
%

resnu=[];
resrho=[];

K=3;
d=3;
MU=1*[1 0 0;0 1 0; 0 0 1];
SIG=[1 1 1
    1 1 1
    1 1 1];
Pi=ones(K,1)/K;
A=[0.6 0.3 0.1
    0.1 0.6 0.3
    0.1 0.3 0.6];
Tt=1000;

for ind=1:1:30
    
    %Initial parameters
    %[T, CL, A, Pi, SIG, MU, d, K]= initialization;
    
    %Signal generation
    [obs,pl0,x,y,centr]= generationsignal (T, CL, A, Pi, SIG, MU, d, K);
    
    %Test
    [obst,pl0t,xt,yt]= generationsignaltest (Tt, A, Pi, SIG, MU, d, K, centr);
    
    %Initial mode generation
    [prior1, transmat1, obsmat1]= initialmodel(K,CL);
    
    %nu
    resnu=nudhmm(T,K,y,yt,x,obs,obst,prior1,transmat1,obsmat1,ind,resnu);
    
    %rho
    resrho=rhodhmm(T,K,y,yt,x,obs,obst,prior1,transmat1,obsmat1,ind,resrho);
    
    disp(ind)
end
%Ploting
%plotres(resnu,resrho)
%save(sprintf('results_dphmm_%dClusters_T%d',CL,T),'resnu','resrho')













