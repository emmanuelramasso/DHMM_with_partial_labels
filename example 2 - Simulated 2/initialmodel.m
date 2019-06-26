function [prior1, transmat1, obsmat1]= initialmodel(K,CL)

prior1 = normalise(rand(K,1));
transmat1 = mk_stochastic(rand(K,K));
obsmat1 = mk_stochastic(rand(K,CL));