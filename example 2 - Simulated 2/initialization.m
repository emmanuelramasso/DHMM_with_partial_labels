function [T, CL, A, Pi, SIG, MU, d, K]=initialization

K=3;
d=3;
MU=2*[1 0 0;0 1 0; 0 0 1];
SIG=[1 1 1
     1 1 1
     1 1 1];
Pi=ones(K,1)/K;
A=[0.6 0.3 0.1
   0.1 0.6 0.3
   0.1 0.3 0.6];
CL=50;
T=100;