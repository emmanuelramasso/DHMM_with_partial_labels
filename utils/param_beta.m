function [a,b]=param_beta(m,v)

% paramètres a et b de la loi Beta en fonction de l'espérance m et d ela
% variance v

a=m.^2*(1-m)/v -m;
b=m*(1-m)/v -1 - a;

