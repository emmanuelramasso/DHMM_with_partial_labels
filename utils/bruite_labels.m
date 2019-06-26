function [pl, y1, pl1] = bruite_labels(rho, y, K)

T = length(y);
if rho==0 % supervised
    perr=zeros(T,1);
elseif rho==1 % unsup
    perr=ones(T,1);
else
    [a,b]=param_beta(rho,(0.2).^2); % calcul des paramètres de la loi beta pour avoir une espérance rho et un écart-type 0.2 (papier de Côme)
    perr=betarnd(a,b,T,1);  % tirage de la probabilité d'erreur
end

[pl,y1,pl1]=add_noise1(y,perr,K); % tirage des labels; y= vrais labels; K = nb de classes; pl = plausibilités, y1= labels bruités, pl1= codage binaire de y1 (à utiliser pour simuler le cas supervisé);

