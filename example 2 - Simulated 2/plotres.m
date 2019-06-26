function plotres(resnu,resrho)

vect=.0:0.1:1;

% ARInu=reshape([resnu(:).ARI],size(resnu));
% figure('name','label imprecision','position',[400,100,800,500]),boxplot(ARInu',vect')
% xlabel('nonspecificity');
% ylabel('adjusted Rand Index');
% set(gca,'YGrid','on');
% 
% ARIrho=reshape([resrho(:).ARI],size(resrho)); % for uncertain
% figure('name','uncertain labels','position',[300,100,800,500]),boxplot(ARIrho',vect')
% xlabel('\rho');
% ylabel('adjusted Rand Index');
% set(gca,'YGrid','on');

ARIrho1=reshape([resrho(:).ARI1],size(resrho)); % noisy case
figure('name','noisy labels','position',[300,100,800,500]),boxplot(ARIrho1',vect')
xlabel('\rho');
ylabel('adjusted Rand Index');
set(gca,'YGrid','on');

%  %TIME
% timenu=reshape([resnu(:).timelapse],size(resnu));
% timerho=reshape([resrho(:).timelapse],size(resrho));
% figure('name','tnu'),boxplot(timenu',vect')
% figure('name','trho'),boxplot(timerho',vect')
%figure('name','time comparison','position',[400,100,800,500]),plot(mean(timeT100CL5rho,2),vect')