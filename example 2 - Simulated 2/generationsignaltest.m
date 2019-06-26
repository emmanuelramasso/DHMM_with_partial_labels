function [obst,pl0t,xt,yt]=generationsignaltest(Tt,A,Pi,SIG,MU,d,K,centr)

% xt=zeros(Tt,d);
% yt=zeros(Tt,1);
% yt(1)=find(mnrnd(1,Pi));
% xt(1,:)=mvnrnd(MU(yt(1),:),SIG(yt(1),:));
% for t=2:Tt;
%    yt(t)= find(mnrnd(1,A(yt(t-1),:)));
%    xt(t,:)=mvnrnd(MU(yt(t),:),SIG(yt(t),:));
% end;
% I=eye(K);
% pl0t=I(yt,:);


xt=zeros(Tt,d);
yt=zeros(Tt,1);
yt(1)=find(mnrnd(1,Pi));
xt(1,:)=samplegeneratorsq(yt(1));
for t=2:Tt;
   yt(t)= find(mnrnd(1,A(yt(t-1),:)));
   xt(t,:)=samplegeneratorsq(yt(t));
end;
I=eye(K);
pl0t=I(yt,:);

obst=infere_fr_kmeans(xt,centr);

