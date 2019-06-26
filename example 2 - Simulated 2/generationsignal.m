function [obs,pl0,x,y,centr]=generationsignal(T,CL,A,Pi,SIG,MU,d,K)

% x=zeros(T,d);
% y=zeros(T,1);
% y(1)=find(mnrnd(1,Pi));
% x(1,:)=mvnrnd(MU(y(1),:),SIG(y(1),:));
% for t=2:T;
%    y(t)= find(mnrnd(1,A(y(t-1),:)));
%    x(t,:)=mvnrnd(MU(y(t),:),SIG(y(t),:));
% end;


x=zeros(T,d);
y=zeros(T,1);
y(1)=find(mnrnd(1,Pi));
x(1,:)=samplegeneratorsq(y(1));
for t=2:T;
   y(t)= find(mnrnd(1,A(y(t-1),:)));
   x(t,:)=samplegeneratorsq(y(t));
end;

%figure,scatter3(x(:,1),x(:,2),x(:,3))
I=eye(K);
pl0=I(y,:);
%plvide=ones(size(pl0));
fact=10;
ind=zeros(fact,1);
lescentres=cell(fact,1);
obst=cell(fact,1);
for u=1:fact
    [obsk, centres]=kmeans(x,CL,'MaxIter',200,'Replicates',5,'Distance','sqEuclidean', 'EmptyAction','singleton');

    ind(u)=davies_bouldinPDHMM(x,obsk);
 
    lescentres{u} = centres;
    obst{u}=obsk;
end

[a,b]=min(ind);
obs=obst{b};
centr = lescentres{b};