function [pl,w1,pl1]=add_noise1(w,perr,r)

[n,p]=size(w);
pl=zeros(n,p,max(r));
pl1=pl;

%rand('twister',5489)

w1=w;
for i=1:n,
    for j=1:p,
        x=rand;
        if x<perr(i,j),
            w1(i,j) = unidrnd(r(j));
        end;
        if isnan(w1(i,j)),
            pl(i,j,1:r(j))=ones(1,r(j));
            pl1(i,j,1:r(j))=ones(1,r(j));
        else,
%           pl(i,j,1:r(j))= perr(i,j)*ones(1,r(j));
%           pl(i,j,w1(i,j))=1;
            pl(i,j,1:r(j))= perr(i,j)*ones(1,r(j))/r(j);
           pl(i,j,w1(i,j))=pl(i,j,w1(i,j))+1-perr(i,j);
           pl1(i,j,w1(i,j))=1;
        end;
    end;
end;

pl=squeeze(pl);
pl1=squeeze(pl1);
