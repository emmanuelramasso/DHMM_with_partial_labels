function ind=davies_bouldin(x,obs)

    
[nr,nc]=size(x);
k=max(obs);
[st,sw,sb,S,Sinter] = valid_sumsqures(x,obs,k);

R = NaN * zeros(k,1);
dbs=zeros(1,k);
  for i = 1:k
    for j = i+1:k
      R(i,j) = (S(i) + S(j))/Sinter(i,j);
    end
    dbs(i) = max(R(i,:));
  end
db=dbs(isfinite(dbs));
ind = mean(db);