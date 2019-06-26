function clusters = infere_fr_kmeans(X, centres)

k=size(centres,1);
n=size(X,1);
D = zeros(k,n);
for l = 1:k
    z=(X- ones(n, 1)*centres(l, :)).^2;
    D(l,:)=(sum(z,2))';
end
[~,clusters]=min(D);
clusters=clusters';