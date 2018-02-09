function [Y] = sorted_EVD(X)

N = size(X,1);
X = X-ones(N,1)*mean(X);
R = (X'*X)/N;
[V,D] = eig(R);
[~,s_idx] = sort(diag(D));
s_idx = flipud(s_idx);
V = V(:,s_idx);
Y = (V'*X')';