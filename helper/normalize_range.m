function [nX] = normalize_range(X)

nX = (X-min(min(X)))./(max(max(X))-min(min(X)));