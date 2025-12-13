% Function y=combin(n,k)
% Calculates C_k^n, ie., the number of combinations n, k to k

function y=combin(n,k)
y=gamma(n+1)./gamma(k+1)./gamma(n-k+1);