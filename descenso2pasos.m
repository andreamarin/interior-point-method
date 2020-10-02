function[W,H] = descenso2pasos(X,k)
% Método de descenso en 2 pasos
% Min ||X - WH||F^2 s.a W >=0, H>=0
% In:
% X: matriz rxp, Xij>=0
% k<<r, k<<p
% Out:
% W: matriz 
% H: matriz

%------------------------------------
% Optimización Numérica
% 8 de octubre de 2020
% ITAM
%------------------------------------
%Parámetros iniciales

tol = 1e-04;
[r, p] = size(X);
W0 = ones(r,k);
H0  = ones(k,p);
k = 0;

crit = 0;

Wk = W0;
Hk = H0;

while crit >= tol
    Q = Wk'*Wk;
    A = Hk;
    c = X'*Wk;
    b = zeros(k);
    [Hk1] = punintpc(Q, A, c, b);
    %[Hk1] = punintpc(Wk'*Wk, Hk, X'*Wk, zeros(k));
    
    Q = 
    A = Wk;    
    [Wk1] = punintpc(Q, A, c, b);
    crit = norm(Hk1 - Hk, 'fro') + norm(Wk1 - Wk, 'fro');
    
    k = k+1;
    
    Wk = Wk1;
    Hk = Hk1;
end

W = Wk1;
H = Hk1;

end