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
kinit = k;
k = 0;

crit = 1;
maxiter = 5;

Wk = W0;
Hk = H0;

while crit >= tol && k < maxiter
    Hk1 = ones(kinit,p);
    A = eye(kinit);
    b = zeros(kinit,1);
    Q = Wk'*Wk; % matriz de kxk
    for j=1:p
        Xj = X(:,j); % vector de rx1
        c = -(Xj'*Wk)'; % vector de 1xk
        [Hk1(:,j), ~, ~] = punintpc(Q, A, c, b);
    end
    
    Wk1 = ones(r,kinit);
    Q = Hk1*Hk1';
    for i=1:r
        Xi = X(i,:);
        c = -(Xi*Hk1')';
        [Wk1(i,:)] = punintpc(Q, A, c, b);
    end
    
    crit = norm(Hk1 - Hk, 'fro')^2 + norm(Wk1 - Wk, 'fro')^2;
    disp('criterio')
    disp(crit)
    
    k = k+1;
    
    Wk = Wk1;
    Hk = Hk1;
end

disp(sprintf('finalizó en %3.0f iteraciones',k))

W = Wk1;
H = Hk1;

end