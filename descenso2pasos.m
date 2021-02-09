%---------------------------------------------------------------------------
% Proyecto 1 Optimizacion Numerica
% 8 de octubre de 2020
% ITAM
%--------------------------------------------------------------------------------
% INTEGRANTES:
% Luis Felipe Landa Lizarralde 158228
% Andrea Perez Vega 154467
% Andrea Marin Alarcon 158999
function[W,H] = descenso2pasos(X,k)
% Método de descenso en 2 pasos (resuleve el problema cuadratico con la
% funcion punintpc)
% Min ||X - WH||F^2 s.a W >=0, H>=0
% In:
% X: matriz rxp, Xij>=0
% k<<r, k<<p
% Out:
% W: matriz 
% H: matriz
% -----
%Parámetros iniciales
tol = 1e-04;
[r, p] = size(X);
W0 = ones(r,k);
H0  = ones(k,p);
kinit = k;
k = 0;

crit = 1;
maxiter = 10;

% inicializa Wk y Hk con W0 y H0
Wk = W0;
Hk = H0;
while crit >= tol && k < maxiter
    % variables necesarias para los problemas cuadraticos
    A = eye(kinit);
    b = zeros(kinit,1);
    
    %--------------------------------------------------------------------
    % Paso 1. Resolvemos para H_k+1
    Hk1 = ones(kinit,p);
    Q = Wk'*Wk; 
    for j=1:p
        % iteramos sobre las columnas y resolvemos el problema cuadrático
        Xj = X(:,j); 
        c = -(Xj'*Wk)'; 
        [Hk1(:,j), ~, ~] = punintpc(Q, A, c, b);
    end
    
    %--------------------------------------------------------------------
    % Paso 2. Resolvemos para W_k+1
    Wk1 = ones(r,kinit);
    Q = Hk1*Hk1';
    for i=1:r
        % iteramos sobre las filas y resolvemos el problema cuadrático
        Xi = X(i,:);
        c = -(Xi*Hk1')';
        [Wk1(i,:)] = punintpc(Q, A, c, b);
    end
    
    %--------------------------------------------------------------------
    % Calculamos el criterio de paro
    crit = norm(Hk1 - Hk, 'fro')^2 + norm(Wk1 - Wk, 'fro')^2;
    disp('criterio')
    disp(crit)
    
    % sumamos una iteracion
    k = k+1;
    
    % actualizamos las matrices
    Wk = Wk1;
    Hk = Hk1;
end

% Guardamos el valor final
W = Wk1;
H = Hk1;

end
