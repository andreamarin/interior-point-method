function [x, mu, y] = punintpc(Q,A,c,b)
% Metodo de punto interior para el problema cuadratico
% Min   (0.5)* x' * Q * x + c * x
% s.a.   A * x >= b
%
% In
% Q.- matriz nxn simetrica y positiva definida
% A.- matriz mxn
% b.- vector columna en R^m .
% c.- vector columna en R^p .
%
%Out
% x.- vector en R^n con la aproximacion del minimo local.
% mu.- vector en R^p con la aproximacion al multiplicador
%          de Lagrange asociado a la restriccion de desigualdad.
% y.- vector en R^p con la variable de holgura en la restriccion
%     de desigualdad.
%---------------------------------------------------------------------------
% Optimizacion Numerica
% 8 de octubre de 2020
% ITAM
%--------------------------------------------------------------------------------
% Parametros iniciales
tol = 1e-08;       % Tolerancia a las condiciones necesarias de 1er orden
maxiter = 250;     % maximo numero de iteraciones permitidas
iter = 0;          % contador de las iteraciones
%-----------------------------------------------------------
n = length(c);     % dimension de la variable principal
m = length(b);     % numero de restricciones de desigualdad
%----------------------------------------------------------
% variables iniciales
x = ones(n,1);
mu = ones(m,1);
y = ones(m,1);
Y = diag(y);
U = diag(mu);
e = ones(m,1);

sigma = 0.5;
eta = sigma*(y'* mu)/m;
%-----------------------------------------------
% vectores para graficacion
cnpo=[]; comp =[];

% Condiciones necesarias de primer orden
F =[Q*x - A'*mu +c; A*x - y - b; Y*U*e];
norma = norm(F);
disp('Iter      CNPO             tau ')
disp('-----------------------------------------')

while(norma > tol & iter < maxiter)
  % Resolvemos sistema lineal de Newton para sigma = 0
  r_x = Q*x - A'*mu + c
  r_y = A*x - y - b
  r_mu = Y*U*e - eta*e

  Y_inv = inv(Y)

  % sistema lineal (5)
  K = Q + A'*Y_inv*U*A
  z = -(r_x + A'*Y_inv*U*r_y + Y_inv*r_mu)
  Dx = K \ z

  % resolvemos para dy y dmu
  Dy = A*Dx + r_y;
  Dmu = - Y_inv*( U*Dy + r_mu);
  %-------------------------------------------------------
  % Recorte de paso
  bt = []; 
  gm = [];
  for k=1:m
      if (Dmu(k) < 0)
        bt = [bt; -(mu(k)/Dmu(k))];
      end

      if(Dy(k) < 0)
          gm = [gm; -(y(k)/Dy(k))];
      end
  end
  
  alfa = min([bt ; gm]);
  alfa =(0.9995)*min([1 alfa]);  
  %-----------------------------------------------------------
  % Actualizamos puntos
  x = x + alfa*Dx;
  mu = mu + alfa*Dmu;
  y = y + alfa*Dy;
  %-------------------------------------------------------  
  % Nueva eta
  eta = sigma*(mu'*y)/m;
  %-------------------------------------------------------  
  % Calculamos otra vez las condiciones necesarias de primer orden
  F =[Q*x - A'*mu +c; A*x - y - b; Y*U*e];
  norma = norm(F);

  % actualizamos número de iteraciones y arreglos
  iter = iter + 1;
  cnpo =[cnpo norma];
  comp = [comp 2*eta];
  disp(sprintf('%3.0f  %2.8f  %2.8f',iter,norma,2*eta)
end

semilogy([1:iter],cnpo,'r',[1:iter],comp,'b')
title('Convergencia de puntos interiores')
legend('CNPO', 'Complementaridad')
    
