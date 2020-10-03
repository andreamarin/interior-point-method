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
x = [1/2; -1]; %ones(n,1);
mu = ones(m,1);
y = ones(m,1);
Y = diag(y);
U = diag(mu);
e = ones(m,1);

sigma = 0.5;
eta = sigma*(y'* mu)/m;
%-----------------------------------------------
% vectores para graficacion
%cnpo=[]; comp =[];

% Condiciones necesarias de primer orden
F = vertcat(Q*x - A'*mu +c, A*x - y - b, Y*U*e);
norma = norm(F);

while(norma > tol && iter < maxiter)
  % Resolvemos sistema lineal de Newton para sigma = 0
  r_x = Q*x - A'*mu + c;
  r_y = A*x - y - b;
  r_mu = Y*U*e - eta*e;

  Y_inv = inv(Y);

  % sistema lineal (5)
  AT_Yinv = A'*Y_inv;
  K = Q + AT_Yinv*U*A;
  z = -(r_x + AT_Yinv*U*r_y + AT_Yinv*r_mu);
  
  Dx = K \ z;

  % resolvemos para dy y dmu
  Dy = A*Dx + r_y;
  Dmu = - Y_inv*( U*Dy + r_mu);
  
  %-------------------------------------------------------
  % Recorte de paso
  bt = ones(m,1); 
  gm = ones(m,1);
  for k=1:m
      if (Dmu(k) < 0)
        bt(k) =  -mu(k)/Dmu(k);
      end

      if(Dy(k) < 0)
          gm(k) = -y(k)/Dy(k);
      end
  end
  
  alfa = min([bt gm]);
  alfa =(0.995)*min([1 alfa]);
  %------------------------------------------------
  % Actualizamos puntos
  x = x + alfa*Dx;
  mu = mu + alfa*Dmu;
  y = y + alfa*Dy;
  
  Y = diag(y);
  U = diag(mu);
  %-------------------------------------------------------  
  % Nueva eta
  eta = sigma*(mu'*y)/m;
  %-------------------------------------------------------  
  % Calculamos otra vez las condiciones necesarias de primer orden
  F = vertcat(Q*x - A'*mu +c, A*x - y - b, Y*U*e);
  norma = norm(F);

  % actualizamos numero de iteraciones
  iter = iter + 1;
end

disp(fprintf('finalizó en %3.0f iteraciones',iter))
    
