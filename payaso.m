%---------------------------------------------------------------------------
% Proyecto 1 Optimizacion Numerica
% 8 de octubre de 2020
% ITAM
%--------------------------------------------------------------------------------
% INTEGRANTES:
% Luis Felipe Landa Lizarralde 158228
% Andrea Perez Vega 154467
% Andrea Marin Alarcon 158999

% cargar la imagen del payaso
load clown
colormap('gray')

% resolvemos el problema 
ks = [5, 20, 30, 60, 80];

% guardar los tiempos 
tiempos = [0, 0, 0, 0, 0];
tiempos_qp = [0, 0, 0, 0, 0];

for i=1:5
    k = ks(i);
    
    % correr nuestro programa con puntos interiores
    tic;
    [W,H] = descenso2pasos(X,k); 
    tEnd = toc;
    
    tiempos(i) = tEnd;
    
    % correr con quadprog
    tic;
    [Wqp,Hqp] = descenso2pasos_qp(X,k); 
    tEnd_qp = toc;
    
    tiempos_qp(i) = tEnd_qp;
    
    figure;
    colormap('gray');
    imshowpair(W*H, Wqp*Hqp, 'montage')
     
    pause(5)
end

disp('k           tiempo punint       tiempo quadprog')
for i=1:5
    disp(sprintf('%3.0f         %6.3f s         %6.3f s', k, tiempos(i), tiempos_qp(i)));
end
