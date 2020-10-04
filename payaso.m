% cargar la imagen del payaso
load clown
colormap('gray')
% resolvemos el problema 
ks = [5, 20, 30, 60, 80];

resultados = [];

disp('k           tiempo punint       tiempo quadprog')
for i=1:2
    k = ks(i);
    
    % correr nuestro programa con puntos interiores
    tic;
    [W,H] = descenso2pasos(X,k); 
    tEnd = toc;
    
    resultados = [resultados; [W,H]];
    
    % correr con quadprog
    tic;
    [Wqp,Hqp] = descenso2pasos_qp(X,k); 
    tEnd_qp = toc;
    
    figure;
    colormap('gray');
    imshowpair(W*H, Wqp*Hqp, 'montage')
    
    disp(sprintf('%3.0f         %6.3f s         %6.3f s', k, tEnd, tEnd_qp))
    pause(5)
end
