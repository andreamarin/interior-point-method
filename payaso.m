% cargar la imagen del payaso
load clown
colormap('gray')
image(X)
surf(X)
% resolvemos el problema 
k = 5;
[W,H] = descenso2pasos(X,k);
%[W,H] = nnmf(X,k);
%surf(W*H)
image(W*H)