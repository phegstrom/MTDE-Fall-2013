
function H = polarimagesc(ang,r,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funci�n para representar en coordenadas polares los valores
% contenidos en una matriz T, en la que cada fila corresponde
% a un �ngulo de exploraci�n distinto, y cada columna a una
% distancia de exploraci�n diferente.
% 
% * ang debe ser un vector fila con las referencias angulares
%   para las filas de la matriz T
% * r debe ser un vector columna con las referencias de distancias
%   para las columnas de la matriz T
% * T es una matriz length(ang)xlength(r) conteniendo los valores
%   a representar por la funci�n, de modo tal que cada fila 
%   de T corresponde a un �ngulo diferente, y cada columna de T
%   corresponde a una distancia de exploraci�n diferente
%
% JAG - �ltima revisi�n 11 de noviembre de 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = figure
polar(0,r(end),'w')
hold on

[X,Y] = pol2cart(ang'*ones(1,length(r)),ones(length(ang),1)*r);
scatter(X(:),Y(:),15,T(:),'filled');
colormap(1-gray)
colorbar
