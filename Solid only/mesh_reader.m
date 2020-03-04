% Copyright (C) 2015 Stefano Boccelli <dainonottambulo -at- gmail.com>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation (version 3 of the License).
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.




function [vertexMat, triangMat, boundaMat, trisup] = mesh_reader(meshname);

% This reads the mesh!!

% meshname = "mesh_sample02.msh";


% Apro il file in lettura
fid = fopen(meshname, 'r');

% Leggo la prima riga
[val] = fscanf(fid,'%d %d %d',3);

% Nvertici, Ntriangoli, Nborder 
nv = val(1);
nt = val(2);
nb = val(3);

% Inizializzo variabili mesh
vertexMat = zeros(nv,3);
triangMat = zeros(nt,3);
boundaMat = zeros(nb,2);
trivol    = zeros(nt,1);

% Ora leggo e mi salvo i vertici degli elementi

for ii = 1:nv % per ogni vertice..
  val = fscanf(fid,'%f %f %f',3);
  vertexMat(ii,1) = val(1);
  vertexMat(ii,2) = val(2);
  vertexMat(ii,3) = val(3);  
end

% Leggo ora gli elementi e salvo l'identificativo

for ii = 1:nt % per ogni elemento..
  val = fscanf(fid,'%d %d %d %d', 4);
  triangMat(ii,1) = val(1);
  triangMat(ii,2) = val(2);
  triangMat(ii,3) = val(3);
end

% Elementi di contorno 

for ii = 1:nb
  val = fscanf(fid,'%d %d %d',3);
  boundaMat(ii,1) = val(1);
  boundaMat(ii,2) = val(2);
end


% ... Done!

fclose(fid);



% ===================================================
% Mesh Processing
% ===================================================

% Calcola l'area di ogni triangolo
for idtri = 1:nt

  % ID elementi
  id1 = triangMat(idtri,1);
  id2 = triangMat(idtri,2);
  id3 = triangMat(idtri,3);
  
  % Coordinate nodi
  x1 = vertexMat(id1,1);
  y1 = vertexMat(id1,2);
  x2 = vertexMat(id2,1);
  y2 = vertexMat(id2,2);
  x3 = vertexMat(id3,1);
  y3 = vertexMat(id3,2);

  % Area
  trisup(idtri) = 0.5*abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)); 

end

return
