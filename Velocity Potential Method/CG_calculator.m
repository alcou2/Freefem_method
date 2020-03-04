% Copyright (C) 2015
% 
% Stefano Boccelli <dainonottambulo -at- gmail.com>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation (version 3 of the License).
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.



function triCG = CG_calculator(triVec,triIDs)

% Calcola i centri di massa date le coordinate

triCG = zeros(size(triIDs,1),2);
Ntri  = size(triIDs,1);


for idtri = 1:Ntri

  % ID nodi
  nodo1 = triIDs(idtri,1);
  nodo2 = triIDs(idtri,2);
  nodo3 = triIDs(idtri,3);

  % coordinate  
  x1 = triVec(nodo1,1);
  y1 = triVec(nodo1,2);
  
  x2 = triVec(nodo2,1);
  y2 = triVec(nodo2,2);
  
  x3 = triVec(nodo3,1);
  y3 = triVec(nodo3,2);

  % CG
  triCG(idtri,1) = (x1+x2+x3)/3;
  triCG(idtri,2) = (y1+y2+y3)/3;

end


return
