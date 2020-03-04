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




function [MatOut,indexMat] = matrix_file_reader_Morse(FileName)

% This function reads the FreeFem++ sparse matrix output,
% loads the element position and values and creates the 
% sparse MatOut matrix via spconvert() function.



% ---- Opening the file ----

fid = fopen(FileName);
fprintf(['Reading file: ', FileName, '...']);



% ---- Skipping 3 lines ----

tmp_line = fgetl(fid);
tmp_line = fgetl(fid);
tmp_line = fgetl(fid);

% ---- Reading Matrix properties ----
tmp_line = fgetl(fid);
num_tmp_line = str2num(tmp_line);
n_mat    = num_tmp_line(1);
m_mat    = num_tmp_line(2);
coef_num = num_tmp_line(4);

% ---- Initializing Data Vectors ----
ijAij = zeros(coef_num+1,3);

ijAij(1,:) = [n_mat m_mat 0]; % to make sure the full matrix is loaded

% ---- Now reading every entry... ----
%while ischar(tmp_line)

formatSpec = '%d %d %f';
A = fscanf(fid,formatSpec,size(ijAij'));
ijAij(2:end,:) = A';

% for coeff = 1:coef_num
% 
%    tmp_line = fgetl(fid);
%    num_tmp_line = str2num(tmp_line);
% 
%    ijAij(coeff+1,:) = num_tmp_line;   
% 
% %   pause();
% 
% end



% ---- Closing File ----

fclose(fid);
fprintf(' Done reading!\n');



% ---- Converting to Sparse Matrix ----

MatOut = spconvert(ijAij);
indexMat = ijAij;






return
