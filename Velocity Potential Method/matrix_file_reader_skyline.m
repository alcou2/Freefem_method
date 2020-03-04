


function [MatOut,indexMat] = matrix_file_reader_skyline(FileName)

% This function reads FreeFem++ skyline symmetric matrices
% and outputs a Matlab sparse matrix. 

% ---- Opening the file ----

fid = fopen(FileName);
fprintf(['Reading file: ', FileName, '...']);

% ---- Reading Matrix properties ----

tmp_line = fgetl(fid);
tmp_line_cut = strsplit(tmp_line);
j=1;
for i = 1:length(tmp_line_cut)
    temp = str2num(tmp_line_cut{i});
    if(~isempty(temp)) 
        num_tmp_line(j) = temp;
        j = j+1;
    end
end

n = num_tmp_line(1);
m = num_tmp_line(2);

% ---- Skipping 3 lines ----

tmp_line = fgetl(fid);
tmp_line = fgetl(fid);
tmp_line = fgetl(fid);


% ---- Verify matric is skyline symmetric ----

tmp_line = fgetl(fid);

if(string(strip(tmp_line))=="skyline symmetric")

fprintf('Matrix is skyline symmetric');

% ---- Reading the matrix ----

count = 1;
for i=1:n %Going through every row
    tmp_line = fgetl(fid);
    tmp_line_cut = strsplit(tmp_line);
    indI = str2num(tmp_line_cut{1})+1; %Row index
    numVal = str2num(tmp_line_cut{2}(2:end-1))+1; %Number of non-zero values on the row
    for j=1:numVal %Going through each non-zero entry of the row
        if(j-numVal~=0)
            indJ = str2num(tmp_line_cut{2*j+1})+1; %Column index
            value = str2num(strip(tmp_line_cut{2*j+2},';')); %Value
        else
            temp = strsplit(tmp_line_cut{2*numVal+1},':');
            indJ = str2num(temp{1})+1;
            value = str2num(temp{2});
        end
        MatIJ(count,:) = [indI indJ value];
        count = count+1;  
    end
end
  
fclose(fid);
fprintf(' Done reading!\n');

% ---- Export sparse matrix ----

MatOut = spconvert(MatIJ);
MatOut=MatOut'+triu(MatOut',1)';
indexMat = MatIJ;   

elseif(string(strip(tmp_line))=="Skyline  non symmetric")

fprintf('Matrix is skyline non symmetric');

% ---- Reading the matrix ----

count = 1;
for i=1:n %Going through every row
    tmp_line = fgetl(fid);
    tmp_line_cut = strsplit(tmp_line);
    indI = str2num(tmp_line_cut{1})+1; %Row index
    numVal = length(tmp_line_cut)-5; %Number of non-zero values on the row
    for j=1:numVal %Going through each non-zero entry of the row
        if(j-numVal~=0)
            indJ = indI-(numVal-j); %Column index
            value = str2num(tmp_line_cut{j+3}); %Value
        else
            indJ = indI;
            value = str2num(tmp_line_cut{j+4});
        end
        MatIJ(count,:) = [indI indJ value];
        count = count+1;  
    end
end
  
fclose(fid);
fprintf(' Done reading!\n');

% ---- Export sparse matrix ----

MatOut = spconvert(MatIJ);
indexMat = MatIJ;   

else
    error('Unsupport matrix type');
end

return
