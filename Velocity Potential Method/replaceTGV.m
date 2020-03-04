function Mat = replaceTGV(MatIn, threshval, finval)

%This function removes the very big values used to
% impose the Dirichlet BC. 

pos  = find(diag(MatIn) >= threshval);
ncut = numel(pos);

Mat = MatIn;


for jj = 1:ncut

  posnow = pos(jj);

  Mat(:,posnow)      = 0;
  Mat(posnow,:)      = 0; 
  Mat(posnow,posnow) = finval; 

end


return