function [MatKs,MatMs,FixedS] = removeTGVcoupled(MatKs,MatMs,threshval)

posKs = find(diag(MatKs) >= threshval);


MatKs( posKs, : ) = [];
MatKs( :, posKs ) = [];
MatMs( posKs, : ) = [];
MatMs( :, posKs ) = []; 
FixedS = posKs;


return