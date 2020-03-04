function [MatKs,MatMs,MatMf,MatMAfs,MatCfs,MatKfs,FixedS,FixedF] = removeTGVcoupled(MatKs,MatMs,MatMf,MatMAfs,MatCfs,MatKfs,threshval)

posKs = find(diag(MatKs) >= threshval);
posMf = find(diag(MatMf) >= threshval);

MatKs( posKs, : ) = [];
MatKs( :, posKs ) = [];
MatMs( posKs, : ) = [];
MatMs( :, posKs ) = []; 
MatMf( posMf, : ) = [];
MatMf( :, posMf ) = [];
MatMAfs( posKs, : ) = [];
MatMAfs( :, posMf ) = [];
MatCfs( posKs, : ) = [];
MatCfs( :, posMf ) = [];
MatKfs( posMf, : ) = [];
MatKfs( :, posKs ) = [];
FixedS = posKs;
FixedF = posMf;

return