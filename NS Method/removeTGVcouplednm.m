function [Ms,Ks,C1s,C2s,Kf,Mf,Df,Af,A2f,Bf,B2f,Ef,Sf,Tf,Ctf,Cf,STB,FixedS,FixedFv,FixedFp] = removeTGVcouplednm(Ms,Ks,C1s,C2s,Kf,Mf,Df,Af,A2f,Bf,B2f,Ef,Sf,Tf,Ctf,Cf,STB,threshval)

FixedS = find(diag(Ms) >= threshval);
FixedFv  = find(diag(Mf) >= threshval);
FixedFp = find(diag(STB) >= threshval);

Ms( FixedS, : ) = [];
Ms( :, FixedS ) = [];

Ks( FixedS, : ) = [];
Ks( :, FixedS ) = [];

C1s( FixedS, : ) = [];
C1s( :, FixedFp ) = [];

C2s( FixedS, : ) = [];
C2s( :, FixedFv ) = [];

Kf( FixedS, : ) = [];
Kf( :, FixedS ) = [];

Mf( FixedFv, : ) = [];
Mf( :, FixedFv ) = [];

Df( FixedFv, : ) = [];
Df( :, FixedFv ) = [];

Af( FixedFv, : ) = [];
Af( :, FixedFp ) = [];

A2f( FixedFv, : ) = [];
A2f( :, FixedFv ) = [];

Bf( FixedFv, : ) = [];
Bf( :, FixedFp ) = [];

B2f( FixedFv, : ) = [];
B2f( :, FixedFv ) = [];

Ef( FixedFv, : ) = [];
Ef( :, FixedS ) = [];

Sf( FixedFp, : ) = [];
Sf( :, FixedFv ) = [];

Tf( FixedFp, : ) = [];
Tf( :, FixedFv ) = [];

Ctf( FixedFp, : ) = [];
Ctf( :, FixedS ) = [];

Cf( FixedFp, : ) = [];
Cf( :, FixedS ) = [];

STB( FixedFp, : ) = [];
STB( :, FixedFp ) = [];

return