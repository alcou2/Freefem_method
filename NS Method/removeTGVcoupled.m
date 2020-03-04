function [As1,As2,Bs1,Bs2,Asf,K,Af1,Af2,Af3,Af4,Af5,Bf,Afs,T,TestP,TestL,FixedFu,FixedFp,FixedFl,FixedSu,FixedSe] = removeTGVcoupled(As1,As2,Bs1,Bs2,Asf,K,Af1,Af2,Af3,Af4,Af5,Bf,Afs,T,TestP,TestL,threshval)

FixedFu  = find(diag(Af1) >= threshval);
FixedFp  = find(diag(TestP) >= threshval);
FixedFl  = find(diag(TestL) >= threshval);
FixedSu = find(diag(Bs1) >= threshval);
FixedSe = find(diag(Bs2) >= threshval);


As1( FixedSe, : ) = [];
As1( :, FixedSu ) = [];

As2( FixedSu, : ) = [];
As2( :, FixedSe ) = [];

Bs1( FixedSu, : ) = [];
Bs1( :, FixedSu ) = [];

Bs2( FixedSe, : ) = [];
Bs2( :, FixedSe ) = [];

Asf( FixedSu, : ) = [];
Asf( :, FixedFl ) = [];

K( FixedSu, : ) = [];
K( :, FixedSe ) = [];

Af1( FixedFu, : ) = [];
Af1( :, FixedFu ) = [];

Af2( FixedFu, : ) = [];
Af2( :, FixedFp ) = [];

Af3( FixedFp, : ) = [];
Af3( :, FixedFu ) = [];

Af4( FixedFu, : ) = [];
Af4( :, FixedFl ) = [];

Af5( FixedFl, : ) = [];
Af5( :, FixedFu ) = [];

Bf( FixedFu, : ) = [];
Bf( :, FixedFu ) = [];

Afs( FixedFl, : ) = [];
Afs( :, FixedSu ) = [];

T( FixedFl, : ) = [];
T( :, FixedSe ) = [];

TestP( FixedFp, : ) = [];
TestP( :, FixedFp ) = [];

TestL( FixedFl, : ) = [];
TestL( :, FixedFl ) = [];

return


