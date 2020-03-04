function [A,B,As,Bs,sizeFu,sizeFp,sizeFl,sizeSu,sizeSe, FixedFu, FixedFp, FixedFl, FixedSu, FixedSe] = runCoupledP2P1(nH,Re,E,nus,rhos)

% function [EigVectOrd,EigValOrd, sizeS, sizeF, FixedS, FixedF] = runCoupled(U,msS,msF,nEV)

% Remove old files
delete('./FFMatrices/*');

% Run FreeFem code
system(['FreeFem++ ./matrix_creation.edp -nH ', num2str(nH),' -Re ' num2str(Re),' -E ', num2str(E),' -nus ', num2str(nus),' - rhos ', num2str(rhos) ]);

%%

%Importing FreeFem matrices

[As1,indAs1] = matrix_file_reader_Morse('./FFMatrices/As1.dat');
[As2,indAs2] = matrix_file_reader_Morse('./FFMatrices/As2.dat');
[Bs1,indBs1] = matrix_file_reader_Morse('./FFMatrices/Bs1.dat');
[Bs2,indBs2] = matrix_file_reader_Morse('./FFMatrices/Bs2.dat');
[Asf,indAsf] = matrix_file_reader_Morse('./FFMatrices/Asf.dat');
[K,indK] = matrix_file_reader_Morse('./FFMatrices/K.dat');
[Af1,indAf1] = matrix_file_reader_Morse('./FFMatrices/Af1.dat');
[Af2,indAf2] = matrix_file_reader_Morse('./FFMatrices/Af2.dat');
[Af3,indAf3] = matrix_file_reader_Morse('./FFMatrices/Af3.dat');
[Af4,indAf4] = matrix_file_reader_Morse('./FFMatrices/Af4.dat');
[Af5,indAf5] = matrix_file_reader_Morse('./FFMatrices/Af5.dat');
[Bf,indBf] = matrix_file_reader_Morse('./FFMatrices/Bf.dat');
[Afs,indAfs] = matrix_file_reader_Morse('./FFMatrices/Afs.dat');
[T,indT] = matrix_file_reader_Morse('./FFMatrices/T.dat');

[TestL,indTestL] = matrix_file_reader_Morse('./FFMatrices/TestL.dat');
[TestP,indTestP] = matrix_file_reader_Morse('./FFMatrices/TestP.dat');

%% Creating matrices for eigenvalue problem

% removes the zero columns and rows in 
[Asf,Af4,Af5,Bf,Afs,T,TestL,zeroL] = removeLagrangianDOF(Asf,Af4,Af5,Bf,Afs,T,TestL);

%removing TGV of Dirichlet boundary conditions
[As1,As2,Bs1,Bs2,Asf,K,Af1,Af2,Af3,Af4,Af5,Bf,Afs,T,TestP,TestL,FixedFu,FixedFp,FixedFl,FixedSu,FixedSe] = removeTGVcoupled(As1,As2,Bs1,Bs2,Asf,K,Af1,Af2,Af3,Af4,Af5,Bf,Afs,T,TestP,TestL,1e20);

sizeFu = size(Af1,1);  
sizeFp = size(TestP,1);
sizeFl = size(TestL,1);  
sizeSu = size(Bs1,1);
sizeSe = size(Bs2,1);


% A = [[Af,T+Afs];
%      [Asf,As+K]];
% B = [[Bf,                   sparse(sizeF,sizeS)];
%      [sparse(sizeS,sizeF),  Bs]];

A = [[Af1,                  -Af2,                   -Af4,                   sparse(sizeFu,sizeSu),  sparse(sizeFu,sizeSe)];
     [-Af3,                 sparse(sizeFp,sizeFp),  sparse(sizeFp,sizeFl),  sparse(sizeFp,sizeSu),  sparse(sizeFp,sizeSe)];
     [-Af5,                 sparse(sizeFl,sizeFp),  sparse(sizeFl,sizeFl),  Afs,                    -T];
     [sparse(sizeSu,sizeFu),sparse(sizeSu,sizeFp),  Asf,                    sparse(sizeSu,sizeSu),  As2-K];
     [sparse(sizeSe,sizeFu),sparse(sizeSe,sizeFp),  sparse(sizeSe,sizeFl),  -As1,                    sparse(sizeSe,sizeSe)]];
 
B = [[Bf,                   sparse(sizeFu,sizeFp),  sparse(sizeFu,sizeFl),  sparse(sizeFu,sizeSu),  sparse(sizeFu,sizeSe)];
     [sparse(sizeFp,sizeFu),sparse(sizeFp,sizeFp),  sparse(sizeFp,sizeFl),  sparse(sizeFp,sizeSu),  sparse(sizeFp,sizeSe)];
     [sparse(sizeFl,sizeFu),sparse(sizeFl,sizeFp),  sparse(sizeFl,sizeFl),  sparse(sizeFl,sizeSu),  sparse(sizeFl,sizeSe)];
     [sparse(sizeSu,sizeFu),sparse(sizeSu,sizeFp),  sparse(sizeSu,sizeFl),  Bs1,                    sparse(sizeSu,sizeSe)];
     [sparse(sizeSe,sizeFu),sparse(sizeSe,sizeFp),  sparse(sizeSe,sizeFl),  sparse(sizeSe,sizeSu),  Bs2]];

 
As = [[sparse(sizeSu,sizeSu),   As2];
     [As1,                      sparse(sizeSe,sizeSe)]];
 
Bs = [[Bs1,                     sparse(sizeSu,sizeSe)];
     [sparse(sizeSe,sizeSu),    Bs2]];
 
nDOF = size(A,1);


%%

% if (sig == 0)
%     [EigVect,EigVal] = eigs(A,B,nEV,'SM');
% end
% if (sig ~= 0)
%     [EigVect,EigVal] = eigs(A,B,nEV,sig);
% end
% 
% EigVal = diag(EigVal);
% 
% 
% [triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');
% [triVecF,triIDsF,dummyF,triSupF] = mesh_reader('./meshes/fluidMesh.msh');
% 
% if(idFSImode == 1)
%     [EigVect,EigVal] = idFSImodes(EigVect,EigVal,sizeS,sizeFv,sizeFp,triVecS,triVecF,FixedS,FixedFv);
% end

return