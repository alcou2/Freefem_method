function [A,B,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupled(U,msS,msF,nuf)

% function [EigVectOrd,EigValOrd, sizeS, sizeF, FixedS, FixedF] = runCoupled(U,msS,msF,nEV)

% Remove old files
delete('./FFMatrices/*');

% Run FreeFem code
% system(['FreeFem++ ./FF_files/NScoupled.edp -U ' num2str(U), ' -mhS ', num2str(msS), ' -mhF ', num2str(msF), ' -nuf ', num2str(nuf)]);
system(['FreeFem++ ./NScoupledP2P1.edp -mhS ', num2str(msS), ' -mhF ', num2str(msF),' -U ' num2str(U), ' -nuf ', num2str(nuf),' -rhof ' num2str(rhof),' -E ' num2str(E),' -nus ' num2str(nus)]);


%%

%Importing FreeFem matrices

%Importing FreeFem matrices

[Ms,indMs] = matrix_file_reader_Morse('./FFMatrices/Ms.dat');
[Ks,indKs] = matrix_file_reader_Morse('./FFMatrices/Ks.dat');
[Dsf,indDsf] = matrix_file_reader_Morse('./FFMatrices/Dsf.dat');
[Mf,indMf] = matrix_file_reader_Morse('./FFMatrices/Mf.dat');
[Cf,indCf] = matrix_file_reader_Morse('./FFMatrices/Cf.dat');
[Af,indAf] = matrix_file_reader_Morse('./FFMatrices/Af.dat');
[Bf,indBf] = matrix_file_reader_Morse('./FFMatrices/Bf.dat');
[Bf2,indBf2] = matrix_file_reader_Morse('./FFMatrices/Bf2.dat');
[Sf,indSf] = matrix_file_reader_Morse('./FFMatrices/Sf.dat');
[Ntfs,indNtfs] = matrix_file_reader_Morse('./FFMatrices/Ntfs.dat');
[Nfs,indNfs] = matrix_file_reader_Morse('./FFMatrices/Nfs.dat');
[STB,indSTB] = matrix_file_reader_Morse('./FFMatrices/STB.dat');

[Nsf,indNsf] = matrix_file_reader_Morse('./FFMatrices/Nsf.dat');
[Df,indDf] = matrix_file_reader_Morse('./FFMatrices/Df.dat');
[Ir,indIr] = matrix_file_reader_Morse('./FFMatrices/Ir.dat');



%% Creating matrices for eigenvalue problem

%removing TGV of Dirichlet boundary conditions
[Ks,Ms,Dsf,Mf,Cf,Af,Bf,Bf2,Sf,Ntfs,Nfs,STB,Nsf,Df,Ir,FixedS,FixedFv,FixedFp] = removeTGVcoupledP1(Ks,Ms,Dsf,Mf,Cf,Af,Bf,Bf2,Sf,Ntfs,Nfs,STB,Nsf,Df,Ir,1e20);

sizeS = size(Ks,1); 
sizeFv = size(Mf,1); 
sizeFp = size(STB,1); 


MM = [[-Ms,sparse(sizeS,sizeFv),sparse(sizeS,sizeFp)];
    [sparse(sizeFv,sizeS),sparse(sizeFv,sizeFv),sparse(sizeFv,sizeFp)];
    [sparse(sizeFp,sizeS),sparse(sizeFp,sizeFv),sparse(sizeFp,sizeFp)]];

CC = [[sparse(sizeS,sizeS),sparse(sizeS,sizeFv),sparse(sizeS,sizeFp)];
    [sparse(sizeFv,sizeS),1i.*Mf            ,sparse(sizeFv,sizeFp)];
    [-1i*Ntfs,          sparse(sizeFp,sizeFv),sparse(sizeFp,sizeFp)]];

KK = [[Ks               ,sparse(sizeS,sizeFv),Dsf];
    [sparse(sizeFv,sizeS),Cf + 1.*Df + Bf2           ,-Af+Bf];
    [-Nfs                ,Sf-1.*Nsf             ,sparse(sizeFp,sizeFp)]];

nDOF = size(MM,1);


%Reduce to a first order problem
B = -[[MM,sparse(nDOF,nDOF)];[sparse(nDOF,nDOF),KK]];
A = [[CC,KK];[-KK,sparse(nDOF,nDOF)]];

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