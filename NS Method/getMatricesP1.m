%% getMatrices
% This function returns a structure containing all the component matrices
% that are not dependant on the flow speed. These remain the same for all
% the steps in the stability analysis. 
%
% Note: The constant matrices may change if Re or Ur remains constant.

function [matrices,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = getMatricesP1(msS,msF,U,nuf,rhof,E,nus)

% function [EigVectOrd,EigValOrd, sizeS, sizeF, FixedS, FixedF] = runCoupled(U,msS,msF,nEV)

% Remove old files
delete('./FFMatrices/*');

% Run FreeFem code
system(['FreeFem++ ./NScoupled.edp -mhS ', num2str(msS), ' -mhF ', num2str(msF),' -U ' num2str(U), ' -nuf ', num2str(nuf),' -rhof ' num2str(rhof),' -E ' num2str(E),' -nus ' num2str(nus)]);

%%

%Importing FreeFem matrices

%Importing FreeFem matrices

% Ms = matrices(1);
% Ks = matrices(2);
% Dsf = matrices(3);
% Mf = matrices(4);
% [Cf,indCf] = matrix_file_reader_Morse('./FFMatrices/Cf.dat');
% Af = matrices(5);
% Bf = matrices(6);
% Sf = matrices(7);
% Ntfs = matrices(8);
% [Nfs,indNfs] = matrix_file_reader_Morse('./FFMatrices/Nfs.dat');
% STB = matrices(9);
% Nsf = matrices(10);
% Df = matrices(11);

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


[Ks,Ms,Dsf,Mf,Cf,Af,Bf,Bf2,Sf,Ntfs,Nfs,STB,Nsf,Df,Ir,FixedS,FixedFv,FixedFp] = removeTGVcouplednm(Ks,Ms,Dsf,Mf,Cf,Af,Bf,Bf2,Sf,Ntfs,Nfs,STB,Nsf,Df,Ir,1e20);

sizeS = size(Ks,1); 
sizeFv = size(Mf,1); 
sizeFp = size(STB,1); 

matrices{1} = Ms;
matrices{2} = Ks;
matrices{3} = Dsf;
matrices{4} = Mf;
matrices{5} = Cf;
matrices{6} = Af;
matrices{7} = Bf;
matrices{8} = Bf2;
matrices{9} = Sf;
matrices{10} = Ntfs;
matrices{11} = Nfs;
matrices{12} = STB;
matrices{13} = Nsf;
matrices{14} = Df;

return