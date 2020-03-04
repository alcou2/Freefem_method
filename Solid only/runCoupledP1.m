function [EigVal, EigVect, sizeS, FixedS] = runCoupledP1(msS,nEV)

delete('./FFMatrices/*');

% Run FreeFem code
system(['FreeFem++ ./simple_plateP1.edp -mhS ', num2str(msS)]);

%Importing FreeFem matrices
[Ms,indMs] = matrix_file_reader_Morse('./FFMatrices/Ms.dat');
[Ks,indKs] = matrix_file_reader_Morse('./FFMatrices/Ks.dat');

[Ks,Ms,FixedS] = removeTGVcoupled(Ks,Ms,1e30);

sizeS = size(Ks,1); %numbers of dof of solid domain

%%

[EigVect,EigVal] = eigs(Ks,Ms,50,'SM');
% 
EigVal = diag(EigVal);

