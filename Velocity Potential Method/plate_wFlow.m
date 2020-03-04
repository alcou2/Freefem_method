close all;
clear all;
clc;

%%

% Remove old files
delete('./FFMatrices/*');

% Run FreeFem code
system('FreeFem++ ./simple_beam_coupled.edp -U 10');

%%

%Importing FreeFem matrices
[Ms,indMs] = matrix_file_reader_Morse('./FFMatrices/Ms.dat');
[Ks,indKs] = matrix_file_reader_Morse('./FFMatrices/Ks.dat');
[MAfs,indMAfs] = matrix_file_reader_Morse('./FFMatrices/MAfs.dat');
[Mf,indMf] = matrix_file_reader_Morse('./FFMatrices/Mf.dat');
[Cfs,indCfs] = matrix_file_reader_Morse('./FFMatrices/Cfs.dat');

[Ks,Ms,Mf,MAfs,Cfs] = removeTGVcoupled(Ks,Ms,Mf,MAfs,Cfs,1e20);

%%

sizeS = size(Ks,1); %numbers of dof of solid domain
sizeF = size(Mf,1); %numbers of dof of fluid domain

%Combining the solid and fluid domains matrices
MM = [[Ms,MAfs];[zeros(sizeF,sizeS),zeros(sizeF,sizeF)]];
CC = [[zeros(sizeS,sizeS),Cfs];[zeros(sizeF,sizeS),zeros(sizeF,sizeF)]];
KK = [[Ks,zeros(sizeS,sizeF)];[-MAfs',Mf]];

nDOF = size(MM,1);

%Reduce to a first order problem
MMM = [[MM,zeros(nDOF,nDOF)];[zeros(nDOF,nDOF),KK]];
KKK = [[CC,KK];[KK,zeros(nDOF,nDOF)]];

%%

EigVal = eigs(KKK,MMM,6,'SM')

EigFreq = sqrt(EigVal)./(2*pi)



%%

% K2 = full(K);
% M2 = full(M);
% EigVal = eig(K2,M2)
% EigFreq = sqrt(EigVal)
% EigVal = eigs(K2,M2,5,'sm')
% EigVal = diag(D)
% EigFreq = sqrt(EigVal)

