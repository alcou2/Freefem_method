function [EigVal, sizeS, sizeF, FixedS, FixedF] = runCoupled(U,mhS,mhF,rho_f,hp,rho_s,RL1,RL2,nEV)

% Remove old files
delete('./FFMatrices/*');

% Run FreeFem code
system(['FreeFem++ ./coupled_potentielVit5.edp -U ' num2str(U), ' -mhS ', num2str(mhS), ' -mhF ', num2str(mhF),' -rho_f ', num2str(rho_f),' -Tpl ', num2str(hp), ' -rho_s ', num2str(rho_s) ,' -RL1 ', num2str(RL1), ' -RL2 ', num2str(RL2)]);

%%

%Importing FreeFem matrices
[Ms,indMs] = matrix_file_reader_Morse('./FFMatrices/Ms.dat');
[Ks,indKs] = matrix_file_reader_Morse('./FFMatrices/Ks.dat');
[MAfs,indMAfs] = matrix_file_reader_Morse('./FFMatrices/MAfs.dat');
[Mf,indMf] = matrix_file_reader_Morse('./FFMatrices/Mf.dat');
[Cfs,indCfs] = matrix_file_reader_Morse('./FFMatrices/Cfs.dat');
[Kfs,indKfs] = matrix_file_reader_Morse('./FFMatrices/Kfs.dat');

[Ks,Ms,Mf,MAfs,Cfs,Kfs,FixedS,FixedF] = removeTGVcoupled(Ks,Ms,Mf,MAfs,Cfs,Kfs,1e20);

%%

sizeS = size(Ks,1); %numbers of dof of solid domain
sizeF = size(Mf,1); %numbers of dof of fluid domain

%Combining the solid and fluid domains matrices
MM = [[-Ms,zeros(sizeS,sizeF)];[zeros(sizeF,sizeS),zeros(sizeF,sizeF)]];
% CC = [[zeros(sizeS,sizeS),1i.*MAfs];[-1i.*MAfs',zeros(sizeF,sizeF)]];
CC = [[zeros(sizeS,sizeS),2i.*MAfs];[-1i.*MAfs',zeros(sizeF,sizeF)]];
% KK = [[Ks,Cfs];[-Kfs',Mf]];
KK = [[Ks,2.*Cfs];[-Kfs,Mf]];


nDOF = size(MM,1);

%Reduce to a first order problem
% MMM = -[[MM,zeros(nDOF,nDOF)];[zeros(nDOF,nDOF),KK]];
MMM = [[MM,zeros(nDOF,nDOF)];[zeros(nDOF,nDOF),KK]];
KKK = [[CC,KK];[-KK,zeros(nDOF,nDOF)]];

%%

% [X,e] = polyeig(KK,CC,MM);

[EigVect,EigVal] = eigs(KKK,MMM,nEV,'SM');

EigVal = diag(EigVal);

return