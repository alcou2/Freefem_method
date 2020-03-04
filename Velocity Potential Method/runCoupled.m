function [EigVect,EigVal, sizeS, sizeF, FixedS, FixedF] = runCoupled(U,msS,msF,nEV)

% function [EigVectOrd,EigValOrd, sizeS, sizeF, FixedS, FixedF] = runCoupled(U,msS,msF,nEV)

% Remove old files
delete('./FFMatrices/*');

% Run FreeFem code
system(['FreeFem++ ./coupled_potentiel.edp -U ' num2str(U), ' -mhS ', num2str(msS), ' -mhF ', num2str(msF)]);

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

[EigVect,EigVal] = eigs(KKK,MMM,nEV,'SM');

EigVal = diag(EigVal);

%%

% [EigVect,EigVal] = polyeig(KK,CC,MM);
% 
% count = 1;
% for i = 1:length(EigVal) 
%     if (abs(real(EigVal(i))) > 1e-3  || abs(imag(EigVal(i))) > 1e-3)  && abs(real(EigVal(i))) ~= Inf && abs(imag(EigVal(i))) ~= Inf
%         EigValClear(count) = EigVal(i);
%         EigVectClear(:,count) = EigVect(:,i);
%         count = count + 1;
%     end
% end
% EigValClear = EigValClear';
% 
% [~,idx] = sort(abs(real(EigValClear)),'ascend','ComparisonMethod','real');
% EigValOrdT = EigValClear(idx);
% EigValOrd = (EigValOrdT(1:nEV)); 
% EigVectOrdT = EigVectClear(:,idx);
% EigVectOrd = EigVectOrdT(:,1:nEV);


return