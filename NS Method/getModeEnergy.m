function [Ef,Es] =  getModeEnergy (EigVect,sizeFu,sizeFp,sizeFl,sizeSu,sizeSe, FixedFu, FixedFp, FixedFl, FixedSu, FixedSe)

[XPosS,YPosS,XPosF,YPosF] = readMesh('./meshes/xs.dat','./meshes/ys.dat','./meshes/xf.dat','./meshes/yf.dat');

solidEigtempU = EigVect(sizeFu+sizeFp+sizeFl+1:sizeFu+sizeFp+sizeFl+sizeSu);
for i = 1:length(FixedSu)
    solidEigtempU = [solidEigtempU(1:(FixedSu(i)-1)); 0; solidEigtempU((FixedSu(i)):end)];
end
solidEigU = [solidEigtempU(1:2:end-1) solidEigtempU(2:2:end)];

solidEigtempE = EigVect(sizeFu+sizeFp+sizeFl+sizeSu+1:sizeFu+sizeFp+sizeFl+sizeSu+sizeSe);
for i = 1:length(FixedSe)
    solidEigtempE = [solidEigtempE(1:(FixedSe(i)-1)); 0; solidEigtempE((FixedSe(i)):end)];
end
solidEigE = [solidEigtempE(1:2:end-1) solidEigtempE(2:2:end)];

fluidEigU = EigVect(1:sizeFu);
for i = 1:length(FixedFu)
    fluidEigU = [fluidEigU(1:(FixedFu(i)-1)); 0; fluidEigU((FixedFu(i)):end)];
end

fluidEigP = EigVect(sizeFu+1:sizeFu+sizeFp);
for i = 1:length(FixedFp)
    fluidEigP = [fluidEigP(1:(FixedFp(i)-1)); 0; fluidEigP((FixedFp(i)):end)];
end

fluidEigL = EigVect(sizeFu+sizeFp+1:sizeFu+sizeFp+sizeFl);
for i = 1:length(FixedFl)
    fluidEigL = [fluidEigL(1:(FixedFl(i)-1)); 0; fluidEigL((FixedFl(i)):end)];
end


VxEig = fluidEigU(1:2:end-1);
VyEig = fluidEigU(2:2:end);


% solidEigtemp = EigVect(1:sizeS);
% for i = 1:length(FixedS)
%     solidEigtemp = [solidEigtemp(1:(FixedS(i)-1)); 0; solidEigtemp((FixedS(i)):end)];
%     solidEig = [solidEigtemp(1:2:end-1) solidEigtemp(2:2:end)];
% end

solidEigX = solidEigU(:,1);
solidEigY = solidEigU(:,2);

% fluidEigV = EigVect(sizeS+1:sizeS+sizeFv);
% for i = 1:length(FixedFv)
%     fluidEigV = [fluidEigV(1:(FixedFv(i)-1)); 0; fluidEigV((FixedFv(i)):end)];
% end
% 
% VxEig = fluidEigV(1:2:end-1);
% VyEig = fluidEigV(2:2:end);


Vx = [real(XPosF) real(YPosF) real(VxEig)];
dlmwrite('Vx.txt',Vx," ");
Vy = [real(XPosF) real(YPosF) real(VyEig)];
dlmwrite('Vy.txt',Vy," ");
Sx = [real(XPosS) real(YPosS) real(solidEigX)];
dlmwrite('Sx.txt',Sx," ");
Sy = [real(XPosS) real(YPosS) real(solidEigY)];
dlmwrite('Sy.txt',Sy," ");

system(['FreeFem++ ./EnergyInt.edp ']);

Ef = dlmread('Ef.dat');
Es = dlmread('Es.dat');