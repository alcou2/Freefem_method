%% idFSImodes
% This code identifies the coupled modes from all modes obtained with the
% coupled formulation. For now, it only works with P1-P1 elements.
%
% It searches for the modes whose biggest fluid velocity vectors are
% situated on the plate surface and perpendicular to it.
%
% This could probably be done more robustly using the ratio of "energy" in
% each domain.


function [EigVectFSI,EigValFSI] = idFSImodes(EigVect,EigVal,sizeS,sizeFv,sizeFp,triVecS,triVecF,FixedS,FixedFv)

if(~diff(size(EigVal)))
    EigVal = diag(EigVal);
end

if(~iscolumn(EigVal)) 
    EigVal = EigVal';
end
%%

solidEigtemp = EigVect(1:sizeS,:);
for i = 1:length(FixedS)
    solidEigtemp = [solidEigtemp(1:(FixedS(i)-1),:); zeros(1,size(solidEigtemp,2)); solidEigtemp((FixedS(i)):end,:)];
end
solidXEig = solidEigtemp(1:2:end-1,:); 
solidYEig = solidEigtemp(2:2:end,:);


fluidVEigtemp = EigVect(sizeS+1:sizeS+sizeFv,:);
for i = 1:length(FixedFv)
    fluidVEigtemp = [fluidVEigtemp(1:(FixedFv(i)-1),:); zeros(1,size(fluidVEigtemp,2)); fluidVEigtemp((FixedFv(i)):end,:)];
end
fluidXEig = fluidVEigtemp(1:2:end-1,:); 
fluidYEig = fluidVEigtemp(2:2:end,:);

fluidEigN = sqrt(fluidXEig.^2+fluidYEig.^2);
%%

posS = triVecS(:,1:2);
posF = triVecF(:,1:2);

fluidsc = fluidYEig./max(fluidXEig);
% fluidsc = fluidEigN./max(fluidEigN);
solidsc = solidYEig./max(solidYEig);
%%
k = 0.9;

% EigValFSI = EigVal(find(abs(real(max(fluidYsc(find(ismember(posF,posS,'rows')),:))))>k)');
EigValFSI = EigVal(find(abs(real(max(fluidsc(find(ismember(posF,posS,'rows')),:))))>k)');
EigVectFSI = EigVect(:,find(abs(real(max(fluidsc(find(ismember(posF,posS,'rows')),:))))>k));





