clear all; clc;
mhS = [1,2,3];

l=1; % Plate length (m)
rho_s = 2000; %Plate density (kg/m3)
hp = 0.05; % Plate thinckness (m)
E = 1.5e8; %Plate Young modulus
nu = 0; %Plate poisson ratio
D = E*hp^3/(12*(1-nu^2));

for i = 1:length(mhS)
    [EigFreq(:,i)] = runCoupled(mhS(i),10);
end

for a = 1:length(mhS)
%     EigFreq_bar(:,a) = EigFreq(:,a)./(sqrt(D/(rho_s*hp))/l^2);
    EigFreq_bar(:,a) = sqrt(EigFreq(:,a));
end

EigFreq_barP2 = EigFreq_bar;

for i = 1:length(mhS)
    [EigFreq(:,i)] = runCoupledP1(mhS(i),10);
end

for a = 1:length(mhS)
%     EigFreq_bar(:,a) = EigFreq(:,a)./(sqrt(D/(rho_s*hp))/l^2);
    EigFreq_bar(:,a) = sqrt(EigFreq(:,a));
end

EigFreq_barP1 = EigFreq_bar;

% EigFreq_bar = EigFreq_barP2(1)/(sqrt(D/(rho_s*hp))/l^2)
EigFreq_bar = EigFreq_barP2(1) / sqrt(E*hp^2/(12*rho_s*l^4))

alpha = [1.875,4.694,7.8];
Freq_theo = alpha.^2 * sqrt(E*hp^2/(12*rho_s*l^4));

figure(); 

subplot(1,2,1);
hold on;
xlabel('nHs');
ylabel('\omega/\omega_{th}');
plot(mhS,real(EigFreq_barP1(1,:)/Freq_theo(1)),'b');
plot(mhS,real(EigFreq_barP2(1,:)/Freq_theo(1)),'k-');
legend('P1 elements','P2 elements');
title('Mode 1');
subplot(1,2,2);
hold on;
xlabel('nHs');
ylabel('\omega/\omega_{th}');
plot(mhS,real(EigFreq_barP1(3,:)/Freq_theo(3)),'b');
plot(mhS,real(EigFreq_barP2(3,:)/Freq_theo(3)),'k-');
legend('P1 elements','P2 elements');
title('Mode 3');
%%
clear all;
mhS = 3;
EigN = 3;
[EigFreq(:), EigVect(:,:), sizeS, FixedS] = runCoupled(mhS,10);
[EigFreq(:), EigVect(:,:), sizeS, FixedS] = runCoupledP1(mhS,10);

l=1; % Plate length (m)
rho_s = 2000; %Plate density (kg/m3)
hp = 0.05; % Plate thinckness (m)
E = 1.5e8; %Plate Young modulus
nu = 0; %Plate poisson ratio
D = E*hp^3/(12*(1-nu^2));

figure();
hold on; axis equal;
xlim([-1 2])
ylim([-1 1])
title(num2str(EigFreq(EigN)))

solidEigtemp = EigVect(1:sizeS,EigN);
for i = 1:length(FixedS)
    solidEigtemp = [solidEigtemp(1:(FixedS(i)-1)); 0; solidEigtemp((FixedS(i)):end)];
    solidEig = [solidEigtemp(1:2:end-1) solidEigtemp(2:2:end)];
end

[XPosS,YPosS] = readMesh('./meshes/xs.dat','./meshes/ys.dat');

[triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');

[Cx,Cy,indC]=intersectDOF(triVecS(:,1),triVecS(:,2),XPosS,YPosS);
solidEig1x = solidEig(indC,1);
solidEig1y = solidEig(indC,2);

scaling = 0.15/max(abs(real(solidEig(:,2)))); % Deflection scaling factor
trimesh(triIDsS, triVecS(:,1)+real(solidEig1x).*scaling, triVecS(:,2)+real(solidEig1y).*scaling,'Color','r');


xlim([-3 10]);
ylim([-1 1]);
axis equal;


% [triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');
% 
% solidEigtemp = EigVect(:,EigN);
% for i = 1:length(FixedS)
%     solidEigtemp = [solidEigtemp(1:(FixedS(i)-1)); 0; solidEigtemp((FixedS(i)):end)];
%     solidEig = [solidEigtemp(1:2:end-1) solidEigtemp(2:2:end)];
% end
% 
% scaling = 1; % Deflection scaling factor
% trimesh(triIDsS, triVecS(:,1)+real(solidEig(:,1)).*scaling, triVecS(:,2)+real(solidEig(:,2)).*scaling,'Color','r');
% 
% k_n = [1.875 4.694 7.855];
% freq_analytic = k_n(EigN)^2*sqrt((D)/(rho_s*hp*l^2))
