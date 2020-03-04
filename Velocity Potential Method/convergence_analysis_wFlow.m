figure();

clc;
clear all;
%close all;

U_bar = 10;
% mhS = [ 1, 2, 3, 4, 5, 6, 7, 8];
% mhF = [ 6, 8,10,12,14,16,18,20];
% mhS = [3,3,3,3,3];
mhF = [1,2,4,6,8,10,12];
mhS = 3;
% mhF = 8;

l=1; % Plate length (m)
% rho_s = 10000; %Plate density (kg/m3)
% hp = 0.01; % Plate thinckness (m)

% hp = 10.^[-1.0:-0.2:-2.0];
hp = 10.^-2.0;
rho_s = 100./hp;

nu = 0; %Plate poisson ratio
h = 1; % Fluid domain's half height
E = 1.5e8; %Plate Young modulus
D = E*hp.^3./(12*(1-nu^2));
% D = 1000;
% E = (D*12*(1-nu^2))./(hp.^3);

RL1 = 5; %Length of domain in front of the plate (in plate length unit)
RL2 = 5; %Length of domain behind of the plate (in plate length unit)

U = U_bar.* ((1/l)*sqrt(D./(rho_s.*hp))); %Dimensionnal flow speed

rho_f = 100.0; %Fluid density (kg/m3)
mu = (rho_f*l)./(rho_s.*hp);
c = h/l;
r = hp./l;

%%

for i = 1:length(mhF)
    [eV(:,i), sizeS(i), sizeF(i), FixedS, FixedF] = runCoupled_CA(U,mhS,mhF(i),rho_f,hp,rho_s,RL1,RL2,50);
    EigFreq(:,i) = eV(:,i);
    
    [triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');
    [triVecF,triIDsF,dummyF,triSupF] = mesh_reader('./meshes/fluidMesh.msh');
    
%     figure();
%     hold on;
%     trimesh(triIDsS, triVecS(:,1), triVecS(:,2),'Color','b');
%     trimesh(triIDsF, triVecF(:,1), triVecF(:,2),'Color','r');
%     grid off
%     axis equal
end

%%

for a = 1:length(mhF)
    for b = 1:length(EigFreq(:,a))
        if real(EigFreq(b,a)) < -0.1
            EigFreqAbs(b,a) = abs(real(EigFreq(b,a))) + 1i*imag(EigFreq(b,a));
        elseif abs(real(EigFreq(b,a))) < 0.1
            EigFreqAbs(b,a) = 1i*imag(EigFreq(b,a));
        else
            EigFreqAbs(b,a) = EigFreq(b,a);
        end
    end
    [EigFreqOrd(:,a), EigOrder] = sort(EigFreqAbs(:,a),'ascend','ComparisonMethod','real');
%     EigFreq_bar(:,a) = EigFreqOrd(:,a)/(sqrt(D(a)/(rho_s(a)*hp(a)))/l^2);
    EigFreq_bar(:,a) = EigFreqOrd(:,a)/(sqrt(D/(rho_s*hp))/l^2);
end

%%
 figure();
subplot(1,2,1);
hold on;
xlabel('n_{hf}');
ylabel('\omega');


plot(mhF,real(EigFreq_bar(1,:)./EigFreq_bar(1,7)),'k');
% plot(mhF,real(EigFreq_bar(3,:)),'k');



%%

clc;
clear all;
%close all;

U_bar = 10;
% mhS = [ 1, 2, 3, 4, 5, 6, 7, 8];
% mhF = [ 6, 8,10,12,14,16,18,20];
% mhS = [3,3,3,3,3];
% mhF = [1,2,4,6,8];
mhS = 3;
mhF = 8;

l=1; % Plate length (m)
% rho_s = 10000; %Plate density (kg/m3)
% hp = 0.01; % Plate thinckness (m)

% hp = 10.^[-1.0:-0.2:-2.0];
hp = 10.^-2.0;
rho_s = 100./hp;

nu = 0; %Plate poisson ratio
h = 1; % Fluid domain's half height
E = 1.5e8; %Plate Young modulus
D = E*hp.^3./(12*(1-nu^2));
% D = 1000;
% E = (D*12*(1-nu^2))./(hp.^3);

RL1 = 5; %Length of domain in front of the plate (in plate length unit)
RL2 = 10.^[-1:0.2:1]; %Length of domain behind of the plate (in plate length unit)

U = U_bar.* ((1/l)*sqrt(D./(rho_s.*hp))); %Dimensionnal flow speed

rho_f = 100.0; %Fluid density (kg/m3)
mu = (rho_f*l)./(rho_s.*hp);
c = h/l;
r = hp./l;

%%

for i = 1:length(RL2)
    [eV(:,i), sizeS(i), sizeF(i), FixedS, FixedF] = runCoupled_CA(U,mhS,mhF,rho_f,hp,rho_s,RL1,RL2(i),50);
    EigFreq(:,i) = eV(:,i);
    
    [triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');
    [triVecF,triIDsF,dummyF,triSupF] = mesh_reader('./meshes/fluidMesh.msh');
    
%     figure();
%     hold on;
%     trimesh(triIDsS, triVecS(:,1), triVecS(:,2),'Color','b');
%     trimesh(triIDsF, triVecF(:,1), triVecF(:,2),'Color','r');
%     grid off
%     axis equal
end

%%

for a = 1:length(RL2)
    for b = 1:length(EigFreq(:,a))
        if real(EigFreq(b,a)) < -0.1
            EigFreqAbs(b,a) = abs(real(EigFreq(b,a))) + 1i*imag(EigFreq(b,a));
        elseif abs(real(EigFreq(b,a))) < 0.1
            EigFreqAbs(b,a) = 1i*imag(EigFreq(b,a));
        else
            EigFreqAbs(b,a) = EigFreq(b,a);
        end
    end
    [EigFreqOrd(:,a), EigOrder] = sort(EigFreqAbs(:,a),'ascend','ComparisonMethod','real');
%     EigFreq_bar(:,a) = EigFreqOrd(:,a)/(sqrt(D(a)/(rho_s(a)*hp(a)))/l^2);
    EigFreq_bar(:,a) = EigFreqOrd(:,a)/(sqrt(D/(rho_s*hp))/l^2);
end

%%

subplot(1,2,2);
hold on;
xlabel('r_{Lo}');
ylabel('\omega');


plot(RL2,real(EigFreq_bar(1,:)),'k');
% plot(RL2,real(EigFreq_bar(3,:)),'k');

