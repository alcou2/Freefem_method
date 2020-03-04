%% potentialMethod.m
% This file is the main script for the potential method. All the different
% paramenters are defined here and it contains the code for the
% parameter sweep. 

%%

clc;
clear all;
%close all;

%% Definition of the simulation parameters
% The dimensional parameters must match those in FreeFEM code. 

% TO DO: ajuster le code FreeFEM pour que Matlab controle les paramètres, la
% géométrie et le mesh.

% Non dimensionnal flow speed (reduced speed)
% U_bar = [(0:2:6),(8:0.5:10),(12:4:20)]; % Non dimensionnal flow speed
%U_bar = (0:0.2:16);
U_bar = 0:5:15;

l=1; % Plate length (m)
rho_s = 10000; %Plate density (kg/m3)
hp = 0.01; % Plate thinckness (m)
h = 1;
E = 1.5e8; %Plate Young modulus
nu = 0; %Plate poisson ratio
D = E*hp^3/(12*(1-nu^2));

U = U_bar.* ((1/l)*sqrt(D/(rho_s*hp))); %Dimensionnal flow speed

rho_f = 100; %Fluid density (kg/m3)
mu = (rho_f*l)/(rho_s*hp);
c = h/l;

%% Extract the eigenvalues and  eigenvectors

for i = 1:length(U)
    [eVct(:,:,i), eV(:,i), sizeS, sizeF, FixedS, FixedF] = runCoupled(U(i),3,8,50);
    EigFreq(:,i) = eV(:,i);
end

% Sort the eigenvalues

% for a = 1:length(U)
%     EigFreqOrd(:,a) = sort(EigFreq(1:2:end,a),'ascend','ComparisonMethod','real');
%     count = 1;
%     for b = 1:length(EigFreqOrd(:,a))
%         if (real(EigFreqOrd(b,a))>0)
%             EigFreq_bar(count,a) = EigFreqOrd(b,a).*(sqrt(D/(rho_s*hp))/l^2);
%             count = count+1;
%         end
%     end
% end

for a = 1:length(U)
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
    EigFreq_bar(:,a) = EigFreqOrd(:,a)./(sqrt(D/(rho_s*hp))/l^2);
    EigVect(:,:,a) = eVct(:,EigOrder,a);
end

%% Display results - Eigenvalues vs reduced speed

figure(); hold on;
xlabel('Flow speed');
ylabel('Eigenfrequency');

plot(U_bar,real(EigFreq_bar(1,:)),'b','LineWidth',1);
plot(U_bar,imag(EigFreq_bar(1,:)),'g:','LineWidth',1.2);
plot(U_bar,imag(EigFreq_bar(2,:)),'g:','LineWidth',1.2);

plot(U_bar,real(EigFreq_bar(3,:)),'k','LineWidth',1);
plot(U_bar,imag(EigFreq_bar(3,:)),'r:','LineWidth',1.2);
plot(U_bar,imag(EigFreq_bar(4,:)),'r:','LineWidth',1.2);
% 
% plot(U_bar,real(EigFreq_bar(5,:)),'g','LineWidth',1);
% plot(U_bar,imag(EigFreq_bar(5,:)),'b:','LineWidth',1.2);
% plot(U_bar,imag(EigFreq_bar(6,:)),'b:','LineWidth',1.2);
% 
% plot(U_bar,real(EigFreq_bar(7,:)),'r','LineWidth',1);
% plot(U_bar,imag(EigFreq_bar(7,:)),'k:','LineWidth',1.2);
% plot(U_bar,imag(EigFreq_bar(8,:)),'k:','LineWidth',1.2);

legend('Real 1','Img 1+','Img 1-','Real 2','Img 2+','Img 2-','location','best')

% legend('Real 1','Img 1+','Img 1-','Real 2','Img 2+','Img 2-','Real 3','Img 3+','Img 3-','Real 4','Img 4+','Img 4-','location','best')

%% Import and display mesh

[triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');
[triVecF,triIDsF,dummyF,triSupF] = mesh_reader('./meshes/fluidMesh.msh');

figure();
hold on;
trimesh(triIDsS, triVecS(:,1), triVecS(:,2),'Color','b');
trimesh(triIDsF, triVecF(:,1), triVecF(:,2),'Color','r');
grid off
axis equal


%% Extracting and displaying eigenvectors
% Mesh must have been imported

EigN = 1; % Mode number to be displayed

figure();
count = 1;
for u = 1:1:length(U)

solidEigtemp = EigVect(1:sizeS,EigN*2-1,u);
for i = 1:length(FixedS)
    solidEigtemp = [solidEigtemp(1:(FixedS(i)-1)); 0; solidEigtemp((FixedS(i)):end)];
    solidEig = [solidEigtemp(1:2:end-1) solidEigtemp(2:2:end)];
end

fluidEig = EigVect(sizeS+1:sizeS+sizeF,EigN*2-1,u);
for i = 1:length(FixedF)
    fluidEig = [fluidEig(1:(FixedF(i)-1)); 0; fluidEig((FixedF(i)):end)];
end

%Display eigenvectors

subplot(length(1:1:length(U(1:end)))/3,3,count)
count = count+1;
hold on; axis equal;

%Show isopotential lines
% trisurfc( triVecF(:,1), triVecF(:,2),real(fluidEig),200);

% Show deformed solid mesh
scaling = 20; % Deflection scaling factor
trimesh(triIDsS, triVecS(:,1)+real(solidEig(:,1)).*scaling, triVecS(:,2)+real(solidEig(:,2)).*scaling,'Color','r');

% trimesh(triIDsS, triVecS(:,1)+real(solidEig(:,1)).*scaling, triVecS(:,2)+real(solidEig(:,2)).*scaling,'Color','r');
trisurf(triIDsF, triVecF(:,1), triVecF(:,2),real(fluidEig),'EdgeColor','none','FaceColor','interp');
% trimesh(triIDsS, triVecS(:,1)+real(solidEig(:,1)).*5000, triVecS(:,2)+real(solidEig(:,2)).*5000,'Color','r');
% title(strcat('Mode: ',num2str(EigN),' - U: ',num2str(U_bar(u)), ' - Freq: ', num2str(EigFreq_bar(EigN,u)) ));
xlim([-1 2])
ylim([-1 1])

end


%% Display eigenvectors normalized

[triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');
[triVecF,triIDsF,dummyF,triSupF] = mesh_reader('./meshes/fluidMesh.msh');

figure();
hold on; axis equal;
xlim([-1 2])
ylim([-1 1])

EigN = 1; % Mode number
U_num = 4; %Index of speed

Vect = EigVect(:,EigN*2-1,U_num);

Vect = sign(real(Vect)).*abs(Vect);
% Vect = Vect.*conj(Vect);
%Vect = real(Vect);

solidEigtemp = Vect(1:sizeS);
for i = 1:length(FixedS)
    solidEigtemp = [solidEigtemp(1:(FixedS(i)-1)); 0; solidEigtemp((FixedS(i)):end)];
    solidEig = [solidEigtemp(1:2:end-1) solidEigtemp(2:2:end)];
end

fluidEig = Vect(sizeS+1:sizeS+sizeF);
for i = 1:length(FixedF)
    fluidEig = [fluidEig(1:(FixedF(i)-1)); 0; fluidEig((FixedF(i)):end)];
end

% trisurf(triIDsF, triVecF(:,1), triVecF(:,2),real(fluidEig),'EdgeColor','none','FaceColor','interp');

% trisurfc( triVecF(:,1), triVecF(:,2),real(fluidEig),200);

[X Y] = meshgrid(-5:0.01:6, -1:0.01:1);
phi = griddata(triVecF(:,1), triVecF(:,2),fluidEig,X,Y);
imagesc([-5,6],[-1,1],phi);
%colormap(bone);

% Show deformed solid mesh
scaling = -20; % Deflection scaling factor
trimesh(triIDsS, triVecS(:,1)+(solidEig(:,1)).*scaling, triVecS(:,2)+(solidEig(:,2)).*scaling,'Color','r');


%% Display eigenvectors

[triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');
[triVecF,triIDsF,dummyF,triSupF] = mesh_reader('./meshes/fluidMesh.msh');

figure();
hold on; axis equal;
xlim([-1 2])
ylim([-1 1])

EigN = 1; % Mode number
U_num = 30; %Index of speed
solidEigtemp = EigVect(1:sizeS,EigN*2-1,U_num);
for i = 1:length(FixedS)
    solidEigtemp = [solidEigtemp(1:(FixedS(i)-1)); 0; solidEigtemp((FixedS(i)):end)];
    solidEig = [solidEigtemp(1:2:end-1) solidEigtemp(2:2:end)];
end

fluidEig = EigVect(sizeS+1:sizeS+sizeF,EigN*2-1,U_num);
for i = 1:length(FixedF)
    fluidEig = [fluidEig(1:(FixedF(i)-1)); 0; fluidEig((FixedF(i)):end)];
end

% trisurf(triIDsF, triVecF(:,1), triVecF(:,2),real(fluidEig),'EdgeColor','none','FaceColor','interp');

% trisurfc( triVecF(:,1), triVecF(:,2),real(fluidEig),200);

[X Y] = meshgrid(-5:0.01:6, -1:0.01:1);
phi = griddata(triVecF(:,1), triVecF(:,2),real(fluidEig),X,Y);
imagesc([-5,6],[-1,1],phi);
%colormap(bone);

% Show deformed solid mesh
scaling = 20; % Deflection scaling factor
trimesh(triIDsS, triVecS(:,1)+real(solidEig(:,1)).*scaling, triVecS(:,2)+real(solidEig(:,2)).*scaling,'Color','r');

%%

figure();
hold on; axis equal;
xlim([-1 2])
ylim([0 1])

[ux,uy] = trigradient(triVecF(:,1), triVecF(:,2),real(fluidEig));
 uxn=ux./sqrt(ux.^2+uy.^2);
 uyn=uy./sqrt(ux.^2+uy.^2);

[mx,my] = meshgrid(-1:0.01:3,-1:0.01:1);
[sx,sy] = meshgrid(-1:0.1:3,-1:0.1:1);
%streamline(stream2(triVecF(:,1), triVecF(:,2),ux,uy),sx,sy);
%streamline(triVecF(:,1), triVecF(:,2),ux,uy)
%quiver(triVecF(:,1), triVecF(:,2),uxn,uyn);


Fx = scatteredInterpolant(triVecF(:,1), triVecF(:,2),ux);
Fy = scatteredInterpolant(triVecF(:,1), triVecF(:,2),uy);

vx = Fx(mx,my);
vy = Fy(mx,my);

streamline(mx,my,vx,vy,sx,sy);

% Show deformed solid mesh
scaling = 50; % Deflection scaling factor
trimesh(triIDsS, triVecS(:,1)+real(solidEig(:,1)).*scaling, triVecS(:,2)+real(solidEig(:,2)).*scaling,'Color','r');