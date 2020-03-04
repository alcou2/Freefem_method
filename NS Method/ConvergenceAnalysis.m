clear all; clc;

l=2; % Plate length (m)
rhos = 25; %Plate density (kg/m3)
hp = 0.02; % Plate thinckness (m)
nus = 0.4; %Plate poisson ratio

h = 0.5; % Fluid domain height
rhof = 1; % Fluid density
U = 1; % Fluid speed
%Re = [10,50,100,500,1000]; % Reynolds number for which to evaluate
Re = 100;
Re_available = [10,20,30,40,50,60,70,80,90,100,125,175,200,225,250,300,350,400,450,500,600,700,800,900,1000];
nuf = (U*rhof*2*h)./Re; % Fluid dynamic viscosity  
% nuf = 1; % Fluid dynamic viscosity  
% U = Re*nuf/(rho_f*l); % Fluid speed

% U_barMin = 1; %Min reduced speed
% U_barMax = 3; % Max reduced speed
% steps = 8; % Number of steps in the speed continuation method
% U_bar = U_barMin:(U_barMax-U_barMin)/steps:U_barMax;
% U_bar = 10;
U_bar = [4,6,8,9,10,10.4,10.8,11,11.2,11.4,11.6,11.7,11.8,11.9,12,12.1,12.2,12.5,13,14,16,18];

E = (12*(1-nus^2)*rhos*l^2*U^2)./(hp^2.*(U_bar.^2)); %Plate Young's modulus
D = E.*(hp^3/(12*(1-nus^2)));

K_b = D/(rhof*U^2*l^2);

mu = rhof*l/(rhos*hp);
c = 2*h/(l);

% nH = [1,2,3,4,5];
nH = 2;

MacCrit = 0.9;

%% Recreate pfister

l=2; % Plate length (m)
rhos = 25; %Plate density (kg/m3)
hp = 0.02; % Plate thinckness (m)
nus = 0.4; %Plate poisson ratio

h = 0.5; % Fluid domain height
rhof = 1; % Fluid density
U = 1; % Fluid speed
%Re = [10,50,100,500,1000]; % Reynolds number for which to evaluate
Re = 100;
Re_available = [10,20,30,40,50,60,70,80,90,100,125,175,200,225,250,300,350,400,450,500,600,700,800,900,1000];
nuf = (U*rhof*2*h)./Re; % Fluid dynamic viscosity 


U_bar = 13.53;
E = (12*(1-nus^2)*rhos*l^2*U^2)./(hp^2.*(U_bar.^2)); %Plate Young's modulus
D = E.*(hp^3/(12*(1-nus^2)));

K_b = D/(rhof*U^2*l^2);

mu = rhof*l/(rhos*hp);
c = 2*h/(l);

% nH = [1,2,3,4,5];
nH = 3;



%%

%system(['FreeFem++ ./NS_base.edp -nH ', num2str(nH(i)) ]);

for i = 1:length(U_bar)

    
    [A,B,As,Bs,sizeFu,sizeFp,sizeFl,sizeSu,sizeSe, FixedFu, FixedFp, FixedFl, FixedSu, FixedSe] = runCoupledP2P1(nH,Re,E(i),nus,rhos); 
    [EigVectTemp,EigValTemp] = eigs(A,-1.0i.*B,50,0);
    EigVal(:,i) = diag(EigValTemp);
    

end



for i = 1:length(U_bar)

    [temp,index] = sort(1j.*EigVal(:,i),'descend','ComparisonMethod','real');
    EigValSort(:,i) = EigVal(index,i);
    

end

%%

figure();
subplot(1,2,1);
plot(U_bar,-imag(EigValSort(1,:)));
subplot(1,2,2);
plot(U_bar,abs(real(EigValSort(1,:))));
%%

temp1 = [EigValSort(1,1:9),EigValSort(3,10:11),EigValSort(5,12:20),EigValSort(7,21),EigValSort(9,22)];

figure();
hold on;
subplot(2,2,3);
plot(U_bar,-imag(temp1));
xlabel('U_R');
ylabel('Growth ratio');
title('(c)');
subplot(2,2,4);
plot(U_bar,abs(real(temp1)));
xlabel('U_R');
ylabel('Frequency');
title('(d)');

temp2 = [EigValSort(3,2:9),EigValSort(1,10:22)] ;
subplot(2,2,1);
hold on;
plot(U_bar(2:22),-imag(temp2));
results = csvread('results_stability.csv');
plot(results(:,1),results(:,2),'xb');
plot(13.53,0,'s');
xlabel('U_R');
ylabel('Growth ratio');
title('(a)');
subplot(2,2,2);
hold on;
plot(U_bar(2:22),abs(real(temp2)));
results2 = csvread('results_frequency.csv');
plot(results2(:,1),results2(:,2),'xb');
plot(13.53,3.5,'s');
xlabel('U_R');
ylabel('Frequency');
title('(b)');

%%

[triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');
[triVecF,triIDsF,dummyF,triSupF] = mesh_reader('./meshes/fluidMesh.msh');

figure();
hold on;
trimesh(triIDsS, triVecS(:,1), triVecS(:,2),'Color','r');
trimesh(triIDsF, triVecF(:,1), triVecF(:,2),'Color','k');
grid off
axis equal

xlim([1.7 2.3])
ylim([-0.15 0.15])

%%
system(['FreeFem++ ./NS_base.edp -nH ', num2str(nH) ]);
%%
% Finding the first mode with P2-P1 elements.
% Because P1-P1 is not converged with the current mesh, we assume the value
% found was 2x larger than the real value.

% f = waitbar(1/(steps+1),strcat('Step ',num2str(1),' out of  ',num2str(steps+1)));

% Finding the first coupled mode with P2 elements, no flow
[A,B,As,Bs,sizeFu,sizeFp,sizeFl,sizeSu,sizeSe, FixedFu, FixedFp, FixedFl, FixedSu, FixedSe] = runCoupledP2P1(nH,Re,E,nus,rhos); 

%% Test solid part only

[T, Asb] = spbalance(As);

Bsb = T\Bs*T;

opts.tol = 1e-15;
opts.maxit = 300;
[EigVectTempSPartial,EigValTempS] = eigs(Asb,-1i.*Bsb,25,0,opts);
EigValTempS = diag(EigValTempS);

EigVectTempS = [sparse(sizeFu+sizeFp+sizeFl,25);EigVectTempSPartial];

dispEigVect(EigVectTempS(:,1),sizeFu,sizeFp,sizeFl,sizeSu,sizeSe, FixedFu, FixedFp, FixedFl, FixedSu, FixedSe,strcat("Mode 1 - ",num2str(EigValTempS(1))," - U: 0"),2);


for i=1:length(EigValTempS)
    EigNormS(i) = norm(As*EigVectTempSPartial(:,i)-1i.*Bs*EigVectTempSPartial(:,i)*EigValTempS(i));
end

%%
% 
% [T, Ab] = spbalance(A);
% Bb = T\B*T;

opts.tol = 1e-20;
opts.maxit = 500;
sigma = -0.5 + 3*i;
% sigma = 0 + 10.0*i;
[EigVectTemp,EigValTemp] = eigs(A,-1.0i.*B,50,0);
% [EigVectTemp,EigValTemp] = eigs(Ab,Bb,25,sigma,opts);
EigValTemp = diag(EigValTemp);

% EigVectTemp = T * EigVectTemp;

dispEigVect(EigVectTemp(:,1),sizeFu,sizeFp,sizeFl,sizeSu,sizeSe, FixedFu, FixedFp, FixedFl, FixedSu, FixedSe,strcat("Mode 1 - ",num2str(EigValTemp(1))," - U: 0"),2);

for i=1:length(EigValTemp)
    EigNorm(i) = norm(A*EigVectTemp(:,i)+1i.*B*EigVectTemp(:,i)*EigValTemp(i));
end


%%
% Just for quick testing
% opts.tol = 1e-45;
[T, Ab] = spbalance(A);
Bb = T\B*T;

opts.tol = 1e-15;
opts.maxit = 300;
% sigma = -0.50 + 4.00*i;
sigma = 0;
[EigVectTemp,EigValTemp] = eigs(Ab,-1i*Bb,25,sigma,opts);
% [EigVectTemp,EigValTemp] = eigs(A,B,25,sigma,opts);
EigValTemp = diag(EigValTemp);

EigVectTemp = T * EigVectTemp;

dispEigVect(EigVectTemp(:,1),sizeFu,sizeFp,sizeFl,sizeSu,sizeSe, FixedFu, FixedFp, FixedFl, FixedSu, FixedSe,strcat("Mode 1 - ",num2str(EigValTemp(1))," - U: 0"),2);

for i=1:length(EigValTemp)
    EigNorm(i) = norm(A*EigVectTemp(:,i)-B*EigVectTemp(:,i)*EigValTemp(i));
end


%%
results = csvread('results_new.csv')

figure(); hold on;

plot(-imag(EigValTemp(1:30)),real(EigValTemp(1:30)),'*k');
plot(results(:,1),results(:,2),'+b');
% plot(-imag(EigVal(:,1)),real(EigVal(:,1)),'*');
% plot(-imag(EigVal(:,2)),real(EigVal(:,2)),'*');
% plot(-imag(EigVal(:,3)),real(EigVal(:,3)),'*');
% plot(-imag(EigVal(:,4)),real(EigVal(:,4)),'*');
legend('Navier-Stokes','Pfister (2019)')
xlabel('Growth rate');
ylabel('Frequency');
xlim([-3.2 1]);

% plot(-imag(EigValTemp(1:12)),real(EigValTemp(1:12)),'or','MarkerSize',12);
% plot(-imag(EigValTemp(15:19)),real(EigValTemp(15:19)),'or','MarkerSize',12);
% plot(-imag(EigValTemp(13:14)),real(EigValTemp(13:14)),'og','MarkerSize',12);
% plot(-imag(EigValTemp(20:21)),real(EigValTemp(20:21)),'og','MarkerSize',12);

%%
% Guesses for the first few modes
freq(1) = [2.7+1.5i] * (sqrt(D/(rhos(1)*hp))/l^2);
freq(2) = [18+1.1i] * (sqrt(D/(rhos(1)*hp))/l^2);
% freq(2) = [0.5027+4.917i] * (sqrt(D/(rho_s(1)*hp))/l^2);
% freq(3) = [0.0123+30.889i] * (sqrt(D/(rho_s(1)*hp))/l^2);

for i = 1:length(freq)
    [EigVectTemp,EigValTemp] = eigs(A,B,1,freq(i));
    EigValTemp = diag(EigValTemp);

    EigVal(1,1,i) = EigValTemp(1);
    EigVect(:,1,1,i) = EigVectTemp(:,1);

    freq(i) = EigValTemp(1);
    freqADIM(i) = freq(i)/(sqrt(D/(rhos(1)*hp))/l^2);

    dispEigVect(EigVectTemp(:,1),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigValTemp(1))," - U: 0"),2);
end

% for i=1:length(EigValTemp)
%     [Ef(i),Es(i)] = getModeEnergy (sizeS,sizeFv,FixedS, FixedFv,EigVectTemp(:,i));
% end

 

%%

for i = 1:30
    [Ef,Es] = getModeEnergy (EigVectTemp(:,i),sizeFu,sizeFp,sizeFl,sizeSu,sizeSe, FixedFu, FixedFp, FixedFl, FixedSu, FixedSe)
    ratio(i) = Es/Ef;
end

%%

figure();
hold on;
xlabel('Mode number');
ylabel('E_s/E_f')

bar(1:12,ratio(1:12),'r');
bar(15:19,ratio(15:19),'r');
bar(22:30,ratio(22:30),'r');

bar(13:14,ratio(13:14),'g');
bar(20:21,ratio(20:21),'g');















