clear all; clc;

l=1; % Plate length (m)
hp = 0.05; % Plate thinckness (m)

E = 13766; % Plate Youngs modulus
nus = 0.4; %Plate poisson ratio
rhos = 25; %Plate density


Re = [100]; % Reynolds number for which to evaluate
% Re = [10,50,100,250,500,1000,2500,5000];
Re_available = [10,20,30,40,50,60,70,80,90,100,125,175,200,225,250,300,350,400,450,500,600,700,800,900,1000,1250,1500,1750,2000,2250,2500,3000,3500,4000,4500,5000];
% Re = [10,30];


% U_barMin = 1;
% U_barMax = 15; % Max reduced speed
% steps = 56; % Number of steps in the speed continuation method
% U_bar = U_barMin:(U_barMax-U_barMin)/steps:U_barMax;
% E = (12*(1-nus^2)*rho_s*l^2*U^2)./(hp^2.*(U_bar.^2)); %Plate Young's modulus
% D = E.*(hp^3/(12*(1-nus^2)));

%mu = rho_f*l/(rho_s*hp);
%c = h/(2*l);

mhS = 1.0; % Solid domain mesh parameter
mhF = 8.0; % Fluid domain mesh parameter

MacCrit = 0.9;


%% Test


clear all; clc;

l=1; % Plate length (m)
rho_s = 20000; %Plate density (kg/m3)
hp = 0.05; % Plate thinckness (m)
nus = 0.3; %Plate poisson ratio

h = 2; % Fluid domain height
rho_f = 1000; % Fluid density
U = 1; % Fluid speed
%Re = [10,50,100,500,1000]; % Reynolds number for which to evaluate
Re = 10;
Re_available = [10,20,30,40,50,60,70,80,90,100,125,175,200,225,250,300,350,400,450,500,600,700,800,900,1000];
% Re = [10,30];
nuf = (U*rho_f*l)./Re; % Fluid dynamic viscosity  
% nuf = 1; % Fluid dynamic viscosity  
% U = Re*nuf/(rho_f*l); % Fluid speed

% U_barMin = 1;
% U_barMax = 3; % Max reduced speed
% steps = 8; % Number of steps in the speed continuation method
% U_bar = U_barMin:(U_barMax-U_barMin)/steps:U_barMax;
U_bar = 1;
E = (12*(1-nus^2)*rho_s*l^2*U^2)./(hp^2.*(U_bar.^2)); %Plate Young's modulus
D = E.*(hp^3/(12*(1-nus^2)));

K_b = D/(rho_f*U^2*l^2);

mu = rho_f*l/(rho_s*hp);
c = h/(2*l);

mhS = 1.0; % Solid domain mesh parameter
mhF = 8.0; % Fluid domain mesh parameter

MacCrit = 0.9;




%%
% Finding the first coupled mode using P1-P1 elements and using the
% idFSImodes method

% % Run Freefem script and get the matrices of the reduced order problem
% [A,B,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupledP1(mhS,mhF,U,nuf(1),rho_f,E(1),nus);
% 
% [EigVectTemp,EigValTemp] = eigs(A,B,1000,'SM');
% % [EigVectTemp,EigValTemp] = eigs(A,B,1,freq1);
% EigValTemp = diag(EigValTemp); 
% 
% [triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');
% [triVecF,triIDsF,dummyF,triSupF] = mesh_reader('./meshes/fluidMesh.msh');
% 
% [EigVectNoFlow,EigValNoFlow] = idFSImodes(EigVectTemp,EigValTemp,sizeS,sizeFv,sizeFp,triVecS,triVecF,FixedS,FixedFv);
% 
% freq1 = EigValNoFlow(1);
% if real(freq1) < 0
%     freq1 = freq1-2*real(freq1);
% end
% 
% % freq2 = EigValNoFlow(3);
% % if real(freq2) < 0
% %     freq2 = freq2-2*real(freq2);
% % end
% 
% % guess = (0.060839+0.017455i)*Re/10;
% % [EigVectTemp,EigValTemp] = eigs(A,B,1,guess);
% % EigVectNoFlow = EigVectTemp;
% % freq1 = EigValTemp(1);
% 
% fprintf('Identified mode with P1 elements: \n')
% fprintf('Mode 1: %f+%fi \n',real(freq1),imag(freq1));
% % fprintf('Mode 2: %f+%fi \n',real(freq2),imag(freq2));
% 
% NbModes = 1;
% 
% % Displaying potential coupled modes
% dispEigVect(EigVectNoFlow(:,1),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(freq1)," - U: 0"),1);
% % dispEigVect(EigVectNoFlow(:,3),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 2 - ",num2str(freq2)," - U: 0"),1);

%%
% i = 69;
% dispEigVect(EigVectTemp(:,i),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigValTemp(i))," - U: 0"),1);

%%

% Finding the first mode with P2-P1 elements.
% Because P1-P1 is not converged with the current mesh, we assume the value
% found was 2x larger than the real value.

% f = waitbar(1/(steps+1),strcat('Step ',num2str(1),' out of  ',num2str(steps+1)));

% Finding the first coupled mode with P2 elements, no flow
[A,B,As,Bs,sizeFu,sizeFp,sizeFl,sizeSu,sizeSe, FixedFu, FixedFp, FixedFl, FixedSu, FixedSe] = runCoupledP2P1(mhS,mhF,Re,E,nus,rhos); 

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
[EigVectTemp,EigValTemp] = eigs(A,-1.0i.*B,25,0);
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
results = csvread('results.csv')

figure(); hold on;

plot(results(:,1),results(:,2),'+');
plot(-imag(EigValTemp),real(EigValTemp),'*');
legend('method','Fernandez')
% plot(real(EigValTemp),imag(EigValTemp),'*');
%%
% Guesses for the first few modes
freq(1) = [2.7+1.5i] * (sqrt(D/(rho_s(1)*hp))/l^2);
freq(2) = [18+1.1i] * (sqrt(D/(rho_s(1)*hp))/l^2);
% freq(2) = [0.5027+4.917i] * (sqrt(D/(rho_s(1)*hp))/l^2);
% freq(3) = [0.0123+30.889i] * (sqrt(D/(rho_s(1)*hp))/l^2);

for i = 1:length(freq)
    [EigVectTemp,EigValTemp] = eigs(A,B,1,freq(i));
    EigValTemp = diag(EigValTemp);

    EigVal(1,1,i) = EigValTemp(1);
    EigVect(:,1,1,i) = EigVectTemp(:,1);

    freq(i) = EigValTemp(1);
    freqADIM(i) = freq(i)/(sqrt(D/(rho_s(1)*hp))/l^2);

    dispEigVect(EigVectTemp(:,1),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigValTemp(1))," - U: 0"),2);
end

% for i=1:length(EigValTemp)
%     [Ef(i),Es(i)] = getModeEnergy (sizeS,sizeFv,FixedS, FixedFv,EigVectTemp(:,i));
% end

 
%%

% Continuation method on Re

f = waitbar(0,'Continuation on Re : Starting');

for i = 2:length(Re)

% 
%     [A,B,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupledP2P1(mhS,mhF,U,nuf(i),rho_f,E(1),nus); 
%     
%     freq1 = freq1ADIM * (sqrt(D(2)/(rho_s*hp))/l^2);
%     [EigVectTemp,EigValTemp] = eigs(A,B,5,freq1);
%     EigValTemp = diag(EigValTemp);
%     
%     for j = 1:length(EigValTemp)
%         mac(j) = MAC(EigVectTemp(:,j),EigVect1(:,i-1));
%     end
%     
%     Mac = mac(1);
%     
    
    x = (log(Re(i))/log(max(Re))-log(min(Re))/log(max(Re)))/(log(min(Re))/log(max(Re)));
    waitbar(x,f,strcat('Continuation on Re : ',num2str(Re(i))));

    [Atemp,Btemp,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupledP2P1(mhS,mhF,U,nuf(i),rho_f,E,nus,rho_s(1)); 

    nuf_temp = nuf(i);
    nuf_past = nuf(i-1);
    EigVectPast = squeeze(EigVect(:,1,i-1,:));
    MacMin = 0;
    
    divMax = 4;
    divNum = 0;
    
    while MacMin < MacCrit || nuf_temp > nuf(i)

        if divNum <= divMax 
            x = (log((U*rho_f*l)/nuf_temp)/log(max(Re))-log(min(Re))/log(max(Re)))/(log(min(Re))/log(max(Re)));
            waitbar(x,f,strcat('Continuation on Re : ',num2str((U*rho_f*l)/nuf_temp),' - MAC: ',num2str(MacMin,'%.3f')));

            if nuf_temp == nuf(i)    
                A = Atemp;
                B = Btemp;
            else
                [A,B,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupledP2P1(mhS,mhF,U,nuf_temp,rho_f,E,nus,rho_s(1)); 
            end

            for k = 1:length(freq)
                freq(k) = freqADIM(k) * (sqrt(D/(rho_s(1)*hp))/l^2);
                [EigVectTemp(:,1,k),EigValTemp(k)] = eigs(A,B,1,freq(k));

    %             dispEigVect(EigVectTemp(:,1),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigValTemp(1))," - U: 0"),2);

                MacT = MAC(EigVectTemp(:,1,k),EigVectPast(:,k));
                MacS = MAC(EigVectTemp(1:sizeS,1,k),EigVectPast(1:sizeS,k));

                Mac(k) = min(MacT,MacS);
            end

            MacMin = min(Mac);

            if MacMin < MacCrit 
                nuf_temp = (nuf_past+nuf_temp)/2;
                divNum = divNum + 1;
            end
            if MacMin > MacCrit && nuf_temp ~= nuf(i)
                divNum = 0;
                nuf_past = nuf_temp;
                nuf_temp = nuf(i);
                EigVectPast(:,:) = squeeze(EigVectTemp(:,1,:));
                freq(:) = EigValTemp(:);
                freqADIM(:) = freq(:)/(sqrt(D/(rho_s(1)*hp))/l^2);
                MacMin = 0;
            end
        end
        if divNum > divMax
            [A,B,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupledP2P1(mhS,mhF,U,nuf_temp,rho_f,E,nus,rho_s(1)); 
            
            for k = 1:length(freq)
                freq(k) = freqADIM(k) * (sqrt(D/(rho_s(1)*hp))/l^2);
                [EigVectTempMult(:,1,k,:),EigValTempMult(k,:)] = eigs(A,B,10,freq(k));
                for j = 1:size(EigValTempMult,2)
                    MacT = MAC(EigVectTempMult(:,1,k,j),EigVectPast(:,k));
                    MacS = MAC(EigVectTempMult(1:sizeS,1,k,j),EigVectPast(1:sizeS,k));
                    MacMult(k,j) = min(MacT,MacS);
                end
                [Mac(k),index] = max(MacMult(k,:));
                EigValTemp = EigValTempMult(:,index);
                EigVectTemp = EigVectTempMult(:,1,k,index);
            end
            [MacMin,index] = min(Mac);
            
            if MacMin > MacCrit
                divNum = 0;
                nuf_past = nuf_temp;
                nuf_temp = nuf(i);
                EigVectPast(:,:) = squeeze(EigVectTemp(:,1,:));
                freq(:) = EigValTemp(:);
                freqADIM(:) = freq(:)/(sqrt(D/(rho_s(1)*hp))/l^2);
                MacMin = 0;
            else
                disp(['Current Re: ',num2str((U*rho_f*l)/nuf_temp)]);
                disp(['Last tried Re: ',num2str((U*rho_f*l)/nuf_past)]);
                disp(['Last resolved Re: ',num2str(Re(i-1))]);
                
                for I=1:size(EigVectTempMult,4)
                    for J=1:size(EigVectTempMult,4)
                        mac(I,J)=MAC(EigVectTempMult(:,1,index,I),EigVectTempMult(:,1,index,J));
                    end
                end
                % plot mac matrix
                figure
                bar3(mac)
                title('MAC')
                
                error(['Continuation method failed at Re: ',num2str(Re(i))]);
                disp
            end
            
        end
%         for j = 1:length(EigValTemp)
%             mac(j) = MAC(EigVectTemp(:,j),EigVect1(:,i-1));
%         end
        
    end
    
    for k = 1:length(freq)
        freq(k) = EigValTemp(k);
        freqADIM(k) = freq(k)/(sqrt(D/(rho_s(1)*hp))/l^2);
        EigVal(1,i,k) = EigValTemp(k);
        EigVect(:,1,i,k) = EigVectTemp(:,1,k);
    end
    
end

close(f);

%%

% Continuation method on UR for every Re

f = waitbar(0,'Continuation on Mu : Starting');

for i = 1:length(Re)

    
    [matrices,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = getMatrices(mhS,mhF,U,nuf(i),rho_f,E(1),nus);
    
    
    
    for j = 2:length(E)
        
        E_temp = E(j);
        E_past = E(j-1);
        EigVectPast = squeeze(EigVect(:,(j-1),i,:));
        EigValPast = squeeze(EigVal(j-1,i,:));
        U_barTemp = U_bar(j);
        U_barPast = U_bar(j-1);
        if j > 2
            EigValPast2 = squeeze(EigVal(j-2,i,:));
            U_barPast2 = U_bar(j-2);
        else
            EigValPast2 = EigValPast;
            U_barPast2 = U_barPast;
        end
        MacMin = 0;

        
        while MacMin < MacCrit || E_temp > E(j)
            
            
            x = ((U_barTemp-U_barMin)/(U_barMax-U_barMin)/length(Re))+((i-1)/length(Re));
            waitbar(x,f,strcat('Continuation on UR : ',num2str(U_barTemp),' - Re : ',num2str(Re(i)),' - MAC: ',num2str(MacMin,'%.3f')));
            
            [A,B,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupledP2P1step(mhS,mhF,U,nuf(i),rho_f,E_temp,nus,matrices,FixedS,FixedFv,FixedFp);
            
            
            for k = 1:length(freq)
                
                if j == 2
                    EigValGuess = EigValPast(k);
                else
                    EigValGuess = ((EigValPast(k)-EigValPast2(k))/(U_barPast-U_barPast2))*(U_barTemp-U_barPast)+EigValPast(k);
                end
                
                [EigVectTemp(:,1,k),EigValTemp(k)] = eigs(A,B,1,EigValGuess);

                MacT = MAC(EigVectTemp(:,1,k),EigVectPast(:,k));
                MacS = MAC(EigVectTemp(1:sizeS,1,k),EigVectPast(1:sizeS,k));

                Mac(k) = min(MacT,MacS);               
                
            end
            
            MacMin = min(Mac);
            
            if MacMin < MacCrit
            E_temp = (E_past+E_temp)/2;
            U_barTemp = (U_barPast+U_barTemp)/2;
            end
            if MacMin > MacCrit && E_temp ~= E(j)
                E_past = E_temp;
                E_temp = E(j);
                U_barPast2 = U_barPast;
                U_barPast = U_barTemp;
                U_barTemp = U_bar(j);
                
                EigVectPast = squeeze(EigVectTemp(:,1,:));
                EigValPast = EigValTemp(:);
                MacMin = 0;
            end
            
        end
        
        for k = 1:length(freq)
            EigVal(j,i,k) = EigValTemp(k);
            EigVect(:,j,i,k) = EigVectTemp(:,1,k);
        end
        
    end
    
end 

close(f);

%%



% for i = 1:length(freq)
%     
%     figure();
%     hold on;
%     
%     for j = 1:length(Re)
%         
%         EigVal_bar = EigVal(:,j,i)'./(sqrt(D/(rho_s*hp))/l^2);
%         plot(U_bar,real(EigVal_bar),'-');
%         plot(U_bar,-imag(EigVal_bar),'.');
%         
%     end
% 
%     legend('Re: 10','','Re: 30','','Re: 35','','Re: 40','','Re: 50','','Re: 100','');
% 
% end



%% Plot MAC matrix

% for I=1:size(EigVect,2)
%     for J=1:size(EigVect,2)
%         mac(I,J)=MAC(EigVect(:,I,3,1),EigVect(:,J,3,1));
%     end
% end
% % plot mac matrix
% figure
% bar3(mac)
% title('MAC')


%% Manually display modes

i = 1 % Mu index
j = 3 % Re index
k = 1 % Freq index

dispEigVect(EigVect(:,i,j,k),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigVal(i,j,k))," - U: 0"),2);


%%

% filename = 'Animated_Re_100.gif';
% 
% j = 4 % Re index
% 
% for i = 1:length(E)
%     
%     f = dispEigVect(EigVect1(:,i,j),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigVal1(i,j))," - U: 0"),2);
%     
%     frame = getframe(f); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     
%     if i == 1 
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%     else 
%         imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%     end 
%     
% end



%%

% % Finding the first coupled mode with P2 elements, with flow
% % waitbar(2/(steps+1),f,strcat('Step ',num2str(2),' out of  ',num2str(steps+1)));
% 
% [A,B,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupledP2P1(mhS,mhF,U,nuf,rho_f,E(2),nus); 
% 
% [EigVectTemp,EigValTemp] = eigs(A,B,5,freq1); 
% EigValTemp = diag(EigValTemp);
% 
% EigVal1(2) = EigValTemp(1);
% EigVect1(:,2) = EigVectTemp(:,1);
% 
% freq1 = EigValTemp(1);
% 
% dispEigVect(EigVectTemp(:,1),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigValTemp(1))," - U: 0"),2);
% 
% [Ef,Es] = getModeEnergy (sizeS,sizeFv,FixedS, FixedFv,EigVectTemp(:,1));
% 
% 
% %%
% 
% 
% %     [A,B,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupledP2P1(0,mhS,mhF,nuf_f); 
% %     [EigVectTemp,EigValTemp] = eigs(A,B,1,EigValTemp(1));
% %     EigValTemp = diag(EigValTemp)
% % 
% %     EigVal1(1) = EigValTemp(1);
% %     EigVect1(:,1) = EigVectTemp(:,1);
% % 
% %     freq1 = EigValTemp(1);
% % 
% %     dispEigVect(EigVectTemp(:,1),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigValTemp(1))," - U: 0"),2);
% % 
% %     [Ef,Es] = getModeEnergy (sizeS,sizeFv,FixedS, FixedFv,EigVectTemp(:,1))
% 
% 
% %%
% 
% % i = 1;
% % 
% % dispEigVect(EigVectTemp(:,i),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigValTemp(i))," - U: 0"),2);
% % 
% % [Ef,Es] = getModeEnergy (sizeS,sizeFv,FixedS, FixedFv,EigVectTemp(:,i))
% 
% 
% 
% %% Progressivly increasing the base flow  speed
% 
% 
% [matrices,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = getMatrices(mhS,mhF,U,nuf,rho_f,E(3),nus);
% 
% %%
% 
% for i = 3:length(E) 
%       
%     waitbar(i/(steps+1),f,strcat('Step ',num2str(i),' out of  ',num2str(steps+1)));
%     
%     [A,B,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupledP2P1step(mhS,mhF,U,nuf,rho_f,E(i),nus,matrices,FixedS,FixedFv,FixedFp);  
% %     [A,B,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupledP2P1(mhS,mhF,U,nuf,rho_f,E(i),nus);
%     
%     if length(EigVal1) == 1
%         EigGuess = EigVal1(i-1); 
%     else
%         EigGuess = ((EigVal1(i-1)-EigVal1(i-2))/(U_bar(i-1)-U_bar(i-2)))*(U_bar(i)-U_bar(i-1))+EigVal1(i-1);
%     end
%     
%     [EigVectTemp,EigValTemp] = eigs(A,B,1,EigGuess);
% %     [EigVectTemp,EigValTemp] = eigs(A,B,1,EigVal1(i-1));
%     EigValTemp = diag(EigValTemp);
%     EigVal1(i) = EigValTemp(1);
%     EigVect1(:,i) = EigVectTemp(:,1);        
%     
% %     dispEigVect(EigVect1(:,i),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigVal1(i))," - U: ",num2str(U(i))),2);
%     
% end  
% 
% %%
% 
% % U_bar = 0:U_barMax/steps:U_barMax;
% EigVal_bar = EigVal1./(sqrt(D/(rho_s*hp))/l^2);
% 
% 
% figure();
% hold on;
% plot(U_bar,real(EigVal_bar));
% plot(U_bar,imag(EigVal_bar));
% legend('Real','Imag');
% 
% %%
% i =  9;
% dispEigVect(EigVect1(:,i),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigVal1(i))," - U: ",num2str(U_bar(i))),2);
% 
% % [Ef,Es] = getModeEnergy (sizeS,sizeFv,FixedS, FixedFv,EigVect1(:,i))
% 
% % % Testing
% % 
% % i = 1;
% % U = 0.0;
% % 
% % [A,B,sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp] = runCoupledP2P1step(U(i),mhS,mhF,matrices,FixedS,FixedFv,FixedFp);  
% % 
% % [EigVectTemp,EigValTemp] = eigs(A,B,1,EigVal1(i));
% % EigValTemp = diag(EigValTemp);
% % 
% % i=1;
% % dispEigVect(EigVectTemp(:,i),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigValTemp(i))," - U: ",num2str(U(i))),2);
% % 
% % [Ef,Es] = getModeEnergy (sizeS,sizeFv,FixedS, FixedFv,EigVect1(:,i))
% % 
% % figure();
% % plot(real(EigValTemp(1)),imag(EigValTemp(1)),'x');
% % hold on;
% % plot(real(EigValTemp(2:end)),imag(EigValTemp(2:end)),'o');
% % 
% % %
% % i=1;
% % dispEigVect(EigVect1(:,i),sizeS,sizeFv,sizeFp, FixedS, FixedFv,FixedFp,strcat("Mode 1 - ",num2str(EigVal1(i))," - U: ",num2str(U(i))),2);
% 
% % i = 1;
% 
% % [Ef,Es] = getModeEnergy (sizeS,sizeFv,FixedS, FixedFv,EigVect1(:,i))
% 
% %%
% 
% for i = 1:length(U_bar)
%     [Ef,Es] = getModeEnergy (sizeS,sizeFv,FixedS, FixedFv,EigVect1(:,i))
%     ratio(i) = Es/Ef;
% end
% 
% figure();
% plot(U_bar,ratio);
% 
% close(f);
% %%
% 





