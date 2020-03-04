function f = dispEigVect(EigVect,sizeFu,sizeFp,sizeFl,sizeSu,sizeSe, FixedFu, FixedFp, FixedFl, FixedSu, FixedSe,titleStr,order)

% Split eigenvectors into solid displacement, fluid velocity and pressure
% Adds zeros at the dirichlet boundary conditons

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


% Import meshes

[XPosS,YPosS,XPosF,YPosF] = readMesh('./meshes/xs.dat','./meshes/ys.dat','./meshes/xf.dat','./meshes/yf.dat');

[triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');
[triVecF,triIDsF,dummyF,triSupF] = mesh_reader('./meshes/fluidMesh.msh');

% f = figure();
% hold on; axis equal;
% 
% %========================================
% %Show isopressure lines
% % trisurfc( triVecF(:,1), triVecF(:,2),real(fluidEigP),200);
% 
% 
% %========================================
% %Show fluid flow vector field
% % scaling = 0.1/max(sqrt(real(VxEig).^2+real(VyEig).^2));
% quiver(XPosF,YPosF,real(VxEig),real(VyEig));
% % quiver(XPosF,YPosF,scaling*real(VxEig),scaling*real(VyEig),0);
% 
% %========================================
% %Show streamlines
% 
% N = 25; % number of seed locations
% % xmin = min(triVecF(:,1));
% % xmax = max(triVecF(:,1));
% % ymin = min(triVecF(:,2));
% % ymax = max(triVecF(:,2));
% % xmin = -1;
% % xmax = 2;
% % ymin = -1;
% % ymax = 1;
% % xstart = xmin + (-xmin+xmax)*rand(N,1); 
% % ystart = ymin + (-ymin+ymax)*rand(N,1);
% % xstart = [0:0.05:1, 0:0.05:1];
% % ystart = [ones([1,size(xstart,2)./2]).*0.03,ones([1,size(xstart,2)./2]).*-0.03];
% % xstart = [0:0.05:1];
% % ystart = [ones([1,size(xstart,2)]).*0.03];
% % 
% % FVx = scatteredInterpolant(XPosF,YPosF,real(VxEig),'linear','none');
% % FVy = scatteredInterpolant(XPosF,YPosF,real(VyEig),'linear','none');
% % 
% % [X2, Y2] = meshgrid(linspace(min(XPosF),max(XPosF),500),...
% %         linspace(min(YPosF),max(YPosF),100));
% %     
% % V2x = FVx(X2,Y2);
% % V2y = FVy(X2,Y2);
% % 
% % streamline(X2,Y2,V2x,V2y,xstart,ystart);
%  
% 
% %FlowP=TriStream(triIDsF,triVecF(:,1)',triVecF(:,2)',real(VxEig)',real(VyEig)',xstart,ystart);
% %PlotTriStream(FlowP,'b');
% 
% 
% 
% 
% 
%     
% %========================================
% % Show pressure colormap
% trisurf(triIDsF, triVecF(:,1), triVecF(:,2),real(fluidEigP)-1,'EdgeColor','none','FaceColor','interp');
% 
% % [Cx,Cy,indC]=intersectDOF(triVecF(:,1),triVecF(:,2),XPosF,YPosF);
% % fluidEigP = fluidEigP(indC);
% % trisurf(triIDsF, triVecF(:,1), triVecF(:,2),real(fluidEigP)-1,'EdgeColor','none','FaceColor','interp');
%     
%     
% %========================================
% % Show deformed solid mesh
%     
% [Cx,Cy,indC]=intersectDOF(triVecS(:,1),triVecS(:,2),XPosS,YPosS);
% solidEig1x = -solidEigE(indC,1);
% solidEig1y = -solidEigE(indC,2);
% 
% % scaling = 0.3/max(abs(real(solidEigE(:,2)))); % Deflection scaling factor
% scaling = 0./max(abs(real(solidEigE(:,2)))); % Deflection scaling factor
% if max(abs(real(solidEigE(:,1).*scaling))) > 0.2
%     scaling = 0.2/max(abs(real(solidEigE(:,1))));
% end
% trimesh(triIDsS, triVecS(:,1)+real(solidEig1x).*scaling, triVecS(:,2)+real(solidEig1y).*scaling,'Color','r');
% trimesh(triIDsS, triVecS(:,1)+imag(solidEig1x).*scaling, triVecS(:,2)+imag(solidEig1y).*scaling,'Color','b');
% 
% 
% 
% xlim([-3 10]);
% ylim([-1 1]);
% axis equal;
% 
% title(titleStr);


%%

f = figure();
hold on; axis equal;

%========================================
% Show streamwise velocity

FVx = scatteredInterpolant(XPosF,YPosF,imag(VxEig),'linear','none');

[X2, Y2] = meshgrid(linspace(min(XPosF),max(XPosF),1000),...
        linspace(min(YPosF),max(YPosF),1000));
    
V2x = FVx(X2,Y2);

s = pcolor(X2,Y2,V2x); 
s.FaceColor = 'interp';   
set(s, 'EdgeColor', 'none','FaceColor','interp');


%========================================
% Show deformed solid mesh
    
[Cx,Cy,indC]=intersectDOF(triVecS(:,1),triVecS(:,2),XPosS,YPosS);
solidEig1x = -solidEigE(indC,1);
solidEig1y = -solidEigE(indC,2);

% scaling = 0.3/max(abs(real(solidEigE(:,2)))); % Deflection scaling factor
scaling = 0.2/max(abs(imag(solidEigE(:,2)))); % Deflection scaling factor
if max(abs(imag(solidEigE(:,1).*scaling))) > 0.2
    scaling = 0.2/max(abs(imag(solidEigE(:,1))));
end
trimesh(triIDsS, triVecS(:,1)+imag(solidEig1x).*scaling, triVecS(:,2)+imag(solidEig1y).*scaling,'Color','r');
% trimesh(triIDsS, triVecS(:,1)+imag(solidEig1x).*scaling, triVecS(:,2)+imag(solidEig1y).*scaling,'Color','b');



xlim([-3 10]);
ylim([-1 1]);
% xlim([-0.5 5]);
% ylim([-0.5 0.5]);
axis equal;

title(titleStr);