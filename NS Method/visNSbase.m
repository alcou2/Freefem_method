xffile = './NSresults/xf.dat';
xf = fopen(xffile, 'r');

XPosF = [];
fgetl(xf);

 while ~feof(xf)
    line = fgetl(xf) ; 
    content = textscan(line, '%f') ; 
    content = cell2mat(content);
    XPosF = [XPosF;content];
 end
 fclose(xf) ;
 
 
 
yffile = './NSresults/yf.dat';
yf = fopen(yffile, 'r');

YPosF = [];
fgetl(yf);

 while ~feof(yf)
    line = fgetl(yf) ; 
    content = textscan(line, '%f') ; 
    content = cell2mat(content);
    YPosF = [YPosF;content];
 end
 fclose(yf) ;
 
fprintf('Done reading mesh\n');



ufile = './NSresults/U10.dat';
u = fopen(ufile, 'r');

U = [];
fgetl(u);

 while ~feof(u)
    line = fgetl(u) ; 
    content = textscan(line, '%f') ; 
    content = cell2mat(content);
    U = [U;content];
 end
 fclose(u) ;
 
 
 
vfile = './NSresults/V10.dat';
v = fopen(vfile, 'r');

V = [];
fgetl(v);

 while ~feof(v)
    line = fgetl(v) ; 
    content = textscan(line, '%f') ; 
    content = cell2mat(content);
    V = [V;content];
 end
 fclose(v) ;
 
pfile = './NSresults/P210.dat';
p = fopen(pfile, 'r');

P = [];
fgetl(p);

 while ~feof(p)
    line = fgetl(p) ; 
    content = textscan(line, '%f') ; 
    content = cell2mat(content);
    P = [P;content];
 end
 fclose(p) ;
 
figure();
axis equal;
quiver(XPosF,YPosF,U,V);

figure();
axis equal;
Vel = sqrt(U.^2+V.^2);
F = scatteredInterpolant(XPosF,YPosF,Vel);
[xq,yq] = meshgrid(linspace(min(XPosF),max(XPosF),500),linspace(min(YPosF),max(YPosF),500));
cq = F(xq,yq);
h = pcolor(xq,yq,cq);
h.EdgeColor = 'none';
axis equal;

figure();
axis equal;
F = scatteredInterpolant(XPosF,YPosF,P);
[xq,yq] = meshgrid(linspace(min(XPosF),max(XPosF),500),linspace(min(YPosF),max(YPosF),500));
cq = F(xq,yq);
h = pcolor(xq,yq,cq);
h.EdgeColor = 'none';
axis equal;

% z = sqrt(U.^2+V.^2);
% minx = min(XPosF);
% maxx = max(XPosF);
% miny = min(YPosF);
% maxy = max(YPosF);
% meanValue = mean(z);
% heatMapImage = meanValue  * ones(100, 100);
% for k = 1 : length(XPosF)
%   column = round( (XPosF(k) - minx) * 100 / (maxx-minx) ) + 1;;  
%   row = round( (YPosF(k) - miny) * 100 / (maxy-miny) ) + 1;
%   heatMapImage(row, column) = z(k);
% end
% imshow(heatMapImage, []);
% colormap('hot');
% colorbar;


figure();
ystart = [-2.95:0.05:2.95]; 
xstart = [ones([1,size(ystart,2)])*-2];

FVx = scatteredInterpolant(XPosF,YPosF,U,'linear','none');
FVy = scatteredInterpolant(XPosF,YPosF,V,'linear','none');

[X2, Y2] = meshgrid(linspace(min(XPosF),max(XPosF),500),...
        linspace(min(YPosF),max(YPosF),100));
    
V2x = FVx(X2,Y2);
V2y = FVy(X2,Y2);

streamline(X2,Y2,V2x,V2y,xstart,ystart);
axis equal;