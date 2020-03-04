function [XPosS,YPosS,XPosF,YPosF] = readMesh(xsfile,ysfile,xffile,yffile)


xs = fopen(xsfile, 'r');

XPosS = [];
fgetl(xs);

 while ~feof(xs)
    line = fgetl(xs) ; 
    content = textscan(line, '%f') ; 
    content = cell2mat(content);
    XPosS = [XPosS;content];
 end
 fclose(xs) ;

 
ys = fopen(ysfile, 'r');

YPosS = [];
fgetl(ys);

 while ~feof(ys)
    line = fgetl(ys) ; 
    content = textscan(line, '%f') ; 
    content = cell2mat(content);
    YPosS = [YPosS;content];
 end
 fclose(ys) ;

 
 
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
 
 
 
 






