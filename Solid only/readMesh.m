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

 
fprintf('Done reading mesh\n');
 
 
 
 






