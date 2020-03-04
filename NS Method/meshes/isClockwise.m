% Fucntion to determine if the vertices of a polygon are in clockwise
% order. Same as Matlab function ispolycw without the need for a toolbox. 

function cw=isClockwise(x,y)

if ~isempty(x)
if nargin>=2
if size(x,1)~=1
x=x';
end
if size(y,1)~=1
y=y';
end
xx=[x x(1)];
yy=[y y(1)]; 
else
if size(x,2)~=2
x=x';
end
xx=[x(:,1); x(1,1)];
yy=[x(:,2); x(1,2)];
end
a=sum((xx(1:end-1)-xx(2:end)).*(yy(1:end-1)+yy(2:end)))/2;
else
a=0;
end

if a > 0
    cw = true;
else
    cw = false;
end