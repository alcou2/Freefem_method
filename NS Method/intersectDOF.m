function [Cx,Cy,indC] = intersectDOF(Ax,Ay,Bx,By)

indC = [];

if length(Ax) < length(Bx)
    for i = 1:length(Ax)
        for j = 1:length(Bx)
            if(abs(Ax(i)-Bx(j)) < 1e-6 && abs(Ay(i)-By(j))< 1e-6)
                indC = [indC j];
                Cx(i) = Ax(i);
                Cy(i) = Ay(i);
            end
        end
    end
end

Cx = Cx';
Cy = Cy';
indC = indC';