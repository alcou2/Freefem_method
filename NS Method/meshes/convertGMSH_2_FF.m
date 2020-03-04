clear all;

meshFile = 'fluidMesh2.mesh';
newMeshFile = 'fluidMeshFF2.msh';
% meshFile = 'solidMesh2.mesh';
% newMeshFile = 'solidMeshFF2.msh';

fid = fopen(meshFile, 'r');

line = '';

%Extract vertices information
while ~strcmp(line,' Vertices')
    line = fgetl(fid);
end

nbVertices = str2num(fgetl(fid));

Vertices = zeros(nbVertices,4);

for i = 1:nbVertices    
    Vertices(i,:) = sscanf(fgetl(fid),'%f',4);    
end

%Extract edges information
while ~strcmp(line,' Edges')
    line = fgetl(fid);
end

nbEdges = str2num(fgetl(fid));

Edges = zeros(nbEdges,3);

for i = 1:nbEdges   
    Edges(i,:) = sscanf(fgetl(fid),'%f',3);    
end

%Extract elements information
while ~strcmp(line,' Triangles')
    line = fgetl(fid);
end

nbElements = str2num(fgetl(fid));

Elements = zeros(nbElements,4);

for i = 1:nbElements    
    Elements(i,:) = sscanf(fgetl(fid),'%f',4);   
end

fclose(fid);


fid = fopen(newMeshFile,'wt');

fprintf(fid,'%i %i %i\n',nbVertices,nbElements,nbEdges);

for i = 1:nbVertices
    fprintf(fid,'%g %g %g\n',Vertices(i,[1,2,4]));
end

for i = 1:nbElements
    % Freefem wants the elements to be defined clockwise
    x = Vertices(Elements(i,1:3),1);
    y = Vertices(Elements(i,1:3),2);
    if isClockwise(x,y)
        fprintf(fid,'%i %i %i %i\n',Elements(i,[1,2,3,4]));
    else
        fprintf(fid,'%i %i %i %i\n',Elements(i,[1,3,2,4]));
    end
end

for i = 1:nbEdges
    fprintf(fid,'%i %i %i\n',Edges(i,:));
end
    
fclose(fid);