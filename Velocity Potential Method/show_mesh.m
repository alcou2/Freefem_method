[triVecS,triIDsS,dummyS,triSupS] = mesh_reader('./meshes/solidMesh.msh');
[triVecF,triIDsF,dummyF,triSupF] = mesh_reader('./meshes/fluidMesh.msh');

figure();
subplot(2,1,1);
hold on;
trimesh(triIDsS, triVecS(:,1), triVecS(:,2),'Color','r');
trimesh(triIDsF, triVecF(:,1), triVecF(:,2),'Color','k');
grid off
axis equal
xlim([-5 6])
ylim([-1 1])
axis off
% title('(a)')


subplot(2,1,2);
hold on;
trimesh(triIDsS, triVecS(:,1), triVecS(:,2),'Color','r');
trimesh(triIDsF, triVecF(:,1), triVecF(:,2),'Color','k');
grid off
axis equal
xlim([0.9 1.1])
ylim([-0.1 0.1])
axis off
% title('(b)')