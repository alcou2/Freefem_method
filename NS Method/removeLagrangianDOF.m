function [Asf,Af4,Af5,Bf,Afs,T,TestL,zeroL] = removeLagrangianDOF(Asf,Af4,Af5,Bf,Afs,T,TestL)


dofRowsToRemove = find(~any(TestL(:,:),1));  %rows
dofColumnsToRemove = find(~any(TestL(:,:),2));  %columns

Asf(:,dofColumnsToRemove) = [];
Af4(:,dofColumnsToRemove) = [];

Afs(dofRowsToRemove,:) = [];
Af5(dofRowsToRemove,:) = [];
T(dofRowsToRemove,:) = [];

TestL(dofRowsToRemove,:) = [];
TestL(:,dofColumnsToRemove) = [];

zeroL = dofRowsToRemove;

% dofRowsToRemove = find(~all(A==0,2);
% 
% dofRowsToRemove = ~any(Af(dofFv+dofFp:end,:),1);
% dofColsToRemove = ~any(Af(:,dofFv+dofFp:end),2);  %columns


end