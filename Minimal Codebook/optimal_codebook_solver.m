function [minimal_codebook,indices]=optimal_codebook_solver(codebook)
%Takes in a codebook of genetic expression for cells (rows) and genes
%(columns) and reports a minimal subset such that no two row has identical
%gene expression.
%INPUT: codebook (n x d) matrix of n cells, d genes, binary
%OUTPUT: minimal_codebook (n x k) matrix of n cells, k genes, binary
%        indices (k x 1) matrix that stores the indices of the genes kept

%Copyright Molly B. Reilly, Cyril Cros, Erdem Varol, Eviatar Yemini, and Oliver Hobert 2020     

[~,d]=size(codebook);
bits=(1:size(codebook,2));

[A,b]=constraint_set(codebook,bits);
ub=ones(d,1);
lb=zeros(d,1);
Aeq=[];
beq=[];
f=ones(numel(bits),1);
intcon=(1:size(codebook,2));
[x,~] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub)   ;
indices=find(x>=0.99);
minimal_codebook=codebook(:,indices);

if constraint_satisfy(codebook,indices)
    disp('Constraints satisfied!')
end
end

function out=constraint_satisfy(codebook,bits)
% Check whether constraints are satisfied
out=(min(pdist(codebook(:,bits),'hamming'))>0); %Check whether all of the codebook entries are unique

end

function [A,b]=constraint_set(codebook,bits)
% Knapsack constraint setting
C=codebook(:,bits);
A=zeros(size(codebook,1)*(size(codebook,1)-1)/2,numel(bits));
for i=1:size(C,2)
    A(:,i)=-pdist(C(:,i),'hamming');
end
b=-ones(size(A,1),1);
end
