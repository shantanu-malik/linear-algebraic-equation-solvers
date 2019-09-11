function X = gaussElimination(A,B)
% Header: function X = gaussElimination(A,B)
% This function solves A*X = B using Gauss Elimination method.

AB=[A B];
%% Partial pivoting 
for i=1:size(B)
    ratios(i)=abs(A(i,1)/A(i,end));
end

[sorted, indexOrder]=sort(ratios);

ab=zeros(length(B),length(B)+1);
for i=1:size(B)
    ab(end-i+1,:)=AB(indexOrder(i),:);
end

%% Elimination
for J=1:size(B)-1
    for i=size(B):-1:J+1
        c=ab(J,J);d=ab(i,J);
        for j=1:size(B)+1
            ab(i,j)=ab(i,j)-(d/c)*ab(J,j);
        end
    end
end

a=ab(:,1:size(B));b=ab(:,end);

%% Calculating unknowns
X=zeros(length(B),1);
for i=size(B):-1:1
    x=a*X;
    X(i)=(b(i)-x(i))/a(i,i);
end
end