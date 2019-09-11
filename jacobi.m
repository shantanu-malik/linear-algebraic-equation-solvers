function [x, R, n] = jacobi(A,B)
% Header: function [x, R, n] = sor(A,B)
% This function solves A*X = B using successive over-relaxaton iteration method.
% x = unknowns
% R = remainders
% n = no. of iterations

t = 0.000001;   % tolerance
s = length(B);

 w=1;
 while w <= s
    for i=1:s
    if sum(abs(A(i,:)))- abs(A(i,w)) < abs(A(i,w))
        A([w i],:)=A([i w],:);
        B([w i],:)=B([i w],:);
        w=w+1;
        break
    elseif sum(abs(A(i,:)))- abs(A(i,w))== abs(A(i,w))
        A([w i],:)=A([i w],:);
        B([w i],:)=B([i w],:);
        if sum(abs(A(i,:)))- abs(A(i,w))< abs(A(i,w))
        A([w i],:)=A([i w],:);
        B([w i],:)=B([i w],:);
        w=w+1;
        break
        end
        if i==s
            w=w+1;
        end

    elseif i==s
        w=w+1;
    else
        continue

    end

    end
 end
 
x = zeros(s,1);
R = x;
n = 0;

while (1)
    flag = 0;
    C = A*x;
    for i=1:s
        R(i) = B(i) - C(i);
        if((abs(R(i)) > t)&&(flag==0))
            flag = 1;
        end
    end
    n=n+1;
    if(flag==0)
        break;
    end
    for(i=1:s)
        x(i)=x(i)+R(i)/A(i,i);
    end
end
end