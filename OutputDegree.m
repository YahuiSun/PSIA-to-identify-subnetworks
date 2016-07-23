% output degree of vertex

function [degree]=OutputDegree(set,N)

degree=zeros(N,1);  % degree of vertex 
for i=1:N
    for j=1:N
        if set(i,j)==1
            degree(i)=degree(i)+1;
        end
    end
end