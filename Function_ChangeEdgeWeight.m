% change the edge weight and node weight
function [L]=Function_ChangeEdgeWeight(L,set,node_weight,node_num,a2,b2,c2,d2,e2)

N=node_num;

% change edge weight
for i=1:N
    for j=i:N
        if set(i,j)==1
            L(i,j)=fix(a2*b2^(c2*L(i,j)+d2)+e2);
            L(j,i)=L(i,j);
        end
    end
end
%

