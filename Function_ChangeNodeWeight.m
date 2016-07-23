% change the edge weight and node weight
function [node_weight]=Function_ChangeNodeWeight(L,set,node_weight,node_num,a2,b2,c2,d2,e2)

N=node_num;

% change node weight
for i=1:N
    node_weight(i)=fix(a2*b2^(c2*node_weight(i)+d2)+e2);
end
%







