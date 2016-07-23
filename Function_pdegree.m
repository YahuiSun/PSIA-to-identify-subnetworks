function  [degree]=Function_pdegree(N,set,pcheck)

degree=zeros(N,1); % unprune degree; connected vertices which have not been pruned

for i=1:N
    for j=1:N
        if set(i,j)==1 & pcheck(j)==0
            degree(i)=degree(i)+1;
        end
    end
end







