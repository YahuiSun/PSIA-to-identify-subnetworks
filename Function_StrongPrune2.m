function [prune_set]=Function_StrongPrune2(pcut_set,NW,N,L,R)

% Strong pruning
prune_set=pcut_set; % pruning set
pnw=NW; % pruning node weight
pcheck=zeros(N,1); % 1 means this node has been pruned
[pdegree]=Function_pdegree(N,prune_set,pcheck); % unprune degree
for i=1:N
    if pdegree(i)==0
        pcheck(i)=1; % At the first step, assume all the unconnected vertices have been pruned
    end
end

while sum(pcheck)<N-1 % there is always a vertex which cannot be pruned
    for i=1:N
        if pdegree(i)==1 && pcheck(i)==0 && i~=R % prune tree end i. if there is no root, the input R is 0
            for j=1:N
                if prune_set(i,j)==1 && pcheck(j)==0 % j is the predecessor of i
                    if L(i,j)>pnw(i)
                        prune_set(i,j)=0; prune_set(j,i)=0; % remove edge(i,j)
                        % remove subtree containing vertex i
                        LL=sparse(prune_set);
                        [S, C]=graphconncomp(LL);
                        for im=1:N
                            for in=im:N
                                if prune_set(im,in)==1 && C(im)==C(i) % remove subtree containing vertex i
                                    prune_set(im,in)=0; prune_set(in,im)=0;
                                end
                            end
                        end
                        %
                    else
                        pnw(j)=pnw(i)+pnw(j)-L(i,j);
                    end
                    pcheck(i)=1; % i has been pruned
                end
            end
        end
    end
    [pdegree]=Function_pdegree(N,prune_set,pcheck); % unprune degree
end