function [FinalPrize,TotalCost]=Function_solution(FinalPrize,TotalCost,Auto,multiple,final_set,L,N,R,node_weight,Terminal,lefttime)


% solution quality 1: maxamize prize minus edge cost
FinalPrize(Auto,multiple)=0;
for i=1:N
    for j=i:N
        FinalPrize(Auto,multiple)=FinalPrize(Auto,multiple)-final_set(i,j)*L(i,j);
    end
end
degree=zeros(N,1);
for i=1:N
    for j=1:N
        if final_set(i,j)==1
            degree(i)=degree(i)+1;
        end
    end
end
for i=1:N
    if degree(i)>0
        FinalPrize(Auto,multiple)=FinalPrize(Auto,multiple)+node_weight(i);
    else
        if Terminal(i)==1
            fprintf(['Terminal not included!\n'])
        end
    end
end
if R~=0
    if degree(R)==0
        FinalPrize(Auto,multiple)=node_weight(R); % the solution only have the root vertex
    end
else
    a=max(node_weight);
    if FinalPrize(Auto,multiple)<=a
        FinalPrize(Auto,multiple)=a; % the solution only have one vertex
    end
end
% solution quality 2: minimize edge cost add not-included prize
sumprize=sum(node_weight);
TotalCost(Auto,multiple)=sumprize-FinalPrize(Auto,multiple);

fprintf(['FinalPrize(',num2str(multiple),')=',num2str(FinalPrize(Auto,multiple)),' Max Best=',num2str(max(FinalPrize(Auto,1:multiple))),'  TotalCost(',num2str(multiple),')=',num2str(TotalCost(Auto,multiple)),' Min Best=',num2str(min(TotalCost(Auto,1:multiple))),'  Left Time=',num2str(lefttime),' hours.\n'])
%