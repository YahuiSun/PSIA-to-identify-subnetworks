function [new_index,new_surp,new_cpn_num,new_active]=Function_merge_cpn(T1,T2,N,index,surp,cpn_num,active)

new_cpn_num=cpn_num-1;
a=1;b=0;
for i=1:cpn_num
    if i~=T1 & i~=T2
       to(i)=a; % cpn i to new_cpn a
       a=a+1; 
    else
        if b==0
            b=1;
            mergeI=a;
            to(i)=mergeI; % T1 to new_cpn mergeI
            a=a+1;
        else
            to(i)=mergeI; % T2 to new_cpn mergeI
        end
    end
end

for i=1:N
    new_index(i)=to(index(i));
end

new_surp=zeros(new_cpn_num,1);
new_active=zeros(new_cpn_num,1);
for i=1:cpn_num
     new_surp(to(i))=new_surp(to(i))+surp(i);
     if i~=T1 & i~=T2
        new_active(to(i))=active(i);
     else
        new_active(to(i))=min(active(T1)+active(T2),1);
     end
end




