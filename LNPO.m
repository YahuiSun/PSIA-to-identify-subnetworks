clear all;
clc;

% LNPO
load C01_4_548_1391_22_0;
dataname='C01_4_548_1391_22_0';
Target=zeros(node_num,1);
node_weight=80*ones(node_num,1);




% change the edge weight and node weight
% a1=15; b1=10; c1=1; d1=-49; e1=0;
%
% [node_weight]=Function_ChangeNodeWeight(L,set,node_weight,node_num,a1,b1,c1,d1,e1);

% [node_weight]=Function_ChangeNodeWeight(L,set,node_weight,node_num,1,0,1,1,0);
%




% optimal=1e6;  % optimal solution
%
Timelimit=3600*672-3600;
FElimit=1e6;
%

% initialization
fe_PO=0; % function evaluation
fit_PO=0;
N=node_num;
edgevalue=ones(N,N);
success=0;
Maxw=max(node_weight);
terminal_num=sum(Terminal);
foodpoint=Terminal;
foodnum=sum(foodpoint);
edgevalue=ones(N,N);
food=zeros(foodnum,1);
r=0;
for i=1:N
    if foodpoint(i)==1
        r=r+1;
        food(r)=i;  %%%%%%  the num r food is the num i vertice
    end
end
%




for i=1:N
    if Target(i)>0
%         fprintf(['Vertex ',num2str(i),' is a target. It has a weight of ',num2str(node_weight(i)),'.\n']);
    end
end
fprintf(['There are ',num2str(sum(Target)),' targets.\n']);
fprintf(['\n'])



runID=001;


% algorithm parameter
I=1e0;
kk=125;  % inner iteration time
cutoff=1e-3;   % cutoff value of D
% evolution of conductivity   some value of alpha and sigma may make the
% whole program stop outputing results in the middle, don't know why
alpha=0.2;
sigma=1;
%

tic;
Rset=set; RL=L; fuseL=0;

%
whole=0;
while whole<5
    whole=whole+1;   % outer iteration time
    
    runID=whole;
    
    
    set=Rset; L=RL; D=set;
    
    randomD=10e3;
    for i=1:N
        for j=i:N
            D(i,j)=randomD*rand(1);
            D(j,i)=D(i,j);
        end
    end
    
    % whether a vertex exist
    exist=ones(N,1); % 1 means exist, 0 means not
    
    % iteration of D and P
    for k=1:kk  % inner iteration time
        NUM=zeros(N,1);
        % whether a vertex exist
        degree=zeros(N,1);  % degree of vertex
        [degree]=OutputDegree(set,N);
        for i=1:N
            if degree(i)==0
                exist(i)=0;
            end
        end
        vertex=sum(exist);  % number of existing vertex
        
        % unequal possibility 3: choose the sink node and the source nodes
        adjacentlength=zeros(foodnum,1);
        for i=1:foodnum
            for j=1:N
                if set(food(i),j)==1
                    adjacentlength(i)=adjacentlength(i)+L(food(i),j);
                end
            end
        end
        [B,ind]=sort(adjacentlength);   %%%  ascending order,  B is the ordered vector, ind is the intex
        TotalAdjacent=0;
        for i=1:foodnum
            TotalAdjacent=TotalAdjacent+B(i);
        end
        sumad=zeros(foodnum,1);
        for i=1:foodnum
            for j=1:i
                sumad(i)=sumad(i)+B(j);
            end
        end
        random=ceil(rand(1)*TotalAdjacent);
        if random<=sumad(1)
            sink=food(ind(foodnum)); luck=ind(foodnum);
        else
            for i=2:foodnum
                if random>sumad(i-1) & random<=sumad(i)
                    sink=food(ind(foodnum+1-i)); luck=ind(foodnum+1-i);
                end
            end
        end
        
        %%%%%%%%%%%%%%   calculate the pressure
        A=zeros(vertex-1);
        g=0;
        for w1=1:N
            if exist(w1)==1   %%% only calculate the existing vertice
                if w1==sink %%% neglect the vertex if it's the sink point
                    ;
                else
                    g=g+1; %% row number of A
                    h=0;
                    NUM(w1)=g; %%%  record the row number's relationship with vertex number, which will be used when using pressure
                    for w3=1:N  %%%%%%%  the vertex corresponds to the columes of A
                        if exist(w3)==1  %%% only calculate the existing point
                            if w3==sink(1)  %%% neglect the vertex if it's the sink point
                                ;
                            else
                                h=h+1; %% colume number of A
                                if w3==w1   %%%  the row vertex is the colume vertex
                                    for i=1:N
                                        if set(w1,i)==1
                                            %                                             A(g,h)=A(g,h)-D(w1,i)/L(w1,i);
                                            A(g,h)=A(g,h)-D(w1,i)/(L(w1,i)-node_weight(w1)/degree(w1)-node_weight(i)/degree(i)+2*Maxw);
                                        end
                                    end
                                else    %%%%   the row vertex is not the colume vertex
                                    if set(w1,w3)==1
                                        %                                         A(g,h)=D(w1,w3)/L(w1,w3);
                                        A(g,h)=D(w1,w3)/(L(w1,w3)-node_weight(w1)/degree(w1)-node_weight(w3)/degree(w3)+2*Maxw);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        X=zeros(vertex-1,1);
        %%%%%%%%%%%%% all foodpoints that are not sink are sources
        for i=1:foodnum
            if i==luck
                ;
            else
                X(NUM(food(i)))=-I;
            end
        end
        
        P=A\X;
        
        
        %%%%calculate the flux and update the conductivity
        for i=1:N
            for j=1:N
                if set(i,j)==1
                    if i==sink
                        Q(i,j)=D(i,j)*(0-P(NUM(j)))/(L(i,j)-node_weight(i)/degree(i)-node_weight(j)/degree(j)+2*Maxw); % equation 4
                        D(i,j)=edgevalue(i,j)*(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff
                            set(i,j)=0; set(j,i)=0;
                        end
                    end
                    if j==sink
                        Q(i,j)=D(i,j)*(P(NUM(i))-0)/(L(i,j)-node_weight(i)/degree(i)-node_weight(j)/degree(j)+2*Maxw); % equation 4
                        D(i,j)=edgevalue(i,j)*(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff
                            set(i,j)=0; set(j,i)=0;
                        end
                    end
                    if i~=sink && j~=sink
                        Q(i,j)=D(i,j)*(P(NUM(i))-P(NUM(j)))/(L(i,j)-node_weight(i)/degree(i)-node_weight(j)/degree(j)+2*Maxw); % equation 4
                        D(i,j)=edgevalue(i,j)*(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff
                            set(i,j)=0; set(j,i)=0;
                        end
                    end
                end
            end
        end
        
        % move out disconnected part, or there may be inaccurate calculation warnnings
        [set]=function_prim(set,foodpoint,N);
        %
        
        fe_PO=fe_PO+1; % function evaluation
        
        TotalPrice=0; % solution
        for i=1:N
            for j=i+1:N
                TotalPrice=TotalPrice-set(i,j)*L(i,j);
            end
        end
        [degree]=OutputDegree(set,N);
        for i=1:N
            if degree(i)>0
                TotalPrice=TotalPrice+node_weight(i);
            end
        end
        fit_PO(fe_PO)=TotalPrice;
        
        % obtain the pruned FE and last time
        [bestsofar indexM]=max(fit_PO(1:fe_PO));
        indexM=mod(indexM,kk);
        FE=indexM*whole;
        if FE>FElimit
            success=3;
            break;
        end
        Time=toc;  %%% record simulation time
        LastTime=(FElimit-FE)/FE*Time/3600; % hours
        %
        
        if TotalPrice==bestsofar
            MaxSet=set;
            for i=1:N
                for j=1:N
                    if set(i,j)==1
                        edgevalue(i,j)=1.2;
                    else
                        edgevalue(i,j)=0.8;
                    end
                end
            end
        end
        
        if k==1
            LNPO_set=set;
            % analysis
            Analysis_set=LNPO_set;
            % GW_set : the solution of GW
            % GW_MSTp_set : the solution of GW + MSTposeProcessing
            % LNPO_set : the solution of LNPO
            % LPPO_set : the solution of LPPO
            % MST_set : the solution of MST + pruning
            fprintf(['kk=',num2str(k),'\n'])
            [degree]=OutputDegree(Analysis_set,N);
            IncludedVertex1=zeros(N,1);
            for i=1:N
                if degree(i)>0
                    IncludedVertex1(i)=1;
                    %         fprintf(['Vertex ',num2str(i),' is included.\n']);
                end
            end
            fprintf(['There are ',num2str(sum(IncludedVertex1)),' out of ',num2str(N),' vertices included.\n']);
            removeV=N-sum(IncludedVertex1);
            fprintf(['There are ',num2str(removeV),' out of ',num2str(N),' vertices removed.\n']);
            a=0;
            for i=1:N
                if IncludedVertex1(i)==1 & Target(i)==1
                    a=a+1;
%                     fprintf(['Vertex ',num2str(i),' has been detected.\n']);
                end
            end
            fprintf(['There are ',num2str(a),' out of ',num2str(sum(Target)),' targets detected.\n']);
            fprintf(['\n'])
%             fprintf(['Initial Target rate: ',num2str(sum(Target)/N),'. Identified target rate: ',num2str(a/sum(IncludedVertex1)),'.\n']);
            improve_rate1=fix((a/sum(IncludedVertex1)-sum(Target)/N)*10000); % /1000
%             fprintf(['Improved target rate: ',num2str(improve_rate1),' /10000.\n']);
            %
            save(['Result_',dataname,'_NW_',num2str(node_weight(1)),...
                '_runID_',num2str(runID),'_FitnessValue_',num2str(TotalPrice),'_k_',num2str(k),...
                '_RemovedVertices_',num2str(removeV),...
                '_DetectedTarget_',num2str(a),...
                '_ImprovedTargetRate_',num2str(improve_rate1)]);
        elseif k==5
            LNPO_set=set;
            % analysis
            Analysis_set=LNPO_set;
            % GW_set : the solution of GW
            % GW_MSTp_set : the solution of GW + MSTposeProcessing
            % LNPO_set : the solution of LNPO
            % LPPO_set : the solution of LPPO
            % MST_set : the solution of MST + pruning
            fprintf(['kk=',num2str(k),'\n'])
            [degree]=OutputDegree(Analysis_set,N);
            IncludedVertex1=zeros(N,1);
            for i=1:N
                if degree(i)>0
                    IncludedVertex1(i)=1;
                    %         fprintf(['Vertex ',num2str(i),' is included.\n']);
                end
            end
            fprintf(['There are ',num2str(sum(IncludedVertex1)),' out of ',num2str(N),' vertices included.\n']);
            removeV=N-sum(IncludedVertex1);
            fprintf(['There are ',num2str(removeV),' out of ',num2str(N),' vertices removed.\n']);
            a=0;
            for i=1:N
                if IncludedVertex1(i)==1 & Target(i)==1
                    a=a+1;
%                     fprintf(['Vertex ',num2str(i),' has been detected.\n']);
                end
            end
            fprintf(['There are ',num2str(a),' out of ',num2str(sum(Target)),' targets detected.\n']);
            fprintf(['\n'])
%             fprintf(['Initial Target rate: ',num2str(sum(Target)/N),'. Identified target rate: ',num2str(a/sum(IncludedVertex1)),'.\n']);
            improve_rate1=fix((a/sum(IncludedVertex1)-sum(Target)/N)*10000); % /1000
%             fprintf(['Improved target rate: ',num2str(improve_rate1),' /10000.\n']);
            %
            save(['Result_',dataname,'_NW_',num2str(node_weight(1)),...
                '_runID_',num2str(runID),'_FitnessValue_',num2str(TotalPrice),'_k_',num2str(k),...
                '_RemovedVertices_',num2str(removeV),...
                '_DetectedTarget_',num2str(a),...
                '_ImprovedTargetRate_',num2str(improve_rate1)]);
        elseif k==25
            LNPO_set=set;
            % analysis
            Analysis_set=LNPO_set;
            % GW_set : the solution of GW
            % GW_MSTp_set : the solution of GW + MSTposeProcessing
            % LNPO_set : the solution of LNPO
            % LPPO_set : the solution of LPPO
            % MST_set : the solution of MST + pruning
            fprintf(['kk=',num2str(k),'\n'])
            [degree]=OutputDegree(Analysis_set,N);
            IncludedVertex1=zeros(N,1);
            for i=1:N
                if degree(i)>0
                    IncludedVertex1(i)=1;
                    %         fprintf(['Vertex ',num2str(i),' is included.\n']);
                end
            end
            fprintf(['There are ',num2str(sum(IncludedVertex1)),' out of ',num2str(N),' vertices included.\n']);
            removeV=N-sum(IncludedVertex1);
            fprintf(['There are ',num2str(removeV),' out of ',num2str(N),' vertices removed.\n']);
            a=0;
            for i=1:N
                if IncludedVertex1(i)==1 & Target(i)==1
                    a=a+1;
%                     fprintf(['Vertex ',num2str(i),' has been detected.\n']);
                end
            end
            fprintf(['There are ',num2str(a),' out of ',num2str(sum(Target)),' targets detected.\n']);
            fprintf(['\n'])
%             fprintf(['Initial Target rate: ',num2str(sum(Target)/N),'. Identified target rate: ',num2str(a/sum(IncludedVertex1)),'.\n']);
            improve_rate1=fix((a/sum(IncludedVertex1)-sum(Target)/N)*10000); % /1000
%             fprintf(['Improved target rate: ',num2str(improve_rate1),' /10000.\n']);
            %
            save(['Result_',dataname,'_NW_',num2str(node_weight(1)),...
                '_runID_',num2str(runID),'_FitnessValue_',num2str(TotalPrice),'_k_',num2str(k),...
                '_RemovedVertices_',num2str(removeV),...
                '_DetectedTarget_',num2str(a),...
                '_ImprovedTargetRate_',num2str(improve_rate1)]);
        elseif k==125
            LNPO_set=set;
            % analysis
            Analysis_set=LNPO_set;
            % GW_set : the solution of GW
            % GW_MSTp_set : the solution of GW + MSTposeProcessing
            % LNPO_set : the solution of LNPO
            % LPPO_set : the solution of LPPO
            % MST_set : the solution of MST + pruning
            fprintf(['kk=',num2str(k),'\n'])
            [degree]=OutputDegree(Analysis_set,N);
            IncludedVertex1=zeros(N,1);
            for i=1:N
                if degree(i)>0
                    IncludedVertex1(i)=1;
                    %         fprintf(['Vertex ',num2str(i),' is included.\n']);
                end
            end
            fprintf(['There are ',num2str(sum(IncludedVertex1)),' out of ',num2str(N),' vertices included.\n']);
            removeV=N-sum(IncludedVertex1);
            fprintf(['There are ',num2str(removeV),' out of ',num2str(N),' vertices removed.\n']);
            a=0;
            for i=1:N
                if IncludedVertex1(i)==1 & Target(i)==1
                    a=a+1;
%                     fprintf(['Vertex ',num2str(i),' has been detected.\n']);
                end
            end
            fprintf(['There are ',num2str(a),' out of ',num2str(sum(Target)),' targets detected.\n']);
            fprintf(['\n'])
%             fprintf(['Initial Target rate: ',num2str(sum(Target)/N),'. Identified target rate: ',num2str(a/sum(IncludedVertex1)),'.\n']);
            improve_rate1=fix((a/sum(IncludedVertex1)-sum(Target)/N)*10000); % /1000
%             fprintf(['Improved target rate: ',num2str(improve_rate1),' /10000.\n']);
            %
            save(['Result_',dataname,'_NW_',num2str(node_weight(1)),...
                '_runID_',num2str(runID),'_FitnessValue_',num2str(TotalPrice),'_k_',num2str(k),...
                '_RemovedVertices_',num2str(removeV),...
                '_DetectedTarget_',num2str(a),...
                '_ImprovedTargetRate_',num2str(improve_rate1)]);
        end
        
        
        
    end  % inner iteration
    
end   %%%%  whole iteration


% Time=toc;  %%% record simulation time
% save(['Result_',dataname,'_LNPO_WholeData_whole_',num2str(whole),'_BestSolution_',num2str(bestsofar)]);  %%%%  save data

