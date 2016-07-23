clear all;  clc;

% the PSIA algorithm to identify subnetworks
load D_01_b; % load the drug similarity network
dataname='D_01_b'; % the label of the saved result/data


% initialization
fe_PO=0; % function evaluation time
fit_PO=zeros(1e7,1); % the fitness evaluation value (the net-weight of the identified subnetwork) after each function evaluation time
N=node_num; % the total number of vertices in the initial graph
Q=ones(N); % the flux matrix
Maxw=max(node_weight); % the biggest node weight
terminal_num=sum(Terminal); % the number of terminals

%  the orders of terminals and vertices
foodpoint=Terminal;
foodnum=sum(foodpoint);
food=zeros(foodnum,1);
r=0;
for i=1:N
    if foodpoint(i)==1
        r=r+1;
        food(r)=i;  %  the No.r terminal is the No.i vertex
    end
end
%

% parameters in PSIA
I=1e0;
kk=125;  % iteration times to identify each subnetwork
cutoff=1e-3;   % cutoff value of conductivity
alpha=0.2; % parameter in the evoluation function of conductivity
sigma=1; % parameter in the evoluation function of conductivity
%


Rset=set; % the connectivity matrix of the initial graph
RL=L; % the edge length matrix of the initial graph


whole=0;
runID=001; % the ID of the running times or the identified subnetwork
while whole<5 % run 5 times to identify 5 subnetworks
    whole=whole+1;   % the running times (an identified subnetwork can be identified in each running time)
    runID=whole; % the ID of the running times or the identified subnetwork
    set=Rset; % the connectivity matrix of the initial graph
    L=RL; % the edge length matrix of the initial graph
    
    % give each edge conductivity a random value
    D=ones(N);
    randomD=10e3;
    for i=1:N
        for j=i:N
            D(i,j)=randomD*rand(1);
            D(j,i)=D(i,j);
        end
    end
    %
    
    exist=ones(N,1); % 1 means a vertex is in the subnetwork, 0 means not
    
    % iteration of D and P
    for k=1:kk  % the iteration times
        NUM=zeros(N,1);
        [degree]=Function_OutputDegree(set,N);  % degree of vertex
        for i=1:N
            if degree(i)==0
                exist(i)=0; % 1 means a vertex is in the subnetwork, 0 means not
            end
        end
        vertex=sum(exist);  % number of vertices in the subnetwork
        
        % unequal possibility 3: choose the sink and source nodes probabilistically
        adjacentlength=zeros(foodnum,1); % the total cost of edges linked to a terminal
        for i=1:foodnum
            for j=1:N
                if set(food(i),j)==1
                    adjacentlength(i)=adjacentlength(i)+L(food(i),j);
                end
            end
        end
        [B,ind]=sort(adjacentlength);   % re-order adjacentlength in the ascending order,  B is the ordered vector, ind is the intex
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
                if random>sumad(i-1) && random<=sumad(i)
                    sink=food(ind(foodnum+1-i)); luck=ind(foodnum+1-i);
                end
            end
        end
        %
        
        % calculate the pressures using the network Poisson equation
        A=zeros(vertex-1);
        g=0;
        for w1=1:N
            if exist(w1)==1   % only calculate the existing vertice
                if w1==sink % neglect the vertex if it's the sink point
                else
                    g=g+1; % row number of A
                    h=0;
                    NUM(w1)=g; %  record the row number's relationship with vertex number, which will be used when using pressure
                    for w3=1:N  % the vertex corresponds to the columes of A
                        if exist(w3)==1  % only calculate the existing point
                            if w3==sink(1)  % neglect the vertex if it's the sink point
                            else
                                h=h+1; % colume number of A
                                if w3==w1   %  the row vertex is the colume vertex
                                    for i=1:N
                                        if set(w1,i)==1
                                            A(g,h)=A(g,h)-D(w1,i)/(L(w1,i)-node_weight(w1)/degree(w1)-node_weight(i)/degree(i)+2*Maxw);
                                        end
                                    end
                                else    % the row vertex is not the colume vertex
                                    if set(w1,w3)==1
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
        for i=1:foodnum
            if i==luck
            else
                X(NUM(food(i)))=-I;
            end
        end
        
        P=A\X; % the pressures
        
        
        % calculate the flux, update the conductivity and cut edges
        for i=1:N
            for j=1:N
                if set(i,j)==1
                    if i==sink
                        Q(i,j)=D(i,j)*(0-P(NUM(j)))/(L(i,j)-node_weight(i)/degree(i)-node_weight(j)/degree(j)+2*Maxw);
                        D(i,j)=(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff
                            set(i,j)=0; set(j,i)=0;
                        end
                    end
                    if j==sink
                        Q(i,j)=D(i,j)*(P(NUM(i))-0)/(L(i,j)-node_weight(i)/degree(i)-node_weight(j)/degree(j)+2*Maxw);
                        D(i,j)=(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff
                            set(i,j)=0; set(j,i)=0;
                        end
                    end
                    if i~=sink && j~=sink
                        Q(i,j)=D(i,j)*(P(NUM(i))-P(NUM(j)))/(L(i,j)-node_weight(i)/degree(i)-node_weight(j)/degree(j)+2*Maxw);
                        D(i,j)=(D(i,j)+alpha*abs(Q(i,j))-sigma*D(i,j));
                        if D(i,j)<cutoff
                            set(i,j)=0; set(j,i)=0;
                        end
                    end
                end
            end
        end
        
        % move out disconnected part
        [set]=Function_prim(set,foodpoint,N);
        %
        
        fe_PO=fe_PO+1; % function evaluation times
        
        TotalPrice=0; % the net-weight of the subnetwork
        for i=1:N
            for j=i+1:N
                TotalPrice=TotalPrice-set(i,j)*L(i,j);
            end
        end
        [degree]=Function_OutputDegree(set,N);
        for i=1:N
            if degree(i)>0
                TotalPrice=TotalPrice+node_weight(i);
            end
        end
        fit_PO(fe_PO)=TotalPrice;
        
    end  % inner iteration
    
    % find the SMT
    [set]=Function_SMTpost(N,set,RL,Rset);
    %
    
    Analysis_set=set;
    fprintf(['kk=',num2str(k),'\n'])
    [degree]=Function_OutputDegree(Analysis_set,N);
    IncludedVertex1=zeros(N,1);
    for i=1:N
        if degree(i)>0
            IncludedVertex1(i)=1;
        end
    end
    fprintf(['There are ',num2str(sum(IncludedVertex1)),' out of ',num2str(N),' vertices included.\n']);
    removeV=N-sum(IncludedVertex1);
    fprintf(['There are ',num2str(removeV),' out of ',num2str(N),' vertices removed.\n']);
    fprintf(['\n'])
    save(['Result_',dataname,'_NW_',num2str(node_weight(1)),...
        '_runID_',num2str(runID),'_FitnessValue_',num2str(TotalPrice),'_k_',num2str(k),...
        '_RemovedVertices_',num2str(removeV)]);
    
end   %  whole iteration

