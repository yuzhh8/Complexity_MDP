%% ###############################################  实际复杂度
%% 实际复杂度--单网格
function [RealTime,CPUtime,ComplexityrealT_Single]=ST_Reality(gamma,epsilon,beta,c,L,delta,Q_0_infty)
epsilon_iter=(1-beta)*epsilon; % 迭代精度
Num_H=c^(L+delta); 

S_Sigle=ones(1,Num_H+1); % 离散化区间端点
for i=1:Num_H+1
    S_Sigle(1,i)=(i-1)/Num_H;
end

num_inter=100;
Q_0=Q_0_infty*ones(Num_H,Num_H);
Q_0_min=ones(Num_H,1);  
h=1/Num_H;
tStart = cputime;
  tic
for t=1:num_inter
    Q=zeros(Num_H,Num_H);
    for i=1:Num_H
    Q_0_min(i,1)=min(Q_0(i,1:Num_H+1-i));
    end
     Matrix_1=tril(1/Num_H*ones(Num_H-3,Num_H-1));
        Matrix_2=diag(5/(6*Num_H)*ones(1,Num_H-3));
        Matrix_3=diag(1/(6*Num_H)*ones(1,Num_H-3));
        
        Matrix_4=Matrix_1+[zeros(1,Num_H-3);Matrix_2;zeros(1,Num_H-3)]'+[zeros(2,Num_H-3);Matrix_3]';
        Matrix_5=[Num_H-2:-1:2]/Num_H;
        Matrix_6=[Matrix_5;Matrix_4']';
           Matrix_P=[Matrix_6;[4/(3*Num_H),1/Num_H*ones(1,Num_H-2),2/(3*Num_H)];zeros(2,Num_H)];
   for i=1:Num_H
        SSi=1-1/(2*Num_H)-S_Sigle(1,i);
        SS_value_matrix = S_Sigle(2:Num_H+1);
        firstJ =[SSi,zeros(1,Num_H-1)]+S_Sigle(1:Num_H)'+repmat( SS_value_matrix,Num_H,1); 
            
   if i==1
        P=[1-1/(6*Num_H),1/(6*Num_H),zeros(1,Num_H-2);1-1/Num_H,5/(6*Num_H),1/(6*Num_H),zeros(1,Num_H-3);Matrix_P([1:Num_H-2],:)];
   elseif i==2
     P=[1-1/Num_H,5/(6*Num_H),1/(6*Num_H),zeros(1,Num_H-3);Matrix_P([1:Num_H-1],:)];
      elseif i>=3&&i<Num_H
        P=[[-(i-3)/Num_H,zeros(1,Num_H-1)]+Matrix_P([1:Num_H-i],:);Matrix_P([Num_H-3],:);zeros(i-1,Num_H)];
    else
   P=[4/(3*Num_H),1/Num_H*ones(1,Num_H-2),2/(3*Num_H);zeros(Num_H-1,Num_H)];
   end
        
        Q_0_min_value =gamma *P* Q_0_min;
        
       Qj= firstJ.*P;
    
        Q(i,:) = sum(Qj,2)'+Q_0_min_value';
    end
    err = max(max(abs(Q - Q_0)));
    Q_0=Q;
    if err < epsilon_iter
        ComplexityrealT_Single = t; % 迭代复杂度
        break
    end
end
RealTime=toc;
  CPUtime = cputime - tStart;
    disp(['RealTime: ',num2str(RealTime)]);
   disp(['SCPUtime: ',num2str(CPUtime)]);
    disp([' ComplexityrealT_Single: ',num2str( ComplexityrealT_Single)]);