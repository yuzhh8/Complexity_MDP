%% ###############################################  实际复杂度
%% 实际复杂度--多网格
function [RealTime,CPUtime,ComplexityrealT_Multi]=MT_Reality(gamma,epsilon,beta,c,L,delta,Q_0_infty)
epsilon_iter=(1-beta)*epsilon; % 迭代精度
Num_H=c^(L+delta); 

%% 离散化--多网格
h=ones(1,L+1); % 不同网格大小
Num_h=ones(1,L+1); % 网格数量
epsilon_itermulti=ones(1,L+1); % 不同误差
for l=1:L+1
    h(1,l)=c^(L+1-l)/Num_H;
    Num_h(1,l)=1/h(1,l);
    epsilon_itermulti(1,l)=c^(L+1-l)*epsilon_iter;
end

num_inter=100;
Complexityrealt_Multi=ones(1,L+1);
%Complexityreal.C_Multi=ones(1,L+1);
Q_0=Q_0_infty*ones(Num_h(1,1),Num_h(1,1));
tStart = cputime;
     tic
for l=1:L+1
      Q=ones(Num_h(1,l),Num_h(1,l));
      Q_0_min=ones(Num_h(1,l),1);  
    for t=1:num_inter
         for i=1:Num_h(1,l)
          Q_0_min(i,1)=min(Q_0(i,1:Num_h(1,l)+1-i));% 取每行最大值
         end
           Matrix_1=tril(1/Num_h(1,l)*ones(Num_h(1,l)-3,Num_h(1,l)-1));
        Matrix_2=diag(5/(6*Num_h(1,l))*ones(1,Num_h(1,l)-3));
        Matrix_3=diag(1/(6*Num_h(1,l))*ones(1,Num_h(1,l)-3));
        
        Matrix_4=Matrix_1+[zeros(1,Num_h(1,l)-3);Matrix_2;zeros(1,Num_h(1,l)-3)]'+[zeros(2,Num_h(1,l)-3);Matrix_3]';
        Matrix_5=[Num_h(1,l)-2:-1:2]/Num_h(1,l);
        Matrix_6=[Matrix_5;Matrix_4']';
        Matrix_P=[Matrix_6;[4/(3*Num_h(1,l)),1/Num_h(1,l)*ones(1,Num_h(1,l)-2),2/(3*Num_h(1,l))];zeros(2,Num_h(1,l))];
    for i=1:Num_h(1,l)
        SSi=1-1/(2*Num_h(1,l))-(i-1)/Num_h(1,l);
        SS_value_matrix = 1/Num_h(1,l)*[1:Num_h(1,l)];
        firstJ =[SSi,zeros(1,Num_h(1,l)-1)]+1/Num_h(1,l)*[0:Num_h(1,l)-1]'+repmat( SS_value_matrix,Num_h(1,l),1); 
            
     if i==1
        P=[1-1/(6*Num_h(1,l)),1/(6*Num_h(1,l)),zeros(1,Num_h(1,l)-2);1-1/Num_h(1,l),5/(6*Num_h(1,l)),1/(6*Num_h(1,l)),zeros(1,Num_h(1,l)-3);Matrix_P([1:Num_h(1,l)-2],:)];
      elseif i==2
     P=[1-1/Num_h(1,l),5/(6*Num_h(1,l)),1/(6*Num_h(1,l)),zeros(1,Num_h(1,l)-3);Matrix_P([1:Num_h(1,l)-1],:)];
      elseif i>=3&&i<Num_h(1,l)
        P=[[-(i-3)/Num_h(1,l),zeros(1,Num_h(1,l)-1)]+Matrix_P([1:Num_h(1,l)-i],:);Matrix_P([Num_h(1,l)-3],:);zeros(i-1,Num_h(1,l))];
      else
        P=[4/(3*Num_h(1,l)),1/Num_h(1,l)*ones(1,Num_h(1,l)-2),2/(3*Num_h(1,l));zeros(Num_h(1,l)-1,Num_h(1,l))];
    end
        
        Q_0_min_value =gamma *P* Q_0_min;
        
       Qj= firstJ.*P;
    
        Q(i,:) = sum(Qj,2)'+Q_0_min_value';
   end
       err = max(max(abs(Q - Q_0)));
       Q_0=Q;
       if err < epsilon_itermulti(1,l)
            Complexityrealt_Multi(1,l) = t; % 迭代复杂度
 %           Complexityreal.C_Multi(1,l)= Complexityreal.t_Multi(1,l)*Num_h(1,l)^3; % 查询复杂度
      [r,k] = size(Q);
       b= repmat(Q,c,1);
       b=reshape(b,r,c*k);
       d=repmat(b',c,1);
       d=reshape(d,c*k,c*r);
       Q_0=d'; % 初始值
            break
        end
    end
end
RealTime=toc;
CPUtime= cputime - tStart;
ComplexityrealT_Multi=Complexityrealt_Multi;
disp(['RealTime: ',num2str(RealTime)]);
   disp(['MCPUtime: ',num2str(CPUtime)]);
   fprintf('ComplexityrealT_Multi: %f',ComplexityrealT_Multi);
