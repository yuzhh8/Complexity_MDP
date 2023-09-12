%% ############################################### 问题建模
%% 参数设置
% 模型参数设置（状态、行动空间，报酬函数，转移函数，折扣因子），设置后不再改变
S=[0,1]; % 状态空间
A=[0,1]; % 行动空间
%P(s,a,s1)=1; % 转移密度函数
gamma=0.7; % 折扣因子
K_r=1; % r Lipschitz常数
K_P=1; % P Lipschitz常数
r_infty=2; % r infty范数

% 可变参数（h,c,L,epsilon）
delta=5;
c=2;
L=4; % 网格层数
Num_H=c^(L+delta); % 离散网格数
epsilon=0.1; % 总误差
beta=0.5;
epsilon_dis=beta*epsilon; % 离散精度
epsilon_iter=(1-beta)*epsilon; % 迭代精度

Q_0_infty=0;

% Delta=(1-gamma)*(1-beta)*epsilon/(r_infty+(1+gamma)*Q_0_infty)*...
% ((1+c)/(1-beta)-1+(gamma*(1-gamma)*K_P*(1+c)*c^(L+1)*beta*epsilon)/((1-gamma)*K_r+K_P*r_infty))^4;
% disp(['Delta: ',num2str(Delta)]);

%% 离散化网格大小
h_max=epsilon_dis*(1-gamma)^2/((1-gamma)*K_r+K_P*r_infty); % 离散网格最大值
disp(['h_max: ',num2str(h_max)]);
if 1/Num_H<=h_max     % 判断所选网格大小是否达到要求
    disp(['离散网格大小: ',num2str(1/Num_H)]);
else
    disp('离散精度不满足要求，Error!!!!!!!!!!!!!');
return; 
end


%% 相关参数 Corollary2 确定的Q，L，\epsilon的零界点，以满足多网格复杂度小于单网格复杂度
% M<S  multigrid<single-grid
F=K_r+K_P*r_infty/(1-gamma); % 参数
x1=gamma*(1-gamma)*K_P*(1+c)*beta*epsilon/((1-gamma)*K_r+K_P*r_infty);
x2=(1+c)/(1-beta)-1;
x3=(1-gamma)*(1-beta)*epsilon/(r_infty+(1+gamma)*Q_0_infty);
x4=(1+c)*beta/(1-beta)+c+gamma*K_P*(1+c)*c^(L+1)*beta*epsilon/F;
Q_min=((1-gamma)*(1-beta)*epsilon*(x1*c^(L+1)+x2)^4-r_infty)/(1+gamma); % 最小初始值Q_min
y1=(1+gamma)*Q_0_infty+r_infty;
y2=y1/((1-gamma)*(1-beta)*epsilon);
L_max=log((y2^(1/4)-x2)/x1)/log(c)-1;   % 最大网格层数L_max
epsilon_max=y1/((1-gamma)*(1-beta)*(x2+x1*c^(L+1))^4); % 最大精度epsilon
disp(['Q_min: ',num2str(Q_min)]);
disp(['L_max: ',num2str(L_max)]);
disp(['epsilon_max: ',num2str(epsilon_max)]);

% S<M
%Q_max=((1-gamma)*(1-beta)*epsilon*(x1*c^(L+1)+x2)^4-r_infty)/(1+gamma);
%% 离散化--多网格
h=ones(1,L+1); % 不同网格大小
Num_h=ones(1,L+1); % 网格数量
epsilon_itermulti=ones(1,L+1); % 不同误差
for l=1:L+1
    h(1,l)=c^(L+1-l)/Num_H;
    Num_h(1,l)=1/h(1,l);
    epsilon_itermulti(1,l)=c^(L+1-l)*epsilon_iter;
end

% 状态行动空间离散化
% S_Multi=ones(L+1,Num_H+1); % 离散化区间端点
% for l=1:L+1
%     for i=1:Num_h(1,l)
%         S_Multi(l,i)=(i-1)/Num_h(1,l);
%     end
% end

%% ############################################### 离散化
%% 离散化--单网格   
S_Sigle=ones(1,Num_H+1); % 离散化区间端点 状态
for i=1:Num_H+1
    S_Sigle(1,i)=(i-1)/Num_H;
end

%% ############################################### 理论值-直接套用公式
%% 单网格
Condition_shouldLargerThan1=(1-gamma)*epsilon_iter/(r_infty+(1+gamma)*Q_0_infty);

% 迭代时间复杂度 t
ComplexityInTheory.t_Single = ceil(log( (1-gamma)*epsilon_iter/(r_infty+(1+gamma)*Q_0_infty))/log(gamma));
disp(['ComplexityInTheory.t_Single: ',num2str(ComplexityInTheory.t_Single)]);

% 查询复杂度
ComplexityInTheory.C_Single = ComplexityInTheory.t_Single*Num_H^3;
disp(['ComplexityInTheory.C_Single: ',num2str(ComplexityInTheory.C_Single)]);

%% 多网格
T=ones(1,L+1);
T(1,1)=ceil(log((1-gamma)*c^L*epsilon_iter/(r_infty+(1+gamma)*Q_0_infty))/log(gamma));
for i=2:L+1
    T(1,i)=ceil(log( (1-gamma)*epsilon_itermulti(1,i)/((h(1,i)+h(1,i-1))*(K_r+K_P*(gamma*epsilon_itermulti(1,i-1)+r_infty/(1-gamma)))+(1-gamma)*epsilon_itermulti(1,i-1)))/log(gamma));
end

ComplexityInTheory.t_Multi=sum(T); % 总迭代复杂度
disp(['ComplexityInTheory.t_Multi: ',num2str(ComplexityInTheory.t_Multi)]);

C=ones(1,L+1);
C(1,1)=T(1,1)*1/(h(1,1)^3);
for i=2:L+1
    C(1,i)=T(1,i)*1/(h(1,i)^3);
end

ComplexityInTheory.C_Multi=sum(C); % 总查询复杂度
disp(['ComplexityInTheory.C_Multi: ',num2str(ComplexityInTheory.C_Multi)]);

%% ###############################################  实际复杂度
%% 实际复杂度--单网格
num_inter1=10;
T=zeros(1,num_inter1);
Q_single=zeros(1,num_inter1);
Q_single1=zeros(1,num_inter1);
Q_0=Q_0_infty*ones(Num_H,Num_H);
Q_0_min=ones(Num_H,1);  
Q_infty=zeros(1,Num_H);
h=1/Num_H;
tStart = cputime;
  tic
for t=1:num_inter1
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
    if gamma* err/(1-gamma)>=epsilon_iter
       T(1,t)=t-1+log((1-gamma)*epsilon_iter/err)/log(gamma);
%         T(1,t)=t+1;
    else
        T(1,t)=t;
    end
    Q_single(1,t)=Q(Num_H/16,Num_H/16);
    Q_single1(1,t)=Q(Num_H/32,Num_H/32);
    Q_0=Q;
    Q_infty(1,t)=max(max(Q));
%     if err < epsilon_iter
%         ComplexityrealT_Single = t; % 迭代复杂度
%         break
%     end
end
toc
  StEnd = cputime - tStart;
   disp(['SCPUtime: ',num2str(StEnd)]);
%    disp(['Q: ',num2str(Q(Num_H/2,:))]);
   
min_Q=zeros(Num_H,1);  
index=zeros(Num_H,1);
for i=1:Num_H
[min_Q(i,1),index(i,1)]=min(Q(i,1:Num_H+1-i));
end
% [min_Q, index]=min(Q,[],2);
index=index'/Num_H;

figure(13); % 最优策略
plot(1/Num_H*(1:Num_H),index(1:Num_H),'r','LineWidth',1)
xlabel('State','FontSize',14)
ylabel('Optimal policy','FontSize',14)
legend_FontSize =legend('$\pi^*(s)$','interpreter','latex','FontSize',14);
set(legend_FontSize,'FontSize',14)
legend('boxoff')

% 多网格
num_inter=100; 
num=0;
Complexityrealt_Multi=ones(1,L+1);
%Complexityreal.C_Multi=ones(1,L+1);
Q_multi=zeros(1,num_inter);
Q_multi1=zeros(1,num_inter);
Q_0=Q_0_infty*ones(Num_h(1,1),Num_h(1,1));
Q_Minfty=zeros(1,Num_H);
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
       Q_multi(1,num+t)=Q(Num_H/16,Num_H/16);
        Q_multi1(1,num+t)=Q(Num_H/32,Num_H/32);
       Q_Minfty(1,num+t)=max(max(Q));
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
     num=num+Complexityrealt_Multi(1,l);
end


figure(15);
plot(1:num_inter1,Q_infty(1:num_inter1),'-*r','LineWidth',1)
hold on
plot(1:sum(Complexityrealt_Multi),Q_Minfty(1:sum(Complexityrealt_Multi)),'--og','LineWidth',1)
xlabel('Numbers of iteration ','FontSize',14)
ylabel('$\|Q\|_\infty$','interpreter','latex','FontSize',14)
grid on
legend_FontSize =legend('Single-grid','Multigrid');
set(legend_FontSize,'FontSize',14)
legend('boxoff') 

figure(16);
plot(1:num_inter1,Q_single(1:num_inter1),'-*r','LineWidth',1)
hold on
plot(1:sum(Complexityrealt_Multi),Q_multi(1:sum(Complexityrealt_Multi)),'--*g','LineWidth',1)
hold on
plot(1:num_inter1,Q_single1(1:num_inter1),'-or','LineWidth',1)
hold on
plot(1:sum(Complexityrealt_Multi),Q_multi1(1:sum(Complexityrealt_Multi)),'--og','LineWidth',1)
xlabel('Numbers of iteration','FontSize',14)
ylabel('$Q$','interpreter','latex','FontSize',14)
grid on
legend_FontSize =legend('Single(1/16,1/16)','Multi(1/16,1/16)','Single(1/32,1/32)','Multi(1/32,1/32)');
set(legend_FontSize,'FontSize',14)
legend('boxoff') 

ComplexityrealT_Multi=Complexityrealt_Multi;
% disp(['RealTime: ',num2str(RealTime)]);
%    disp(['MCPUtime: ',num2str(CPUtime)]);
   fprintf('ComplexityrealT_Multi: %f',ComplexityrealT_Multi);
   disp(['SQ*: ',num2str(max(min_Q))]);