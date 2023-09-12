function ComplexityInTheoryT_Multi=MT_Theory(gamma,epsilon,beta,c,L,delta,Q_0_infty)
K_r=1; % r Lipschitz常数
K_P=1; % P Lipschitz常数
r_infty=2; % r无穷范数
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

T=ones(1,L+1);
T(1,1)=ceil(log((1-gamma)*c^L*epsilon_iter/(r_infty+(1+gamma)*Q_0_infty))/log(gamma));
for i=2:L+1
    T(1,i)=ceil(log( (1-gamma)*epsilon_itermulti(1,i)/((h(1,i)+h(1,i-1))*(K_r+K_P*(gamma*epsilon_itermulti(1,i-1)+r_infty/(1-gamma)))+(1-gamma)*epsilon_itermulti(1,i-1)))/log(gamma));
end

ComplexityInTheoryT_Multi=T; % 总迭代复杂度