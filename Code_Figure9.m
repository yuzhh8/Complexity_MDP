delta=5;
c=2;
L=4; % 网格层数
Num_H=c^(L+delta); % 离散网格数
Num_h=ones(1,L+1); % 网格数量
Num_h3=ones(1,L+1);
for l=1:L+1
    Num_h(1,l)=Num_H/c^(L+1-l);
    Num_h3(1,l)=(Num_H/c^(L+1-l))^3;
end

%% single-grid
ST_T=[];
TR2=[];
TR3=[];
ST_R=[];
%% 初始值的影响
    for Q_0_infty=0:0.2:4
T=Real_time(0.7,0.1,0.5,2,4,5,Q_0_infty);
[SRealTime,SCPUTime,ComplexityrealT_Single]=ST_Reality(0.7,0.1,0.5,2,4,5,Q_0_infty);
ST_R(end+1)=ComplexityrealT_Single;
ST_T(end+1)=ST_Theory(0.7,0.1,0.5,Q_0_infty);
TR2(end+1)=T(1,2);
TR3(end+1)=T(1,5);
    end

figure(13);
plot(0:0.2:4,ST_R(:),'-*r','LineWidth',1)
hold on
plot(0:0.2:4,ST_T(:),'--*g','LineWidth',1)
hold on
plot(0:0.2:4,TR2(:),'--og','LineWidth',1)
hold on
plot(0:0.2:4,TR3(:),'--dg','LineWidth',1)
xlabel('Initial value $Q_0(s,a)=q_0$','interpreter','latex','FontSize',14)
ylabel('Iteration complexity','FontSize',14)
grid on
legend_FontSize =legend('Actual $\mathcal{T}^{(sg)}$','Theoretical $\hat{\mathcal{T}}^{(sg)}(0)$','Theoretical $\hat{\mathcal{T}}^{(sg)}(2)$','Theoretical $\hat{\mathcal{T}}^{(sg)}(5)$','interpreter','latex');
set(legend_FontSize,'FontSize',14)
legend('boxoff')

%% multigrid
TR1=[];
TR2=[];
TR3=[];
MT_T=[];
MT_R=[];
%% 初始值的影响
    for Q_0_infty=0:0.2:4
MT=Multi_Real_time(0.1,0.5,2,4,5,Q_0_infty);
[SRealTime,SCPUTime,ComplexityrealT_Multi]=MT_Reality(0.7,0.1,0.5,2,4,5,Q_0_infty);
ComplexityInTheoryT_Multi=MT_Theory(0.7,0.1,0.5,2,4,5,Q_0_infty);
MT_T(end+1)=sum(ComplexityInTheoryT_Multi);
MT_R(end+1)=sum(ComplexityrealT_Multi);
TR1(end+1)=MT(1,1);
TR2(end+1)=MT(1,2);
TR3(end+1)=MT(4,2);
    end

figure(14);
plot(0:0.2:4,MT_R(:),'-*r','LineWidth',1)
hold on
plot(0:0.2:4,MT_T(:),'--*g','LineWidth',1)
hold on
plot(0:0.2:4,TR1(:),'--.g','LineWidth',1)
hold on
plot(0:0.2:4,TR2(:),'--og','LineWidth',1)
hold on
plot(0:0.2:4,TR3(:),'--dg','LineWidth',1)
xlabel('Initial value $Q_0(s,a)=q_0$','interpreter','latex','FontSize',14)
ylabel('Iteration complexity','FontSize',14)
grid on
legend_FontSize =legend('Actual $\mathcal{T}^{(mtg)}$','Theoretical $\hat{\mathcal{T}}^{(mtg)}(0,0)$','Theoretical $\hat{\mathcal{T}}^{(mtg)}(0,1)$','Theoretical $\hat{\mathcal{T}}^{(mtg)}(0,2)$','Theoretical $\hat{\mathcal{T}}^{(mtg)}(3,2)$','interpreter','latex');
set(legend_FontSize,'FontSize',14)
legend('boxoff')