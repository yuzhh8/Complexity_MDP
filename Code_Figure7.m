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

%% beta的影响
%ST_T=[];
ST_R=[];
%SC_T=[];
SC_R=[];
%MT_T=[];
MT_R=[];
%MC_T=[];
MC_R=[];
SReal_T=[];
MReal_T=[];
% SCPU_T=[];
% MCPU_T=[];
for beta=0.5:0.05:0.95
    %% ###############################################  理论复杂度
%% 单网格
% ComplexityInTheoryT_Single=ST_Theory(0.7,0.1,beta,0); 
% ST_T(end+1)=ComplexityInTheoryT_Single;
% 
% SC_T(end+1)=ComplexityInTheoryT_Single*Num_H^3;
%% 多网格
% ComplexityInTheoryT_Multi=MT_Theory(0.7,0.1,beta,2,4,5,0);
% MT_T(end+1)=sum(ComplexityInTheoryT_Multi);
% 
% MC_T(end+1)=ComplexityInTheoryT_Multi*Num_h3';
%% ###############################################  实际复杂度
%% 单网格
[SRealTime,SCPUTime,ComplexityrealT_Single]=ST_Reality(0.7,0.1,beta,2,4,5,0);
SReal_T(end+1)=SRealTime;
%SCPU_T(end+1)=SCPUTime;
ST_R(end+1)=ComplexityrealT_Single;
SC_R(end+1)=ComplexityrealT_Single*Num_H^3;
%% 多网格
[MRealTime,MCPUTime,ComplexityrealT_Multi]=MT_Reality(0.7,0.1,beta,2,4,5,0);
MReal_T(end+1)=MRealTime;
%MCPU_T(end+1)=MCPUTime;
MT_R(end+1)=sum(ComplexityrealT_Multi);
MC_R(end+1)=ComplexityrealT_Multi*Num_h3';
end

figure(7);
% plot(0.5:0.05:0.95,ST_T(:),'-or','LineWidth',1)
% hold on
plot(0.5:0.05:0.95,ST_R(:),'-*r','LineWidth',1)
hold on
% plot(0.5:0.05:0.95,MT_T(:),':og','LineWidth',1)
% hold on
plot(0.5:0.05:0.95,MT_R(:),'--og','LineWidth',1)
xlabel('Weight coefficient $\beta$','interpreter','latex','FontSize',14)
ylabel('Iteration complexity','FontSize',14)
grid on
legend_FontSize =legend('$\mathcal{T}^{(sg)}$','$\mathcal{T}^{(mtg)}$','interpreter','latex');
set(legend_FontSize,'FontSize',14)
legend('boxoff')

figure(8);
% plot(0.5:0.05:0.95,SC_T(:),'-or','LineWidth',1)
% hold on
plot(0.5:0.05:0.95,SC_R(:),'-*r','LineWidth',1)
hold on
% plot(0.5:0.05:0.95,MC_T(:),':og','LineWidth',1)
% hold on
plot(0.5:0.05:0.95,MC_R(:),'--og','LineWidth',1)
xlabel('Weight coefficient $\beta$','interpreter','latex','FontSize',14)
ylabel('Query complexity','FontSize',14)
grid on
legend_FontSize =legend('$\mathcal{C}^{(sg)}$','$\mathcal{C}^{(mtg)}$','interpreter','latex');
set(legend_FontSize,'FontSize',14)
legend('boxoff')

figure(9);
plot(0.5:0.05:0.95,SReal_T(:),'-*r','LineWidth',1)
hold on
% plot(0.5:0.05:0.95,SCPU_T(:),'-*r','LineWidth',1)
% hold on
plot(0.5:0.05:0.95,MReal_T(:),'--og','LineWidth',1)
%hold on
% plot(0.5:0.05:0.95,MCPU_T(:),'--og','LineWidth',1)
xlabel('Weight coefficient $\beta$','interpreter','latex','FontSize',14)
ylabel('Running times','FontSize',14)
grid on
legend_FontSize =legend('$T^{(sg)}$','$T^{(mtg)}$','interpreter','latex');
set(legend_FontSize,'FontSize',14)
legend('boxoff')
