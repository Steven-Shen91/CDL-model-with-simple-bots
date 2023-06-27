%%replicator dynamics
clear;
clc;
r= 0.2; %%temptation to defect 1+r;
K1=-5:1:3;%%imitation strength
delta=0.3; %%loner's payoff
% mu=power(10,-3);
mu=0;%%mutation
rho=0:0.2:1.0;%%fractions of bots.
type=3;%% 0,1,2,3 represent Z_C,Z_D,Z_L as wll as Z_ALL;
MCS=500;
tspan=[0:1:MCS]; 
% init_f=(0.01:0.01:0.99)';
% y0=[(1-init_f)/3 (1-init_f)/3 (1-init_f)/3 init_f];
% y0=[0.96 0.01 0.01 0.01];
% y0=[0.01 0.96 0.01 0.01];
% y0=[0.01 0.01 0.96 0.01];
% y0=[0.01 0.01 0.01 0.96];
try_time=100;
y0_temp=rand(3,try_time);
y0=y0_temp./sum(y0_temp,1);%%initial fraction of actors
try_time=size(y0,2);
Y_final=[];
% y0=[0.1,0.8,0.1];
frac_C_temp=zeros(length(rho),length(K1));
frac_D_temp=zeros(length(rho),length(K1));
frac_L_temp=zeros(length(rho),length(K1));
% for i=1:1:length(r)
for i=1:1:length(K1)
    for j=1:1:length(rho)
        tempC=0;
        tempD=0;
        tempL=0;
        for k=1:1:size(y0,2)
            options=odeset('RelTol',1e-12,'AbsTol',1e-12);
            [T,Y]=ode45(@(t,y) odefun(t,y,r,rho(j),delta,type,mu,power(10,K1(i))),tspan,y0(:,k),options); 
%             Y_final=vertcat(Y_final,Y);
%             tempC=round(Y(end,1),3);
%             tempD=round(Y(end,2),3);
%             tempCP=round(Y(end,3),3);
%             tempDP=round(Y(end,4),3);
        tempC=tempC+round(Y(end,1),3);
        tempD=tempD+round(Y(end,2),3);
        tempL=tempL+round(Y(end,3),3);
%         ternaryplot(Y(:,1),Y(:,2),Y(:,3))
% figure(2)
%         plot(T,Y(:,1),'k',T,Y(:,2),'r',T,Y(:,3),'g');
% % %         legend('C','D','L')
%         hold on;
% %             frac_C_temp(k,i)=tempC;
%             frac_D_temp(k,i)=tempD;
%             frac_CP_temp(k,i)=tempCP;
%             frac_DP_temp(k,i)=tempDP;
        end
            frac_C_temp(j,i)=tempC/try_time;
            frac_D_temp(j,i)=tempD/try_time;
            frac_L_temp(j,i)=tempL/try_time;
    end
end

frac_C=round(frac_C_temp,2);
frac_D=round(frac_D_temp,2);
frac_L=round(frac_L_temp,2);

figure(1)
subplot(1,3,1)
C=pcolor(r,rho,frac_C);
colormap(jet)
set(C,'linestyle','none');
xlabel('Dilemma strength')
ylabel('fraction of bots')
title('Cooperation')

subplot(1,3,2)
C=pcolor(r,rho,frac_D);
colormap(jet)
set(C,'linestyle','none');
xlabel('Dilemma strength')
ylabel('fraction of bots')
title('Defection')
subplot(1,3,3)
C=pcolor(r,rho,frac_L);
colormap(jet)
set(C,'linestyle','none');
xlabel('Dilemma strength')
ylabel('fraction of bots')
title('Loners')