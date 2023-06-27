function dy=odefun(t,y,r,rho,delta,type,mu,k)
if type==0 %%Zealots C
    dy=zeros(3,1);
    %% average payoff for C, D AND L
    P_C=(y(1)-y(2)*r+y(3)*delta+rho)/(1+rho);
    P_D=((1+r)*(y(1)+rho)+delta*y(3))/(1+rho);
    P_L=delta;

    %% probabilities for strategy transition.
    F_D_C=1/(1+exp(-k*(P_C-P_D)));
    F_L_C=1/(1+exp(-k*(P_C-P_L)));
    %% C replaces D and L
    F_C_D=1/(1+exp(-k*(P_D-P_C)));
    F_L_D=1/(1+exp(-k*(P_D-P_L)));
    %% D replaces C and L
    F_C_L=1/(1+exp(-k*(P_L-P_C)));
    F_D_L=1/(1+exp(-k*(P_L-P_D)));
    %% L replaces D and C

    dy(1)=((y(1)+rho)*y(2)*F_D_C+(y(1)+rho)*y(3)*F_L_C-y(1)*y(2)*F_C_D-y(1)*y(3)*F_C_L)*2/(1+rho);
    dy(2)=(y(1)*y(2)*F_C_D+y(2)*y(3)*F_L_D-(y(1)+rho)*y(2)*F_D_C-y(2)*y(3)*F_D_L)*2/(1+rho);
    dy(3)=(y(1)*y(3)*F_C_L+y(2)*y(3)*F_D_L-(y(1)+rho)*y(3)*F_L_C-y(2)*y(3)*F_L_D)*2/(1+rho);
  
elseif type==1 %%Zealots D
    dy=zeros(3,1);
        P_C=(y(1)-(y(2)+rho)*r+y(3)*delta)/(1+rho);
        P_D=((1+r)*y(1)+delta*y(3))/(1+rho);
        P_L=delta;

        %% probabilities for strategy transition.
        F_D_C=1/(1+exp(-k*(P_C-P_D)));
        F_L_C=1/(1+exp(-k*(P_C-P_L)));
        %% C replaces D and L
        F_C_D=1/(1+exp(-k*(P_D-P_C)));
        F_L_D=1/(1+exp(-k*(P_D-P_L)));
        %% D replaces C and L
        F_C_L=1/(1+exp(-k*(P_L-P_C)));
        F_D_L=1/(1+exp(-k*(P_L-P_D)));
        %% L replaces D and C

        dy(1)=(y(1)*y(2)*F_D_C+y(1)*y(3)*F_L_C-y(1)*(y(2)+rho)*F_C_D-y(1)*y(3)*F_C_L)*2/(1+rho);
        dy(2)=(y(1)*(y(2)+rho)*F_C_D+(y(2)+rho)*y(3)*F_L_D-y(1)*y(2)*F_D_C-y(2)*y(3)*F_D_L)*2/(1+rho);
        dy(3)=(y(1)*y(3)*F_C_L+y(2)*y(3)*F_D_L-y(1)*y(3)*F_L_C-(y(2)+rho)*y(3)*F_L_D)*2/(1+rho);
elseif type==2   %%Zealots L
     dy=zeros(3,1);
        P_C=(y(1)-y(2)*r+(y(3)+rho)*delta)/(1+rho);
        P_D=((1+r)*y(1)+delta*(y(3)+rho))/(1+rho);
        P_L=delta;

        %% probabilities for strategy transition.
        F_D_C=1/(1+exp(-k*(P_C-P_D)));
        F_L_C=1/(1+exp(-k*(P_C-P_L)));
        %% C replaces D and L
        F_C_D=1/(1+exp(-k*(P_D-P_C)));
        F_L_D=1/(1+exp(-k*(P_D-P_L)));
        %% D replaces C and L
        F_C_L=1/(1+exp(-k*(P_L-P_C)));
        F_D_L=1/(1+exp(-k*(P_L-P_D)));
        %% L replaces D and C

        dy(1)=(y(1)*y(2)*F_D_C+y(1)*y(3)*F_L_C-y(1)*y(2)*F_C_D-y(1)*(y(3)+rho)*F_C_L)*(1-mu)*2/(1+rho)-mu*y(1)+y(2)*mu/2+y(3)*mu/2;
        dy(2)=(y(1)*y(2)*F_C_D+y(2)*y(3)*F_L_D-y(1)*y(2)*F_D_C-y(2)*(y(3)+rho)*F_D_L)*(1-mu)*2/(1+rho)-mu*y(2)+y(1)*mu/2+y(3)*mu/2;
        dy(3)=(y(1)*(y(3)+rho)*F_C_L+y(2)*(y(3)+rho)*F_D_L-y(1)*y(3)*F_L_C-y(2)*y(3)*F_L_D)*(1-mu)*2/(1+rho)-mu*y(3)+y(1)*mu/2+y(2)*mu/2;
       
else %%Zealots C/D/L
    dy=zeros(3,1);
        P_C=(y(1)+rho/3-(y(2)+rho/3)*r+(y(3)+rho/3)*delta)/(1+rho);
        P_D=((1+r)*(y(1)+rho/3)+delta*(y(3)+rho/3))/(1+rho);
        P_L=delta;

        %% probabilities for strategy transition.
        F_D_C=1/(1+exp(-k*(P_C-P_D)));
        F_L_C=1/(1+exp(-k*(P_C-P_L)));
        %% CN replaces D and AP
        F_C_D=1/(1+exp(-k*(P_D-P_C)));
        F_L_D=1/(1+exp(-k*(P_D-P_L)));
        %% D replaces CN and AP
        F_C_L=1/(1+exp(-k*(P_L-P_C)));
        F_D_L=1/(1+exp(-k*(P_L-P_D)));
        %% AP replaces D and CN

        dy(1)=((y(1)+rho/3)*y(2)*F_D_C+(y(1)+rho/3)*y(3)*F_L_C-y(1)*(y(2)+rho/3)*F_C_D-y(1)*(y(3)+rho/3)*F_C_L)*2/(1+rho);
        dy(2)=(y(1)*(y(2)+rho/3)*F_C_D+(y(2)+rho/3)*y(3)*F_L_D-(y(1)+rho/3)*y(2)*F_D_C-y(2)*(y(3)+rho/3)*F_D_L)*2/(1+rho);
        dy(3)=(y(1)*(y(3)+rho/3)*F_C_L+y(2)*(y(3)+rho/3)*F_D_L-(y(1)+rho/3)*y(3)*F_L_C-(y(2)+rho/3)*y(3)*F_L_D)*2/(1+rho);
     
end
end