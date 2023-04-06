clc
clear all 
clear vars
%% The inputs
% NACA Geometry Inputs
y_mc=0; % Max Camber
x_mc=0; % Position of Max Camber
t_m=0.12; % Max Thickness

%% Flow Parameters
mu = 1.7894 * 10^-5; % Constant
rho = 1.225; % Constant
V_inf = 100; % Constant
alpha=0*(pi/180);
%% Panel Method
n=92;  % number of node points (must be even)
d_theta=(2*pi)/(n-1);

syms f(z)
f11(z) = (y_mc/x_mc^2)*(2*x_mc*z-z^2); % differentiating y_cl before x_mc to sub in sin(theta_cl)
dy_cl_dx_1 = diff(f11,z);

syms f(zz)
f22(zz) = (y_mc/(1-x_mc)^2)*(1-2*x_mc+2*x_mc*zz-zz^2); % differentiating y_cl after x_mc to sub in sin(theta_cl)
dy_cl_dx_2 = diff(f22,zz);


for i=1:n/2
    x(i)=0.5*(1-cos((i-0.5)*d_theta)); %position of nodes
    if x(i)<=x_mc
        y_cl(i)=(y_mc/x_mc^2)*(2*x_mc*x(i)-x(i).^2); % y-coordinate of the any point on the camber line BEFORE x_mc
        dy_cl_dx(i)=double(dy_cl_dx_1(x(i)));    %computing the value of the derivitve
        sin_theta_cl(i)=dy_cl_dx(i)/sqrt(1+(dy_cl_dx(i))^2); % sub in sin eqution
        cos_theta_cl(i)=1/sqrt(1+(dy_cl_dx(i))^2); % sub in cos eqution
    else
        y_cl(i)=(y_mc/(1-x_mc)^2)*(1-2*x_mc+2*x_mc*x(i)-x(i).^2); % y-coordinate of the any point on the camber line AFRER x_mc
        dy_cl_dx(i)=double(dy_cl_dx_2(x(i)));
        sin_theta_cl(i)=dy_cl_dx(i)/sqrt(1+(dy_cl_dx(i))^2);
        cos_theta_cl(i)=1/sqrt(1+(dy_cl_dx(i))^2);
    end
    t=(t_m/2)*(2.969*((x).^0.5)-(1.260*x)-(3.516*x.^2)+(2.843*x.^3)-(1.015*x.^4)); % airfoil thickness at any given x
end

x_u=x-(t.*sin_theta_cl); %x coordinte of upper nodes
x_l=x+(t.*sin_theta_cl);   %x coordinte of lower nodes
y_u=y_cl+(t.*cos_theta_cl);   %y coordinte of upper nodes
y_l=y_cl-(t.*cos_theta_cl);   %y coordinte of lower nodes

for i=1:n/2     % combining upper and lower into one array
    x_N((n/2)+i)=x_u(i);
    x_N((n/2)+1-i)=x_l(i);
    y_N((n/2)+i)=y_u(i);
    y_N((n/2)+1-i)=y_l(i);
end
x_N(1)=1;
x_N(n)=1;
y_N(1)=0;
y_N(n)=0;

for j=1:n-1
    x_c(j)=0.5*(x_N(j)+x_N(j+1));
    y_c(j)=0.5*(y_N(j)+y_N(j+1));
    l_j(j)=sqrt((x_N(j+1)-x_N(j))^2+(y_N(j+1)-y_N(j))^2);
    theta_i(j)=atan2((y_N(j+1)-y_N(j)),(x_N(j+1)-x_N(j)));
    RHS(j)= sin(theta_i(j)- alpha);
end
RHS(n) = 0;

for i=1:n-1

    for j=1:n-1
        if i==j
            C_n1_ij(i,j)=-1;
            C_n2_ij(i,j)=1;
            C_t1_ij(i,j)=0.5*pi;
            C_t2_ij(i,j)=0.5*pi;
        else
            A=-(x_c(i)-x_N(j)).*cos(theta_i(j))-(y_c(i)-y_N(j)).*sin(theta_i(j));
            B=(x_c(i)-x_N(j)).^2+(y_c(i)-y_N(j)).^2;
            C=sin(theta_i(i)-theta_i(j));
            D=cos(theta_i(i)-theta_i(j));
            E=(x_c(i)-x_N(j)).*sin(theta_i(j))-(y_c(i)-y_N(j)).*cos(theta_i(j));
            F=log(1+((l_j(j)).^2+2.*A.*l_j(j))./B);
            G=atan2((E.*l_j(j)),(B+A.*l_j(j)));
            P=(x_c(i)-x_N(j)).*sin(theta_i(i)-2*theta_i(j))+(y_c(i)-y_N(j)).*cos(theta_i(i)-2.*theta_i(j));
            Q=(x_c(i)-x_N(j)).*cos(theta_i(i)-2*theta_i(j))-(y_c(i)-y_N(j)).*sin(theta_i(i)-2.*theta_i(j));
            a=6;
            C_n2_ij(i,j)=D+0.5.*Q.*F./l_j(j)-(A.*C+D.*E).*G./l_j(j);
            C_n1_ij(i,j)=0.5.*D.*F+C.*G-C_n2_ij(i,j);
            C_t2_ij(i,j)=C+0.5.*P.*F./l_j(j)+(A.*D-C.*E).*G/l_j(j);
            C_t1_ij(i,j)=0.5.*C.*F-D.*G-C_t2_ij(i,j);

        end
    end
end

A_n=zeros(n,n);
A_t=zeros(n,n);

for i=1:n-1
    A_n(i,1)=C_n1_ij(i,1);
    A_n(i,n)=C_n2_ij(i,n-1);
    A_t(i,1)=C_t1_ij(i,1);
    A_t(i,n)=C_t2_ij(i,n-1);
    for j=2:n-1
        A_n(i,j)=C_n1_ij(i,j)+C_n2_ij(i,j-1);
        A_t(i,j)=C_t1_ij(i,j)+C_t2_ij(i,j-1);
    end
end
A_n(n,1)=1;
A_n(n,n)=1;
for j=2:n-1
    A_n(n,j)=0; 
end
A_n_inv=inv(A_n);
gama=(A_n_inv*RHS')';


for i=1:n-1
    for j=1:n-1
        K(i,j)=C_t1_ij(i,j)*gama(j)+C_t2_ij(i,j)*gama(j+1);
    end
    KK(i)=sum(K(i,:));
end
for i=1:n-1
        V(i)=cos(theta_i(i)-alpha)+KK(i);
end
V=abs(V);
%% Rearrangment
x_c=[x_c(n/2:n-1) x_c(1:n/2-1)];
y_c=[y_c(n/2:n-1) y_c(1:n/2-1)];
V=[V(n/2:n-1) V(1:n/2-1)];
%% The Derivatives
V_dash=zeros(1,n-1);
V_dash(1:n-2)=diff(V)./diff(x_c);
V_dash(n-1)=(V(1)-V(n-1))/(x_c(1)-x_c(n-1));
V_ddash=((V_dash(2)-V_dash(1))/(x_c(2)-x_c(1)))*V_inf;
V=V*V_inf;
V_dash=V_dash*V_inf;
V_bar=V/V_inf;

%% (a)&(b) Fully Laminar and Taking Transition in consideration
Gamma =zeros(1,n-1);
Gamma(1) = 7.05232; % Assume Gamma from data sheet pg5
Delta1_bar = (3/10)-(Gamma(1)/120);
Delta2_bar = (37/315)-(Gamma(1)/945)-(Gamma(1)^2/9072);
K(1) = (Delta2_bar^2)*Gamma(1);
f1(1) = Delta1_bar/Delta2_bar;
f2(1) = (2+(Gamma(1)/6))*Delta2_bar;
F(1) = 2*f2(1)-4*K(1)-2*K(1)*f1(1); % for iteration (0) F should equal zero
Z_i(1) = K(1)./V_dash(1); 
G(1) = -0.0652*V_ddash(1)./V_dash(1).^2;
Z_i(2) = Z_i(1)+ G(1).*(x_c(2)- x_c(1));
Delta(1) = sqrt((Gamma(1)*mu)/(V_dash(1) * rho));
Delta1(1) = Delta1_bar*Delta(1);
Delta2(1) = Delta2_bar*Delta(1);

for i = 2:n-1
        syms JJ
        eqn = Z_i(i).* V_dash(i) == ((37/315)-(JJ/945)-(JJ^2/9072))^2*JJ; % solve for Gamma
        S = solve(eqn,JJ);
        S=vpa(S);
    for j = 1:5
        if S(j)>-12 && S(j)<12
            Gamma(i)=S(j);
        end
    end
    if Gamma(i)==0
        x_sep=x_c(i-1)
        break
    end
        Delta1_bar = (3/10)-(Gamma(i)/120);
        Delta2_bar = (37/315)-(Gamma(i)/945)-(Gamma(i)^2/9072);
        K(i) = Delta2_bar^2*Gamma(i);
        f1(i) = Delta1_bar/Delta2_bar;
        f2(i) = (2+(Gamma(i)/6))*Delta2_bar;
        F(i) = 2*f2(i)-4*K(i)-2*K(i)*f1(i);
        G(i)=F(i)/V(i);
        Z_i(i+1) = Z_i(i)+G(i).*(x_c(i+1)-x_c(i));
        Delta(i) = sqrt((Gamma(i)*mu)/(V_dash(i)*rho));
        Delta1(i) = Delta1_bar*Delta(i);
        Delta2(i) = Delta2_bar*Delta(i);
        Re(i) = (rho*V(i)*x_c(i)/mu);
        T_w(i) = mu*f2(i)*(V(i)/Delta2(i));
        Cf(i) = T_w(i)/(0.5*rho*V_inf^2);
        d(i-1) = (T_w(i) + T_w(i-1))*(x_c(i)-x_c(i-1));
        
    if  abs(log10(Re(i-1))-( -40.4457+64.8066*f1(i-1)-26.7538*f1(i-1)^2+3.3819*f1(i-1)^3)) < 0.008
        x_tr=x_c(i-1) 
        kk=i;
    end
end
Drag_FL = sum(d)
for ttt=1:n-1
    if x_c(ttt)<=x_tr
        Drag_0(ttt)=d(ttt);
        Delta2_Tran(ttt)=Delta2(ttt);
        C_f_x_Tran(ttt)=Cf(ttt);
        Delta1_Tran(ttt)=Delta1(ttt);
    else
        Drag_L=sum(Drag_0)
        break
    end
end
%% Turbulant
Delta2_T=zeros(1,n-1);
DeltaT=zeros(1,n-1);
Delta2_T_bar=zeros(1,n-1);
Red2=zeros(1,n-1);
c_f_x_t=zeros(1,n-1);
d_Delta2_bar_T_dx=zeros(1,n-1);
H=zeros(1,n-1);
H1=zeros(1,n-1);
num1=zeros(1,n-1);
num2=zeros(1,n-1);
dem1=zeros(1,n-1);
dem2=zeros(1,n-1);
cc=zeros(1,n-1);
T_w_T=zeros(1,n-1);
d_T=zeros(1,n-1);
ReL=rho*V_inf/mu;
for i=kk-1:n-1
    if x_c(i)==x_tr
        Delta2_T(i)=Delta2(i);
        DeltaT(i)=1.4*Delta(i);
        Delta2_T_bar(i)=Delta2_T(i)/DeltaT(i);
        H(i)=fzero(@(AA) ((DeltaT(i)/Delta2_T(i))-AA-3.3-1.5501*(AA-0.6718)^(-3.064)),1.2);
        Red2(i)=rho*V(i)*Delta2_T(i)/mu;
        c_f_x_t(i)=0.246*10^(-0.678*H(i))*(Red2(i)^(-0.268));
        
        % for i+1
        d_Delta2_bar_T_dx(i)=c_f_x_t(i)/2-(Delta2_T_bar(i)/V(i))*V_dash(i)*(H(i)+2);
        Delta2_T_bar(i+1)=Delta2_T_bar(i)+d_Delta2_bar_T_dx(i)*(x_c(i+1)-x_c(i));
        H1(i+1)=3.3+0.8234*(H(i)-1.1)^(-1.287);
       
        if H1(i+1) <3.3
            H(i+1)=3;
        elseif H1(i+1)<=5.3 && H1(i+1)>3.3
            H(i+1)=0.6778+1.1536*(H1(i+1)-3.3)^(-0.326);
        elseif H1(i+1)>5.3
            H(i+1)=1.1+0.86*(H1(i+1)-3.3)^(-0.777);
        end
        if  H(i+1)>=1.8 && H(i+1)<=2.8
            x_sep_T=x_c(i);
            break
        end
        T_w_T(i)=0.5*rho*V_inf^2*c_f_x_t(i);
        d_T(i-1)=(T_w_T(i) + T_w_T(i-1))*(x_c(i)-x_c(i-1));
    elseif x_c(i)>x_tr
        Red2(i)=V_bar(i)*Delta2_T_bar(i)*ReL;
        c_f_x_t(i)=0.246*10^(-0.678*H(i))*(Red2(i)^(-0.268));
        T_w_T(i)=0.5*rho*V_inf^2*c_f_x_t(i)*a;
        d_T(i-1)=(T_w_T(i) + T_w_T(i-1))*(x_c(i)-x_c(i-1));
        % for i+1
        if i==n/2
            break
        end
        syms ddd
        Delta2_T_bar(i+1)=Delta2_T_bar(i)+((0.246*10^(-0.678*H(i))*((V_bar(i)*Delta2_T_bar(i)*ReL)^(-0.268)))/2-(Delta2_T_bar(i)/V(i))*V_dash(i)*(H(i)+2))*(x_c(i+1)-x_c(i));
        num1(i)=V_bar(i+1)*0.0306*(H1(i)-3)^(-0.6169);
        num2(i)=(V_bar(i+1)*Delta2_T_bar(i+1)*H1(i))/(x_c(i+1)-x_c(i));
        dem1(i)=((V_bar(i+1)*Delta2_T_bar(i+1))-(V_bar(i)*Delta2_T_bar(i)))/(x_c(i+1)-x_c(i));
        dem2(i)=(V_bar(i+1)*Delta2_T_bar(i+1))/(x_c(i+1)-x_c(i));
        H1(i+1)=(num1(i)+num2(i))/(dem1(i)+dem2(i));

        if H1(i+1) <3.3
            H(i+1)=3;
        elseif H1(i+1)<=5.3 && H1(i+1)>3.3
            H(i+1)=0.6778+1.1536*(H1(i+1)-3.3)^(-0.326);
        elseif H1(i+1)>5.3
            H(i+1)=1.1+0.86*(H1(i+1)-3.3)^(-0.777);
        end
        if H(i+1)>=1.8 && H(i+1)<=2.8
            x_sep_T=x_c(i)
            break
        end
    end
end
Drag_T=sum(d_T)
Drag_Total=Drag_T+Drag_L
CD=Drag_Total/(0.5*rho*V_inf^2)
%% Figure
% Fully Laminar
for i=1:n/2
    if x_c(i)<=x_sep
        x_fl(i)=x_c(i);
    end
    if x_c(i)<=x_c(n/2)
        x_turb(i)=x_c(i);
    end
    if x_c(i)>=x_tr
        x_fx(i)=x_c(i);
        
    end
    if x_c(i)<=x_tr
        x_l_t(i)=x_c(i); 
    end
end
for i=17:n/2
    V_delta2(i)=V(i);
end
V_delta2=nonzeros(V_delta2)';
Red2delta2=nonzeros(Red2)';
Delta2_Turb=nonzeros(Red2delta2.*(mu./(rho*V_delta2)));
Delta2_Turb(30)=0;
Delta2_Turb=nonzeros(Delta2_Turb)';
H1_n=nonzeros(H1)';
H1_nn=H1_n;
H1_nn(29)=0;
H1_nn(28)=0;
H1_nn=nonzeros(H1_nn)';
H_n=nonzeros(H)';
H_n(30)=0;
H_n=nonzeros(H_n)';
Delta1_Turb=H_n.*Delta2_Turb;
H_n(29)=0;
H_n(28)=0;
H_n=nonzeros(H_n)';


Delta1_Total=[Delta1_Tran Delta1_Turb];
Delta2_Total=[Delta2_Tran Delta2_Turb];
c_f_x_Turb=nonzeros(c_f_x_t)';
c_f_x_Turb(1)=C_f_x_Tran(17);
x_Fx=nonzeros(x_fx)';

figure
plot(x_N,y_N,"Linewidth",2)
grid on
xlim([-0.5 1.5])
ylim([-0.5 0.5])
title('Airfoil Profile','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')

figure
plot(x_c,V_bar,"Linewidth",2)
grid on
title('Velocity Profile','Interpreter','latex')
xlabel('$\bar{x}$','Interpreter','latex')
ylabel('$\bar{U}$','Interpreter','latex')

figure
plot(x_fl,Delta1,"Linewidth",2)
grid on
axis tight
title('Fully Laminar ${\delta_1}$ vs $\bar{x}$ ','Interpreter','latex')
xlabel('$\bar{x}$','Interpreter','latex')
ylabel('${\delta_1}$','Interpreter','latex')

figure
plot(x_fl,Delta2,"Linewidth",2)
grid on
axis tight
title('Fully Laminar ${\delta_2}$ vs $\bar{x}$','Interpreter','latex')
xlabel('$\bar{x}$','Interpreter','latex')
ylabel('${\delta_2}$','Interpreter','latex')

figure
plot(x_fl,Cf,"Linewidth",2)
grid on
axis tight
title('Fully Laminar $C_{fx}$ vs $\bar{x}$','Interpreter','latex')
xlabel('$\bar{x}$','Interpreter','latex')
ylabel('$C_{fx}$','Interpreter','latex')

figure
plot(x_turb,Delta1_Total,"Linewidth",2)
grid on

title('Laminar-Transition-Turbulant ${\delta_1}$ vs $\bar{x}$','Interpreter','latex')
xlabel('$\bar{x}$','Interpreter','latex')
ylabel('${\delta_1}$','Interpreter','latex')

figure
plot(x_turb,Delta2_Total,"Linewidth",2)
grid on
axis tight
title('Laminar-Transition-Turbulant ${\delta_2}$ vs $\bar{x}$','Interpreter','latex')
xlabel('$\bar{x}$','Interpreter','latex')
ylabel('${\delta_2}$','Interpreter','latex')

figure
plot(x_l_t,C_f_x_Tran,x_Fx,c_f_x_Turb,"Linewidth",2)
grid on
axis tight
title('Laminar-Transition-Turbulant $C_{fx}$ vs $\bar{x}$','Interpreter','latex')
xlabel('$\bar{x}$','Interpreter','latex')
ylabel('$C_{fx}$','Interpreter','latex')


figure
plot(H_n,H1_nn,"Linewidth",2)
grid on
axis tight
title('Laminar-Transition-Turbulant $H_{1}$ vs $H$','Interpreter','latex')
xlabel('$H$','Interpreter','latex')
ylabel('$H_{1}$','Interpreter','latex')



