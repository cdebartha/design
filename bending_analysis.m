clc; clear all;

%% 
load('loads.mat') ;
a = load('n63415-il.dat') ;
c = 1.6190*3.28084 ;
z = -a(:,1)*c*0.01;
y = a(:,2)*c*0.01;
spar1_ytop = 0.1517*3.28084;
spar1_ybot = -0.08855*3.28084;
spar1_z = -0.3*c*3.28084;
spar2_ytop = 0.09546*3.28084;
spar2_ybot = -0.031975*3.28084;
spar2_z = -0.7*c;
b1 = 0.10*(spar1_ytop - spar1_ybot);
b2 = 0.10*(spar2_ytop - spar2_ybot);

%% Thicknesses and Modulus
tw = 0.004*3.28084; %ft
tf = tw/2; % ft
tsc = 5*0.002*3.28084; %ft
hw1 = (spar1_ytop-spar1_ybot-2*(tf+tsc));
hw2 = (spar2_ytop-spar2_ybot-2*(tf+tsc));
E_al = 1.5351*10^9*32.1522; % lbf/ft^2 
E_steel = 4.2*10^9*32.1522; % lbf/ft^2 
E_star = E_steel/E_al;
ybar = [(spar1_ytop+spar1_ybot)/2,(spar2_ytop+spar2_ybot)/2];
zbar = [spar1_z,spar2_z];

%% Areas of individual parts
A1 = [b1*tsc,b1*tsc,b1*tf,b1*tf,tw*hw1];
A2 = [b2*tsc,b2*tsc,b2*tf,b2*tf,tw*hw2];
E = [E_steel,E_steel,E_al,E_al,E_al];
E_star_spar = (1/E_al).*E ;

%% Finding spar centroid
A_star_i = [dot(A1,E_star_spar),dot(A2,E_star_spar)];
Sigma_A_star = sum(A_star_i) ;
y_cent = dot(A_star_i,ybar)/Sigma_A_star ;
z_cent = dot(A_star_i,zbar)/Sigma_A_star;

%% Finding Moment of Inertia for spars
Izz1 = 2*((E_star*(b1*tsc^3/12))+(A1(1)*E_star*((tsc+2*tf+hw1)/2)^2))+ 2*(((b1*tf^3/12))+(A1(3)*((tw/4)+(hw1/2))^2))+ (tw*hw1^3/12);
Iyy1 = 2*E_star*(tsc*b1^3/12)+2*(tf*b1^3/12)+ (hw1*tw^3/12);
Izz2 = 2*((E_star*(b2*tsc^3/12))+(A2(1)*E_star*((tsc+2*tf+hw1)/2)^2))+ 2*(((b2*tf^3/12))+(A2(3)*((tw/4)+(hw2/2))^2))+ (tw*hw2^3/12);
Iyy2 = 2*E_star*(tsc*b2^3/12)+2*(tf*b2^3/12)+ (hw2*tw^3/12);

%% Total MOI stringers + spars
Iyy_total = Iyy1+ (A_star_i(1))*(zbar(1)-z_cent)^2 + Iyy2+ (A_star_i(2))*(zbar(2)-z_cent)^2  ;
Izz_total = Izz1+ (A_star_i(1))*(ybar(1)-y_cent)^2 + Izz2+ (A_star_i(2))*(ybar(2)-y_cent)^2  ;    
Iyz_total = ((A_star_i(1))*(ybar(1)-y_cent)*(zbar(1)-z_cent) + (A_star_i(2))*(zbar(2)-z_cent)*(ybar(2)-y_cent)) ;

%% Finding V_xx, W_xx
A = [-E_al*Iyz_total,-E_al*Iyy_total ; E_al*Izz_total,E_al*Iyz_total];
M = [M_y_due_to_D(1),M_z_due_to_L(1)];
%M=[855620,7376500]
v = A\M';

%% Finding Sigma_xx
F_S = 1.5 ;
epsilon_xx= (-(y-y_cent)*v(1)-(z-z_cent)*v(2));  %conversion from meter to feet
sigma_xx= F_S*E_al*epsilon_xx;  %% in lbf/ft^2
sigma_yield= 32.1522.*[6804000,6804000]; %32.1522 %lbf/ft^2 
z1=[-5.3117,0];

%% Plotting Sigma_xx
figure(5)
plot(-z,sigma_xx,'-*', 'Linewidth', 1.2)
hold on ; 
plot(-z1,sigma_yield, 'Linewidth', 1.2);
grid on ;
ax = gca;
ax.FontSize = 14;
ax.XAxis.LineWidth = 1.2;
ax.YAxis.LineWidth = 1.2;
set(gcf, 'Position',  [100, 100, 900, 600]);
title('\sigma_{xx} distribution on the surface of the airfoil');
xlabel('x(ft)');
ylabel('\sigma_{xx}(in lbf/ft^2)');

