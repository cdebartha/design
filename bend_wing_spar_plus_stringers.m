clc; clear all ;
%%
load('loads.mat') ;
a=load('n63415-il.dat') ;
c=1.6190 ;
z=-a(:,1)*c*3.28084*0.01;
y= a(:,2)*c*3.28084*0.01;
spar1_ytop=0.1517*3.28084;
spar1_ybot=-0.08855*3.28084;
spar1_z=-0.3*c*3.28084;
spar2_ytop=0.09546*3.28084;
spar2_ybot=-0.031975*3.28084;
spar2_z=-0.7*c*3.28084;
b1=(spar1_ytop-spar1_ybot);
b2=(spar2_ytop-spar2_ybot);
%% thicknesses and Modulus
tw=0.004*3.28084;
tf=tw/2; 
tsc=0.002*3.28084;
hw1=(spar1_ytop-spar1_ybot-2*(tf+tsc));
hw2=(spar2_ytop-spar2_ybot-2*(tf+tsc));
E_al= 1.5351*10^9*32.1522;
E_steel= 4.2*10^9*32.1522;
E_star=E_steel/E_al;
ybar=[(spar1_ytop+spar1_ybot)/2,(spar2_ytop+spar2_ybot)/2];
zbar=[spar1_z,spar2_z];
%% Areas of individual parts
A1=[b1*tsc,b1*tsc,b1*tf,b1*tf,tw*hw1];
A2=[b2*tsc,b2*tsc,b2*tf,b2*tf,tw*hw2];
E=[E_steel,E_steel,E_al,E_al,E_al];
E_star_spar= (1/E_al).*E ;

%% finding spar centroid
A_star_i=[dot(A1,E_star_spar),dot(A2,E_star_spar)];
Sigma_A_star = sum(A_star_i) ;
y_cent_spar= dot(A_star_i,ybar)/Sigma_A_star ;
z_cent_spar= dot(A_star_i,zbar)/Sigma_A_star;

%% adding stigers
no_of_stringer = 4 ;
z_stinger=[z(12),z(9),z(40),z(43)];
y_stinger=[y(12),y(9),y(40),y(43)];

%% stringer centroid 
A_stringer = 0.0052 ; % in ft^2
E_star_stgr = E_steel/E_al ;
EA_stgr = A_stringer*E_star_stgr ;
Sigma_EA_stgr = no_of_stringer*EA_stgr ;
y_cent_stgr = EA_stgr*(sum(y_stinger))/(Sigma_EA_stgr) ;
z_cent_stgr = EA_stgr*(sum(z_stinger))/(Sigma_EA_stgr) ;

%% centroid of spar and stringer
y_cent = (y_cent_spar*Sigma_A_star + y_cent_stgr*Sigma_EA_stgr)/(Sigma_A_star+Sigma_EA_stgr) ;
z_cent = (z_cent_spar*Sigma_A_star + z_cent_stgr*Sigma_EA_stgr)/(Sigma_A_star+Sigma_EA_stgr) ;


%% Finding Moment of Inertia for spars
Izz1=2*((E_star*(b1*tsc^3/12))+(A1(1)*E_star*((tsc+2*tf+hw1)/2)^2))+ 2*(((b1*tf^3/12))+(A1(3)*((tw/4)+(hw1/2))^2))+ (tw*hw1^3/12);
Iyy1=2*E_star*(tsc*b1^3/12)+2*(tf*b1^3/12)+ (hw1*tw^3/12);
Izz2=2*((E_star*(b2*tsc^3/12))+(A2(1)*E_star*((tsc+2*tf+hw1)/2)^2))+ 2*(((b2*tf^3/12))+(A2(3)*((tw/4)+(hw2/2))^2))+ (tw*hw2^3/12);
Iyy2=2*E_star*(tsc*b2^3/12)+2*(tf*b2^3/12)+ (hw2*tw^3/12);

%% MOI for stinger about neutral axis
Iyy_stingers = EA_stgr*( norm(z_stinger - z_cent*linspace(1,1,length(z_stinger))) )^2 ;
Izz_stingers = EA_stgr*( norm(y_stinger - y_cent*linspace(1,1,length(y_stinger))) )^2 ;
Iyz_stringer = EA_stgr*( dot((z_stinger - z_cent*linspace(1,1,length(z_stinger))),(y_stinger - y_cent*linspace(1,1,length(y_stinger)))) ) ;

%% total MOI stringers + spars
Iyy_total=Iyy1+ (A_star_i(1))*(zbar(1)-z_cent)^2 + Iyy2+ (A_star_i(2))*(zbar(2)-z_cent)^2  + Iyy_stingers ;
Izz_total=Izz1+ (A_star_i(1))*(ybar(1)-y_cent)^2 + Izz2+ (A_star_i(2))*(ybar(2)-y_cent)^2  + Izz_stingers ;    
Iyz_total=((A_star_i(1))*(ybar(1)-y_cent)*(zbar(1)-z_cent) + (A_star_i(2))*(zbar(2)-z_cent)*(ybar(2)-y_cent))  + Iyz_stringer;
%% finding V_xx, W_xx
A=[-E_al*Iyz_total,-E_al*Iyy_total ; E_al*Izz_total,E_al*Iyz_total];
M=[M_y_due_to_D(1),M_z_due_to_L(1)];

v=A\M';
%% finding sigma_xx
F_S = 1.5 ;
epsilon_xx= (-(y-y_cent)*v(1)-(z-z_cent)*v(2));  %conversion from meter to feet
sigma_xx= F_S*E_al*epsilon_xx;
sigma_yield= 32.1522.*[6804000,6804000]; %32.1522 %psi
z1=[-5.3117,0];
figure(2)
plot(-z,sigma_xx)
hold on ; plot(-z1,sigma_yield);


