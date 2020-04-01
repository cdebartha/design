clc; clear all;

%% 
load('loads.mat') ;
a = load('n63415-il.dat') ;
c = 1.6190*3.28084 ;
z = -a(:,1)'*c*0.01;
y = a(:,2)'*c*0.01;
spar1_ytop = 0.1517*3.28084;
spar1_ybot = -0.08855*3.28084;
spar1_z    = -0.3*c ;
spar2_ytop = 0.09546*3.28084;
spar2_ybot = -0.031975*3.28084;
spar2_z    = -0.7*c ;
b1         = 0.1*(spar1_ytop - spar1_ybot);
b2         = 0.1*(spar2_ytop - spar2_ybot);

%% Thicknesses and Modulus
tw = 0.002*3.28084; %ft
tf = tw/2; % ft
tsc = 2*0.001*3.28084; %ft
hw1 = (spar1_ytop-spar1_ybot-2*(tf+tsc));
hw2 = (spar2_ytop-spar2_ybot-2*(tf+tsc));
E_al = 1.5351*10^9*32.1522; % lbf/ft^2 
E_steel = 4.2*10^9*32.1522; % lbf/ft^2 
E_star = E_steel/E_al;
ybar = [(spar1_ytop+spar1_ybot)/2,(spar2_ytop+spar2_ybot)/2];
zbar = [spar1_z,spar2_z];

%% for skin
t_skin = 1e-3*3.28084 ; %         ft
E_skin = E_al;              % Skin Material 2024T3 Aluminium
E_ref = E_al ;
E_star_skin = E_skin/E_ref ;
% z = wingdata(:,1)*chord/100;
% y = wingdata(:,2)*chord/100;

for j=1 : (length(y)-1)
    s = sqrt((y(j)-y(j+1))^2 + ((z(j)-z(j+1))^2));
    EdA_star_skin(j) = s*t_skin*(E_star_skin);
    ydA_star(j) = EdA_star_skin(j)*0.5*(y(j)+y(j+1));
    zdA_star(j) = EdA_star_skin(j)*0.5*(z(j)+z(j+1));
%     y2dA_star(j) = ydA_star(j)*0.5*(y(j)+y(j+1));
%     z2dA_star(j) = zdA_star(j)*0.5*(z(j)+z(j+1));
%     yzdA_star(j) = zdA_star(j)*0.5*(y(j)+y(j+1));
end

Sigma_EA_skin = sum(EdA_star_skin);
y_cent_skin = sum(ydA_star)/Sigma_EA_skin;
z_cent_skin = sum(zdA_star)/Sigma_EA_skin;
%% Areas of individual parts of spar
A1 = [b1*tsc,b1*tsc,b1*tf,b1*tf,tw*hw1];
A2 = [b2*tsc,b2*tsc,b2*tf,b2*tf,tw*hw2];
E = [E_steel,E_steel,E_al,E_al,E_al];
E_star_spar = (1/E_al).*E ;

%% Finding spar centroid
EA_star_i = [dot(A1,E_star_spar),dot(A2,E_star_spar)];
Sigma_EA_star_spar = sum(EA_star_i) ;
y_cent_spar = dot(EA_star_i,ybar)/Sigma_EA_star_spar ;
z_cent_spar = dot(EA_star_i,zbar)/Sigma_EA_star_spar;

%% Centroid of spar and Skin
y_cent = (y_cent_spar*Sigma_EA_star_spar + y_cent_skin*Sigma_EA_skin)/(Sigma_EA_star_spar+Sigma_EA_skin ) ;
z_cent = (z_cent_spar*Sigma_EA_star_spar + z_cent_skin*Sigma_EA_skin)/(Sigma_EA_star_spar+Sigma_EA_skin) ;


%% Finding Moment of Inertia for spars
Izz1 = 2*((E_star*(b1*tsc^3/12))+(A1(1)*E_star*((tsc+2*tf+hw1)/2)^2))+ 2*(((b1*tf^3/12))+(A1(3)*((tw/4)+(hw1/2))^2))+ (tw*hw1^3/12);
Iyy1 = 2*E_star*(tsc*b1^3/12)+2*(tf*b1^3/12)+ (hw1*tw^3/12);
Izz2 = 2*((E_star*(b2*tsc^3/12))+(A2(1)*E_star*((tsc+2*tf+hw1)/2)^2))+ 2*(((b2*tf^3/12))+(A2(3)*((tw/4)+(hw2/2))^2))+ (tw*hw2^3/12);
Iyy2 = 2*E_star*(tsc*b2^3/12)+2*(tf*b2^3/12)+ (hw2*tw^3/12);

%% MOI of Skin about neutral axis
z2_wrt_cent = ((z(1:length(z)-1) - z_cent*linspace(1,1,length(z)-1)).^2) ;
y2_wrt_cent = ((y(1:length(z)-1) - y_cent*linspace(1,1,length(z)-1)).^2);
yz_wrt_cent = ( ((z(1:length(z)-1) - z_cent*linspace(1,1,length(z)-1))).*(y(1:length(y)-1) - y_cent*linspace(1,1,length(y)-1)) ) ;
Iyy_skin = dot (EdA_star_skin, z2_wrt_cent) ;
Izz_skin = dot (EdA_star_skin, y2_wrt_cent) ;
Iyz_skin = dot (EdA_star_skin, yz_wrt_cent);

%% MOI spars about neutral axis
Iyy_spar = Iyy1+ (EA_star_i(1))*(zbar(1)-z_cent)^2 + Iyy2+ (EA_star_i(2))*(zbar(2)-z_cent)^2  ;
Izz_spar = Izz1+ (EA_star_i(1))*(ybar(1)-y_cent)^2 + Izz2+ (EA_star_i(2))*(ybar(2)-y_cent)^2  ;    
Iyz_spar = ((EA_star_i(1))*(ybar(1)-y_cent)*(zbar(1)-z_cent) + (EA_star_i(2))*(zbar(2)-z_cent)*(ybar(2)-y_cent)) ;

%% Total MOI skin + spars about neutral axis
Iyy_total= Iyy_spar + Iyy_skin;
Izz_total= Izz_spar + Izz_skin;    
Iyz_total= Iyz_spar + Iyz_skin;


%% Finding V_xx, W_xx
A = [-E_al*Iyz_total,-E_al*Iyy_total ; E_al*Izz_total,E_al*Iyz_total];
M = [M_y_due_to_D(1),M_z_due_to_L(1)];
%M=[855620,7376500]
v = A\M';

%% Finding Sigma_xx
F_S = 1.5 ;
epsilon_xx= (-(y-y_cent)*v(1)-(z-z_cent)*v(2));  %conversion from meter to feet
sigma_xx= F_S*E_al*epsilon_xx;  %% in lbf/ft^2
sigma_yield_pos= 32.1522.*[6804000,6804000]; %32.1522 %lbf/ft^2 
sigma_yield_neg = -32.1522.*[6804000,6804000]; %lbf/ft^2 
z1=[-5.3117,0];


width_val = 1.5 ;
figure(5)
plot(-z,sigma_xx,'k','linewidth',width_val)
hold on ; 
plot(-z1,sigma_yield_pos,-z1,sigma_yield_neg,'r','linewidth',width_val);
title('Stress distribution over the Root section');
xlabel('x (ft)');
ylabel('\bf \sigma_{xx} (lbf/ft^2)');
grid on
%% Plotting Sigm_xx
figure(6)
plot(-z,sigma_xx,'-*', 'Linewidth', 1.2)
hold on ; 
plot(-z1,sigma_yield_pos,-z1,sigma_yield_neg,'Linewidth', 1.2);grid on ;

ax = gca;
ax.FontSize = 14;
ax.XAxis.LineWidth = 1.2;
ax.YAxis.LineWidth = 1.2;
set(gcf, 'Position',  [100, 100, 900, 600]);
title('\sigma_{xx} distribution on the surface of the airfoil');
xlabel('x(ft)');
ylabel('\sigma_{xx}(in lbf/ft^2)');


