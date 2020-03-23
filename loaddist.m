clc, clear all,close all ;

%%
rho_h = 0.963*(2.205/(3.28084)^3) ; % in slug/ft^3
rho_sl = 1.225*(2.205/(3.28084)^3) ;
V_cruise = 75*3.28084; % in ft/s
g = 9.8*3.28084; % in ft/s^2
w = 1800*0.0685218; % in slug
w_cruise = 0.97*0.985*w;
wbys = 82.0114*0.00636588; % in slug/ft^2
lambda = 1 ;

%%
pressure = 10.9174 ;% in lb/ft^2
density = 18.685e-4 ; % in slug/ft^3
viscosity = 3.58e-7  ;% in lb.s/ft^2

%%
S = 225.7235;                                     % area of the wing in ft2
b = 42.4944;                                     % total span of the wing in ft
a = b/2 ;                                       % span of one wing in ft                          
c = 1.6190*3.28084 ;
AR = 8 ;                                        % aspect ratio of the wing
MAC = 5.3117 ;                                   % Mean Aerodynamic chord of the wing in ft
%%
n1 = 3.1 ;                                   % maximum load factor
W = 3968.3;                                       % gross take off weight in lbs

Lift = n1*W*g ;                                   % total lift 

L = Lift/2 ;                                    % lift loads bore by one wing
%%
syms x_sw  % x_spanwise 
x1 = linspace(0,a,100) ; % x position varing from 0 --> b/2
ellip_spw_cd_vari = (4*S/pi*b)*sqrt(1 - (x_sw/a)^2);      %  elliptic chord distribution
lift_ellip_spw_vari = (4*Lift/(pi*b))*sqrt(1 - (x1/a).^2) ;
%L = L0*linspace(0,a,100);                       % elliptic lift distribution

trap_spw_cd_vari = (1 - ((1/a).*x1)*(1-lambda)) ;
total_trap_ar = trapz(x1,trap_spw_cd_vari) ;
L_trap_at_x0 = L/total_trap_ar ;
lift_trap_spw_vari = L_trap_at_x0.*(1 - ((1/a).*x1)*(1-lambda)) ;

lift_shrenk_spw_vari = (lift_trap_spw_vari + lift_ellip_spw_vari)/2 ;

width_val = 1.5 ;
figure(1)
plot (x1,lift_shrenk_spw_vari, 'k','linewidth',width_val) ;
title('Spanwise lift distribution using Schrenks Approximation');
xlabel('x(ft)');
ylabel('Lift per unit length (lb/ft)');
grid on
ax = gca;
ax.FontSize = 14;
ax.XAxis.LineWidth = 1.2;
ax.YAxis.LineWidth = 1.2;

%% Drag Distribution
Cd0_wing =0.00636;  % 0.009189 ;
parasite_drag = Cd0_wing ;                   % total parasite drag

KK = 0.0568 ;    %check
CL_max = 1.337 ; %check
induced_drag = KK*(CL_max)^2 ;

CD_total =  parasite_drag + induced_drag ;

Drag = (0.5)*rho_h*V_cruise^2*S*CD_total ;                  % in lbs
D0_pt_unit_span = Drag/b ;

x_d_total = linspace(0,a,100) ;
%x_d1 = (0:step_s:0.8*a) ; %linspace(0,.8*a,80) ;
D_First_prt = 0.95*(Drag/2)/(0.8*a) ;

%x_d2 = (0.8*a+step_s:step_s:a) ; % :linspace(x_d1(2)-x_d1(1)+0.8*a,a,20) ;
D_second_prt = 0.05*(Drag/2)/(0.2*a) ;

%x_d_total = [x_d1, x_d2] ;

D_dist_spw = [D_First_prt*linspace(1,1,80), D_second_prt*linspace(1,1,20) ]; 

figure(2)
plot(x_d_total,D_dist_spw,'k','linewidth',width_val);
ax = gca;
ax.FontSize = 14;
ax.XAxis.LineWidth = 1.2;
ax.YAxis.LineWidth = 1.2;
title('Spanwise drag distribution using Schrenks Approximation');
xlabel('x(ft)');
ylabel('Drag per unit length (lb/ft)');
grid on

%% Bending Moment Dist
for i = 1 : length(x_d_total)
    x_end_st(i) = x_d_total(length(x_d_total)-i+1) ;
    F_d_per_len_end_st(i) = D_dist_spw(length(x_d_total)-i+1) ;
    F_l_per_len_end_st(i) = lift_shrenk_spw_vari(length(x_d_total)-i+1) ;
end
V_z_end_st = cumtrapz(x_end_st,F_d_per_len_end_st) ; 
M_y_end_st = cumtrapz(x_end_st,V_z_end_st) ; 

for i = 1 : length(x_d_total)
    M_y_due_to_D(i) = M_y_end_st(length(x_d_total)-i+1) ;
end
figure(3)
plot(x_d_total,M_y_due_to_D,'b','linewidth',width_val) ; 
ax = gca;
ax.FontSize = 14;
ax.XAxis.LineWidth = 1.2;
ax.YAxis.LineWidth = 1.2;
title('Spanwise M_y distribution');
xlabel('x(ft)');
ylabel('M_y (lb.ft)');
grid on

%% bending due to lift
V_y_end_st = cumtrapz(x_end_st,F_l_per_len_end_st) ;
M_z_end_st = cumtrapz(x_end_st,V_y_end_st) ; 
for i = 1 : length(x_d_total)
    M_z_due_to_L(i) = M_z_end_st(length(x_d_total)-i+1) ;
end

figure(4)
plot(x_d_total,M_z_due_to_L,'b','linewidth',width_val) ; 
ax = gca;
ax.FontSize = 14;
ax.XAxis.LineWidth = 1.2;
ax.YAxis.LineWidth = 1.2;
title('Spanwise M_z distribution');
xlabel('x(ft)');
ylabel('M_z (lb.ft)');
grid on
%% Twisting Moment dist
%% testing of lift and drag dist
vlaue_lift =  trapz(x1,lift_shrenk_spw_vari) ;
value_drag = trapz(x_d_total,D_dist_spw) ;
%% Save distributions
save('loads.mat','lift_shrenk_spw_vari','D_dist_spw','M_y_due_to_D','M_z_due_to_L')