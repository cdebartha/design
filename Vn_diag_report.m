clc , clear all
%%
rho_h = 0.963*(2.205/(3.28084)^3) ; % in lb/ft^3
rho_sl = 1.225*(2.205/(3.28084)^3) ;
cl_pos_max = 1.337 ; % per radian
cl_neg_max = 0.9329 ; % per radian
weight = 1800*(2.20462262)*9.81*(3.28084) ; % in lb
S_area = 20.9704*(3.28084^2)  ; % m^2 wing span area 
n_pos_max = 3.1 ;
n_neg_max = 1.25 ;
V_pos_stall = sqrt(2*n_pos_max*weight/(rho_h*S_area*cl_pos_max)) ; % ft/s
V_neg_stall = sqrt(2*n_neg_max*weight/(rho_h*S_area*cl_neg_max)) ; % ft/s
V_pos_stall_mph = V_pos_stall/1.466 ; %mph
V_neg_stall_mph = V_neg_stall/1.466 ; %mph
V_cruse = 75*(3.28084) ; %ft/s
V_cruse_mph = 75*(3.28084)/1.466 ; % mph

%% v-n diagram
vel_upto_clmax_mph = linspace(0,V_pos_stall_mph,35);
n_plus_upto_clmax = (1/weight)*(.5*1.466^2*rho_h.*(vel_upto_clmax_mph.^2)*S_area*cl_pos_max) ;
vel_aft_clmax_mph = linspace(V_pos_stall_mph+.01,1.4*V_cruse_mph,40) ;
n_plus_aft_clmax = n_plus_upto_clmax(end).*linspace(1,1,length(vel_aft_clmax_mph));

vel_neg_upto_clmax_mph = linspace(0,V_neg_stall_mph,35);
n_neg_upto_clmax = -(1/weight)*(.5*1.466^2*rho_h.*(vel_neg_upto_clmax_mph.^2)*S_area*cl_neg_max) ;
vel__neg_aft_clmax_mph = linspace(V_neg_stall_mph+.01,1.4*V_cruse_mph,40) ;
n_neg_aft_clmax = n_neg_upto_clmax(end).*linspace(1,1,length(vel__neg_aft_clmax_mph));
%% assembly
n = [n_plus_upto_clmax,n_plus_aft_clmax];
vel = [vel_upto_clmax_mph,vel_aft_clmax_mph] ;
n_neg = [n_neg_upto_clmax,n_neg_aft_clmax] ;
vel_neg = [vel_neg_upto_clmax_mph,vel__neg_aft_clmax_mph] ;
%% ploting
plot(vel,n,'k',vel_neg,n_neg,'b',[vel(end),vel_neg(end)],[n(end),n_neg(end)],'g','linewidth',1'); grid on;
xlabel('\bf V (mph)') ; ylabel('\bf n'); title('V-n diagram');
hold on
%% adding the corner point in V_n diagram for low mid and high

%%


%% low gust
wgh = 1800*(2.20462262)*9.81*(3.28084) ; % in lb ; % aircraft weight in lbs     
%S_area = 225.7235 ; % wing area in ft^2
%cl_pos_max = 1.337 ; % 
W_F_S = wgh/S_area ;

u_gust = [25, 56, 66] ; % in ft/s   low gust speed considered in upward dirn
V_cruise = 75*(3.28084)  ; % in ft/s cruze speed in right-hori
V_stall = sqrt(2*n_pos_max*wgh/(rho_h*S_area*1.54))  ; %in mph
V_plot = [1.4*V_cruise, V_cruise, 1.1*V_stall] ; 
cl_alpha3d = 4.7163 ; % per radian

Weight = 3968.321 ; % aircraft weight in lbs     
MAC = 5.3117 ; % in ft
W_m_S = Weight/S_area ; %v in lb/ft^2 wing loading

mu_g = (2*W_m_S)/(rho_sl*(0.031081)*MAC*cl_alpha3d) ;
K_g = 0.88*mu_g/(5.3 + mu_g) ;

dn_low = rho_h*(0.031081)*K_g*u_gust(1)*V_plot(1)*cl_alpha3d/(2*W_m_S) ;
dn_med = rho_h*(0.031081)*K_g*u_gust(2)*V_plot(2)*cl_alpha3d/(2*W_m_S) ;
dn_hig = rho_h*(0.031081)*K_g*u_gust(3)*V_plot(3)*cl_alpha3d/(2*W_m_S) ;

dn = [dn_low, dn_med, dn_hig] 
n_gust = [1+dn_low, 1+dn_med, 1+dn_hig] ;
n_neg_gust = [1-dn_low, 1-dn_med, 1-dn_hig] ;
V_plot_mph = (1/1.466).*V_plot ;
plot([0,V_plot_mph(1)],[1, n_gust(1)],'--k',[0,V_plot_mph(2)],[1,n_gust(2)],'--b',[0,V_plot_mph(3)],[1,n_gust(3)],'--r','linewidth',1) ;
grid on;
hold on;
plot([0,V_plot_mph(1)],[1,n_neg_gust(1)],'--k',[0,V_plot_mph(2)],[1,n_neg_gust(2)],'--b',[0,V_plot_mph(3)],[1,n_neg_gust(3)],'--r') ;
% legend('low gust','med gust','high gust');
ylim([-2.5,4.5]);
xlim([0,250]);
ax = gca;
ax.FontSize = 14;
ax.XAxis.LineWidth = 1.2;
ax.YAxis.LineWidth = 1.2;

