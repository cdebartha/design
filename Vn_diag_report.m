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
xlabel('\bf V (mph)') ; ylabel('\bf n'); title('v-n diagram');
legend('blue line');