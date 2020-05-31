clc,clear all
%% spar buckling -->
spar1_ytop = 0.1517*3.28084;
spar1_ybot = -0.08855*3.28084;
ycent      = 0.5*(spar1_ytop + spar1_ybot);

b1         = 0.1*(spar1_ytop - spar1_ybot);
half_span = 42.4944/2 ;
%% Thicknesses and Modulus
tw = 0.008*3.28084;      % ft
tf = tw/2;              % ft
tsc = 0.01*3.28084;     % ft
spr_cp_ar    = b1*tsc ; % ft^2
spr_flng_ar  = b1*tf  ; % ft^2

hw1 = (spar1_ytop-spar1_ybot-2*(tf+tsc));

% flexural rigidity of the composite beam E_st*I_sprCap + E_al*I_flange 
E_al = 1.5351*10^9*32.1522; % lbf/ft^2 
E_st = 4.2*10^9*32.1522;    % lbf/ft^2 

I_zz_sprCap = ( ((b1*tsc^3/12))+(spr_cp_ar*((tsc+2*tf+hw1)/2)^2) ) ;
I_zz_flange = ( ((b1*tf^3/12))+(spr_flng_ar*((tf/2)+(hw1/2))^2) )   ;


y_boom = ( E_al*spr_flng_ar*(tf/2) + E_st*spr_cp_ar*(tf+(tsc/2)) )/( (E_al*spr_flng_ar)+(E_st*spr_cp_ar) ) ;

I_buc_sprCap = ((b1*tsc^3/12)) + spr_cp_ar*( tf+(tsc/2) - y_boom )^2 ;
I_buc_flange = ((b1*tf^3/12)) + spr_flng_ar*( (tf/2) - y_boom )^2    ;
len_panel   = half_span/11 ;

Mz = (.8/2)*1941400 ; % lbf*ft

boom_area   = spr_cp_ar + spr_flng_ar ;
sgma_xx_max = 3.363e7 ;%Mz*E_st*(spar1_ytop-ycent)/( E_st*I_Sigmazz_SC + E_al*(I_Sigmazz_FL+I_Sigmazz_web) ) ;
%% for panel at root
F_comp = sgma_xx_max*boom_area 
P_crit = ( pi^2*(E_al*I_buc_sprCap + E_st*I_buc_flange) )/len_panel^2  

%% for spar at root
Fspar_comp = sgma_xx_max*boom_area 
Pspar_crit = ( pi^2*(E_al*I_zz_flange + E_st*I_zz_sprCap ))/(1.12*half_span)^2  

%% for plate(web)
V_rootSh = 10580 ; % lbf shear force at root

Kc = 3.6 ; % a/b = 2.77
Ks = 5.3 ; % 
%
I_zz_flange = I_zz_flange + (tw*hw1^3/12) ;

F_xxC    =  sgma_xx_max %Mz*E_al*(spar1_ytop-tsc)/( E_st*I_zz_sprCap + E_al*I_zz_flange ) 
F_ul_Ccric =  Kc*E_al*(tw/hw1)^2    % critical compression lbf/ft2

%

h1 = hw1/2 ;
h2 = hw1/2 + tf + tsc ;
%
F_xySh = (V_rootSh/2)*(E_al*( (h1+tf)^2 - h1^2 ) + E_st*( h2^2 - (h1+tf)^2 ))/(E_st*I_zz_sprCap + E_al*I_zz_flange)
F_ul_Scric =  Ks*E_al*(tw/hw1)^2    % critical shear lbf/ft2

%% -----------------XXXXXXXXXXXXXXXXXXXXX========== spar --2 ====<<<<<<<

% %% spar buckling -->
% spar2_ytop = 0.09546*3.28084;
% spar2_ybot = -0.031975*3.28084;
% ycent     = 0.5*(spar2_ytop + spar2_ybot);
% 
% b2         = 0.2*(spar2_ytop - spar2_ybot);
% half_span = 42.4944/2 ;
% %% Thicknesses and Modulus
% tw2 = 0.002*3.28084;      % ft
% tf2 = tw2/2;              % ft
% tsc2 = 0.001*3.28084;     % ft
% spr_cp_ar2    = b2*tsc2 ; % ft^2
% spr_flng_ar2  = b2*tf2  ; % ft^2
% 
% hw2 = (spar2_ytop-spar2_ybot-2*(tf2+tsc2));
% 
% % flexural rigidity of the composite beam E_st*I_sprCap + E_al*I_flange 
% E_al = 1.5351*10^9*32.1522; % lbf/ft^2 
% E_st = 4.2*10^9*32.1522;    % lbf/ft^2 
% 
% I_zz_sprCap2 = ( ((b2*tsc2^3/12))+(spr_cp_ar2*((tsc2+2*tf2+hw2)/2)^2) ) ;
% I_zz_flange2 = ( ((b2*tf2^3/12))+(spr_flng_ar2*((tf2/2)+(hw2/2))^2) ) + (tw2*hw2^3/12)   ;
% 
% 
% y_boom2 = ( E_al*spr_flng_ar2*(tf2/2) + E_st*spr_cp_ar2*(tf2+(tsc2/2)) )/( (E_al*spr_flng_ar2)+(E_st*spr_cp_ar2) ) ;
% 
% I_buc_sprCap2 = ((b2*tsc2^3/12)) + spr_cp_ar2*( tf2+(tsc2/2) - y_boom2 )^2 ;
% I_buc_flange2 = ((b2*tf2^3/12)) + spr_flng_ar2*( (tf2/2) - y_boom2 )^2    ;
% len_panel   = half_span/15 ;
% 
% Mz = (.8/2)*1941400 ; % lbf*ft
% 
% boom_area2   = spr_cp_ar2 + spr_flng_ar2 ;
% sgma_xx_max = 1.78e8 ;%Mz*E_st*(spar1_ytop-ycent)/( E_st*I_Sigmazz_SC + E_al*(I_Sigmazz_FL+I_Sigmazz_web) ) ;
% %% for panel at root
% F_comp2 = sgma_xx_max*boom_area2 
% P_crit2 = ( pi^2*(E_al*I_buc_sprCap2 + E_st*I_buc_flange2) )/len_panel^2  
% 
% %% for spar at root
% Fspar_comp2 = sgma_xx_max*boom_area2 
% Pspar_crit2 = ( pi^2*(E_al*I_zz_flange2 + E_st*I_zz_sprCap2 ))/(1.12*half_span)^2  
% 
% %% for plate(web)
% V_rootSh = 10580 ; % lbf shear force at root
% 
% Kc = 3.6 ;
% Ks = 5.3 ;
% %
% I_zz_flange2 = I_zz_flange2 + (hw2*tw2^3/12);
% 
% F_xxC2    =  sgma_xx_max ; %Mz*E_al*(spar1_ytop-tsc)/( E_st*I_zz_sprCap + E_al*I_zz_flange ) 
% F_ul_Ccric2 =  Kc*E_al*(tw2/hw2)^2    % critical compression lbf/ft2
% 
% %
% 
% h1 = hw2/2 ;
% h2 = hw2/2 + tf2 + tsc2 ;
% %
% F_xySh2 = (V_rootSh/2)*(E_al*( (h1+tf2)^2 - h1^2 ) + E_st*( h2^2 - (h1+tf2)^2 ))/(E_st*I_zz_sprCap2 + E_al*I_zz_flange2)
% F_ul_Scric2 =  Ks*E_al*(tw2/hw2)^2    % critical shear lbf/ft2
