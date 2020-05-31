clc,clear all
%% spar buckling -->
spar1_ytop = 0.1517*3.28084;
spar1_ybot = -0.08855*3.28084;
ycent      = 0.5*(spar1_ytop + spar1_ybot);

b1         = 0.2*(spar1_ytop - spar1_ybot);
half_span = 42.4944/2 ;
%% Thicknesses and Modulus
tw = 0.01*3.28084;      % ft
tf = tw/2;              % ft
tsc = 0.02*3.28084;     % ft
spr_cp_ar    = b1*tsc ; % ft^2
spr_flng_ar  = b1*tf  ; % ft^2

hw1 = (spar1_ytop-spar1_ybot-2*(tf+tsc));

% flexural rigidity of the composite beam E_st*I_sprCap + E_al*I_flange 
E_al = 1.5351*10^9*32.1522; % lbf/ft^2 
E_st = 4.2*10^9*32.1522;    % lbf/ft^2 

I_zz_sprCap = ( ((b1*tsc^3/12))+(spr_cp_ar*((tsc+2*tf+hw1)/2)^2) ) ;
I_zz_flange = ( ((b1*tf^3/12))+(spr_flng_ar*((tf/2)+(hw1/2))^2) ) + (tw*hw1^3/12)   ;

I_Sigmazz_SC = 2*( ((b1*tsc^3/12))+(spr_cp_ar*((tsc+2*tf+hw1)/2)^2) ) ;
I_Sigmazz_FL = 2*( ((b1*tf^3/12))+(spr_flng_ar*((tf/2)+(hw1/2))^2) )    ;
I_Sigmazz_web  = (tw*hw1^3/12) ;

I_buc_sprCap = ((b1*tsc^3/12)) ;
I_buc_flange = ((b1*tf^3/12))    ;
len_panel   = half_span/15 ;

Mz = (.8/2)*1941400 ; % lbf*ft

boom_area   = spr_cp_ar + spr_flng_ar ;
sgma_xx_max = 1.85e7 ;%Mz*E_st*(spar1_ytop-ycent)/( E_st*I_Sigmazz_SC + E_al*(I_Sigmazz_FL+I_Sigmazz_web) ) ;
%% for panel at root
F_comp = sgma_xx_max*boom_area 
P_crit = ( pi^2*(E_al*I_buc_sprCap + E_st*I_buc_flange) )/len_panel^2  

%% for spar at root
Fspar_comp = sgma_xx_max*boom_area 
Pspar_crit = ( pi^2*(E_al*I_zz_flange + E_st*I_zz_sprCap ))/(1.12*half_span)^2  

%% for plate(web)
V_rootSh = 10580 ; % lbf shear force at root

Kc = 3.6 ;
Ks = 5.3 ;
%
I_zz_flange = I_zz_flange + (hw1*tw^3/12);

F_xxC    =  1.85e7 %Mz*E_al*(spar1_ytop-tsc)/( E_st*I_zz_sprCap + E_al*I_zz_flange ) 
F_ul_Ccric =  Kc*E_al*(tw/hw1)^2    % critical compression lbf/ft2

%

h1 = hw1/2 ;
h2 = hw1/2 + tf + tsc ;
%
F_xySh = (V_rootSh/2)*(E_al*( (h1+tf)^2 - h1^2 ) + E_st*( h2^2 - (h1+tf)^2 ))/(E_st*I_zz_sprCap + E_al*I_zz_flange)
F_ul_Scric =  Ks*E_al*(tw/hw1)^2    % critical shear lbf/ft2
