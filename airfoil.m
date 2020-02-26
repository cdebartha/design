%% parameters
e = 0.7;
AR = 8;
v_cruise = 75;
cruise_rho = 0.963; % cruise altitude is 8000ft
g = 9.8;
w = 1800;
w_cruise = 0.97*0.985*w;
kapa = 1/(pi*e*AR) ;
%% data
a = load('xf-n63415-il-1000000.dat');
alpha = a(:,1);
cl = a(:,2);
cd = a(:,3);
cm = a(:,5);
clbycd = cl./cd;

%% finding clalpha2d
alpha_rad = alpha*pi/180;
p = polyfit(alpha_rad(37:103),cl(37:103),1);
clalpha2d = p(1);
cl_fit = polyval(p,alpha_rad);

%% finding clalpha3d
clalpha3d = clalpha2d/(1+clalpha2d/(pi*e*AR));
cl0 = cl(alpha_rad==0)*clalpha3d/clalpha2d;
coeff = [clalpha3d cl0];
cl3d = polyval(coeff,alpha_rad);

%% selecting cl design
alpha_design = 0.5;
cl_design = cl3d(alpha==alpha_design);

%% w/s calculation
wbys = 0.5*cruise_rho*v_cruise^2*cl_design/g;

%% s,b and c
s = w_cruise/wbys;
b = sqrt(s*AR);
c = s/b;

% plot
figure(1)
plot(alpha, cl3d);
hold on;
plot(alpha, cl);
hold on;
plot(alpha, cl_fit);
grid on;
ylabel('c_L');
xlabel('\alpha');
title('c_L vs \alpha');
figure(2)
plot(alpha,cd);
hold on;
plot(alpha_design,cd(alpha==alpha_design),'*');
grid on;
ylabel('c_D');
xlabel('\alpha');
title('c_D vs \alpha');
figure(3)
plot(cl,cd);
hold on;
plot(cl(alpha==alpha_design),cd(alpha==alpha_design),'*');
grid on;
ylabel('c_D');
xlabel('c_L');
title('c_L vs c_D');
figure(4)
plot(alpha,clbycd);
hold on;
plot(alpha_design,clbycd(alpha==alpha_design),'*');
grid on;
ylabel('c_L/c_D');
xlabel('\alpha');
title('c_L/c_D vs \alpha');

%% finding cm0
SM = 0.12;
slope = -SM;
cm0 = slope*(cl0-cl_design);

%% estimating lt, Vh, st, AR_t
r = .5 ;
lt = 0.49*b;
Vh = 0.7;
st = Vh*s*c/lt;
AR_t = 4.2;
b_t = sqrt(st*AR_t);
%c_t = st/b_t;
% tapered
c_t_root = (2*st)/(b_t*(1+r)) ;
c_t_mean = (2/3)*c_t_root*((1+r+r^2)/(1+r));
c_t_tip = c_t_root/2 ;
%% estimating clalpha3d_t
e_t = 1;
clalpha2d_t = 2*pi;
clalpha3d_t = clalpha2d_t/(1+clalpha2d_t/(pi*e_t*AR_t));

%% estimating cg, it
cgbyc = 0.5246;
cm0_w = cm(alpha==0)+cl0*(cgbyc-0.25);
iw = 0;
eps0 = 2*cl0/(pi*AR);
eta = 1;
it_rad = iw + eps0 - (cm0-cm0_w)/(eta*Vh*clalpha3d_t);
cm0_t = eta*Vh*clalpha3d_t*(iw + eps0 - it_rad);
it = it_rad*180/pi;

%% estimating np
eps_alpha = 2*clalpha3d/(pi*AR);
np = c*(0.25 + eta*Vh*clalpha3d_t/clalpha3d*(1-eps_alpha));
npbyc = np/c;
SM_est = npbyc - cgbyc;

%% calculating vertical tail dimensions
Vv = 0.04;
stv = Vv*s*b/lt;
% ht = stv/c_t;
% tapered
h_vt_tap = 2*stv/((1+r)*c_t_root) ;
Ar_vert = h_vt_tap^2/stv ;
%% calculating dimensions of control surfaces
% Elevator Sizing
sebyst = 0.4;
se = sebyst*st;
%c_e = se/b_t;
r_e = .6;
% c_e_root = (c_t_root/2.5) ;
% c_e_tip = c_e_root/2 ;
b_e = b_t ; %(2*se)/((1+r_e)*c_e_root);
c_e = se/b_e;
% Rudder Sizing
str = 0.5*stv;
%c_r = 0.5*c_t;
c_r_root = (2/(1+r))*str/h_vt_tap ;
c_r_tip = c_r_root/2 ;
% Aileron Sizing
c_a = 0.25*c;
b_a = 0.4*(b/2);
% Flap Sizing
c_f = 0.3*c;
b_f = 0.6*(b/2);

%% estimating cd0
lf = 0.7*b;
s_f = (0.150*lf)*(0.182*lf);
s_lg = 0.3*0.15;
l_strut = 0.2; 
s_strut = 2*0.07*2.915;
cd_pi = [0.0098, 1.2*0.063, 0.0043, 0.5, 0.74, 0.1*l_strut];
A_pi = [s, s_f, (st+stv), s_lg, 2*s_lg, s_strut];
f = cd_pi.*A_pi;
f_total = sum(f);
cd0 = 1.176*f_total/s;

%% estimating POWER required to climb @5m/s at 8000ft
%cd0 = 0.0311;
ROC = 5;
gamma = asin(ROC/v_cruise);
cl_climb = cl_design*cos(gamma);
cd = cd0 + 1/(pi*e*AR)*cl_climb^2;
P_r = cd*0.5*cruise_rho*v_cruise^3*s;
P_a = w_cruise*ROC + P_r;
P_e = 1.3*P_a; 

%% take-off distance and ground roll distance
% Engine Selected is of 300 hp
max_P = 300*745.7;
P_to = 0.7*max_P;

%% fuel tank positioning
w_f = 0.149*w;
rho_f = 720 ; % fuel used is Avgas
v_f = w_f/rho_f;
coordinates = load('n63415-il.dat');
x = .01*c*coordinates(:,1);
y = .01*c*coordinates(:,2);
b_ft = 0.4*(b);
A_ft = v_f/b_ft;
b1 = 0.115;
b2 = -0.04;
l2 = 0.9714 ;
l1 = l2 - A_ft/(b1-b2);
figure(5)
plot(x,y);
xlim([0 1]*c);
ylim([-0.5 0.5]*c);
hold on
plot([l1,l2,l2,l1,l1],[b1,b1,b2,b2,b1]);
grid on

%% derivatives
%% alpha
l_fuse = .7*b ;
cl_al_w = clalpha3d ;
% foe fuse
k2_k1 = .79 ;
So = 0.182*0.150*l_fuse^2 ;
Vb = 137.95 ;%% s1 = .112l* l_f * .136 lf  and s2 = So  h = .25lf -->   s3 = 0 
calpha_fuse = (2*(k2_k1)*So)/Vb^(2/3);
cl_al_t_cont = clalpha3d_t*(1- (2*clalpha3d)/(pi*AR))*(st/s)*1 ;
Clalpha_total = cl_al_w + calpha_fuse + cl_al_t_cont ;
Cdalpha_total = (2*cl_design*clalpha3d)*kapa ;
% Cmalpha
cmalpha_w = (1 + 2*cl_design*kapa*alpha_design*(pi/180) + cd/clalpha3d )*(cgbyc-0.25)*clalpha3d + ...
              (2*cl_design*kapa*(-alpha_design*(pi/180))*cl_design/clalpha3d)*0*clalpha3d ;
cmalpha_t = -(clalpha3d_t*(1 -(2*clalpha3d)/(pi*AR))*(st/s)*(lt/c)*1 ) ;
pos_wing_ac_in_fuse = (.35 - (cgbyc-0.25)*(c/lf)) ;
kf = .6 ;
wf = .15*lf ;
cmalpha_fuse = kf*wf^2*lf/(s*c) ;
Cmalpha_total = cmalpha_w +cmalpha_t +cmalpha_fuse ;
%% alpha dot
lt_dash = lt + (cgbyc-0.25)*c ;
Clalpha_dot = 2*clalpha3d_t*(2*clalpha3d)/(pi*AR)*(lt_dash/c)*(st/s)*1 ;
Cdalpha_dot = 0 ;
Cmalpha_dot = -Clalpha_dot*(lt/c) ;
%% q pitch rate
x_dash = -(cgbyc-0.25)*c ;
clq_w = 2*cl_al_w*(x_dash/c);
clq_t = 2*clalpha3d_t*(lt/c)*(st/s)*1 ;
clq_total = clq_w + clq_t ;
Cdq_total = 0 ;
cmq_t = -(lt/c)*clq_t ;
cmq_w = -(abs(x_dash)/c)*clq_w ;
Cmq_total = cmq_w + cmq_t ;
%% delta_e
dlpha_t_delta_t = .6 ;
cldelta_e = clalpha3d_t*.5*(st/s) ;
cddelta_e = 0 ;
cmdelta_e = -(lt/c)*cldelta_e ;
%% beta
Ae =1.55*(h_vt_tap^2/stv) ;
clalpha3d_v_emp =  .057 ;
AR_vert = h_vt_tap^2/stv ;
clalpha3d_vert1  = clalpha2d_t/(1+clalpha2d_t/(pi*e_t*AR_vert));
kk = 1;
cl_alpha_fuse = .0525 ;
kkf = 1.8;

cybeta_t = -kk*clalpha3d_vert1*(.724 +1.53*(stv/s) + .009*AR_vert)*(stv/s) ;
cybeta_fuse = -kkf*cl_alpha_fuse*(Vb^(2/3)/s) ;
Cybeta_total = cybeta_t + cybeta_fuse ;

clbeta_w = cl_design*(-1*(.21*1 +.29)/AR + 0.05) + 0;
zv = h_vt_tap - 1.054 ;
clbeta_vert = -clalpha3d_vert1*(stv/s)*(zv/b)*1 ;
Clbeta_total = clbeta_w + clbeta_vert - .0006 + .00016 ;

cnbeta_w = cl_design^2*kapa/4 ;
Ss = 8.34 ;
cnbeta_fuse = -(.96*.148/57.3)*(8.34/s)*(lt/b)*(.5)^(1/2+1/3) ;
cnbeta_vert  = clalpha3d_v_emp*(stv/s)*(lt/b)*1 ;
Cnbeta = cnbeta_vert + cnbeta_w + cnbeta_fuse + .0002 - 0.0005 ;
%Cnbeta1 = -() - cybeta_t*(lt/b)
%% p rolling velocity
cyp_w = cl_design/AR ;
cyp_t = 0.0;
Cyp_total = 0 ;

clp_w = (-0.51*((AR + 4)/((2*pi*AR)/cl0 + 4)) - (1/8)*cd0) ;
clp_t = -.33 ;
clp_vt = 2*(zv/b)^2*cybeta_t ;
Clp_total = clp_w + clp_t + clp_vt ;

cnp_by_cl = -.06 ;
cnp_w = cnp_by_cl*cl_design ;
cnp_vt = 57.3*clalpha3d_v_emp*(stv/s)*(2/b^2)*(lt*zv);
Cnp_total = cnp_w + cnp_vt ;

%% r
cyr_w = 0 ;
cyr_v = -2*(lt/b)*cybeta_t  ; 
Cyr_total = cyr_v ;

clr_w = cl_design/3 ;
clr_t = -2*(lt/b)*(zv/b)*cybeta_t ;
Clr_total = clr_w + clr_t ;

dcnr_cl2 = -0.020;
dcnr_cdo = -0.30;
cnr_w = cl_design^2*dcnr_cl2 + cd0*dcnr_cdo ;
cnr_t = 2*(stv/s)*(lt/b)*cybeta_t ;
Cnr_total = cnr_w + cnr_t ;
%% aileron
Cydelta_a = 0;
Cldelta_a = (2*clalpha3d*0.45/(s*b))*(c*.6*(b)^2)/8 ;
K = -.18 ;
Cndelta_a =2*K*cl_design*Cldelta_a ;
%% rudder
Cydelta_r = clalpha3d_v_emp*0.66*(stv/s) ;
Cldetla_r = Cydelta_r*(zv/b) ;
Cndelta_r = -Cydelta_r*(lt/b) ;

%% take off distance
v_stall = 35  ;
M=w ;
gr = 9.81;
mu = .02 ;
alpha_stall = 15*pi/180 ;
clmax  = clalpha3d*alpha_stall ;
rho = 1.225 ;
v_to = 1.2*v_stall ;
acc = 2 ;
takeoff_dist = v_to^2/(2*acc) ;
Th = M*acc + (mu*M*gr + (0.0278 + kapa*clmax^2 - mu*clmax)*.5*rho*(.7*v_to)^2*s ) ;
power = Th*(v_to)/746