spar1_y1=0.1517
spar1_y2=-0.08855
spar1_z=0.3*c
spar2_y1=0.09546
spar2_y2=-0.031975
spar2_z=0.7*c
l1=(spar1_y1-spar1_y2)
l2=(spar2_y1-spar2_y2)
v=[0,0]
%% thicknesses and 
tw=0.004
tf=tw/2
tsc=0.002
hw1=(spar1_y1-spar1_y2-tf-tsc)
hw2=(spar2_y1-spar2_y2-tf-tsc)
E_al= 69*10^9
E_steel= 200*10^9
E_star=E_steel/E_al
ybar=[(spar1_y1+spar1_y2)/2,(spar2_y1+spar2_y2)/2]
zbar=[spar1_z,spar2_z]
%% Areas of individual parts
A1=[l1*tsc,l1*tsc,l1*tf,l1*tf,tw*hw1]
A2=[l2*tsc,l2*tsc,l2*tf,l2*tf,tw*hw2]
E=[E_steel,E_steel,E_al,E_al,E_al];
%% finding centroid
sigmaEA=[dot(A1,E),dot(A2,E)]
y_cent= dot(sigmaEA,ybar)/sum(sigmaEA);
z_cent= dot(sigmaEA,zbar)/sum(sigmaEA);
%% Finding Moment of Inertia
Izz1=2*((E_star*(l1*tsc^3/12))+(A1(1)*E_star*((tsc+2*tf+hw1)/2)^2))+ 2*(((l1*tf^3/12))+(A1(3)*((tw/4)+(hw1/2))^2))+ (tw*hw1^3/12)
Iyy1=2*E_star*(tsc*l1^3/12)+2*(tf*l1^3/12)+ (hw1*tw^3/12)
Izz2=2*((E_star*(l2*tsc^3/12))+(A2(1)*E_star*((tsc+2*tf+hw1)/2)^2))+ 2*(((l2*tf^3/12))+(A2(3)*((tw/4)+(hw2/2))^2))+ (tw*hw2^3/12)
Iyy2=2*E_star*(tsc*l2^3/12)+2*(tf*l2^3/12)+ (hw2*tw^3/12)

Iyy_total=Iyy1+ sigmaEA(1)*(zbar(1)-z_cent)^2 + Iyy2+ sigmaEA(2)*(zbar(2)-z_cent)^2;
Izz_total=Izz1+ sigmaEA(1)*(ybar(1)-y_cent)^2 + Izz2+ sigmaEA(2)*(ybar(2)-y_cent)^2    
Iyz_total= sigmaEA(1)*(ybar(1)-y_cent)*(zbar(2)-z_cent) + sigmaEA(2)*(zbar(2)-z_cent)*(ybar(2)-y_cent)
%% finding V_xx, W_xx
A=[-E_star*Iyz_total,-E_star*Iyy_total ; E_star*Izz_total,E_star*Iyz_total]
M=[M_y_due_to_D(1),M_z_due_to_L(1)]
Ainv=inv(A);
v=Ainv/M';