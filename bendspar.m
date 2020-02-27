c=1.6190
spar1_y1=0.1517*3.28084;
spar1_y2=-0.08855*3.28084;
spar1_z=-0.3*c*3.28084;
spar2_y1=0.09546*3.28084;
spar2_y2=-0.031975*3.28084;
spar2_z=-0.7*c*3.28084;
b1=(spar1_y1-spar1_y2);
b2=(spar2_y1-spar2_y2);
%% thicknesses and Modulus
tw=0.004*3.28084;
tf=tw/2;
tsc=0.002*3.28084;
hw1=(spar1_y1-spar1_y2-tf-tsc);
hw2=(spar2_y1-spar2_y2-tf-tsc);
E_al= 69*10^9*0.000145038;
E_steel= 200*10^9*0.000145038;
E_star=E_steel/E_al;
ybar=[(spar1_y1+spar1_y2)/2,(spar2_y1+spar2_y2)/2];
zbar=[spar1_z,spar2_z];
%% Areas of individual parts
A1=[b1*tsc,b1*tsc,b1*tf,b1*tf,tw*hw1];
A2=[b2*tsc,b2*tsc,b2*tf,b2*tf,tw*hw2];
E=[E_steel,E_steel,E_al,E_al,E_al];
%% finding centroid
sigmaEA=[dot(A1,E),dot(A2,E)];
y_cent= dot(sigmaEA,ybar)/sum(sigmaEA);
z_cent= dot(sigmaEA,zbar)/sum(sigmaEA);
%% Finding Moment of Inertia
Izz1=2*((E_star*(b1*tsc^3/12))+(A1(1)*E_star*((tsc+2*tf+hw1)/2)^2))+ 2*(((b1*tf^3/12))+(A1(3)*((tw/4)+(hw1/2))^2))+ (tw*hw1^3/12);
Iyy1=2*E_star*(tsc*b1^3/12)+2*(tf*b1^3/12)+ (hw1*tw^3/12);
Izz2=2*((E_star*(b2*tsc^3/12))+(A2(1)*E_star*((tsc+2*tf+hw1)/2)^2))+ 2*(((b2*tf^3/12))+(A2(3)*((tw/4)+(hw2/2))^2))+ (tw*hw2^3/12);
Iyy2=2*E_star*(tsc*b2^3/12)+2*(tf*b2^3/12)+ (hw2*tw^3/12);

Iyy_total=Iyy1+ (sigmaEA(1)/E_al)*(zbar(1)-z_cent)^2 + Iyy2+ (sigmaEA(2)/E_al)*(zbar(2)-z_cent)^2;
Izz_total=Izz1+ (sigmaEA(1)/E_al)*(ybar(1)-y_cent)^2 + Izz2+ (sigmaEA(2)/E_al)*(ybar(2)-y_cent)^2;    
Iyz_total= (sigmaEA(1)/E_al)*(ybar(1)-y_cent)*(zbar(1)-z_cent) + (sigmaEA(2)/E_al)*(zbar(2)-z_cent)*(ybar(2)-y_cent);
%% finding V_xx, W_xx
A=[-E_star*Iyz_total,-E_star*Iyy_total ; E_star*Izz_total,E_star*Iyz_total];
M=[M_y_due_to_D(1),M_z_due_to_L(1)];
v=A\M';
%% finding sigma_xx
epsilon_xx= (y*v(1) +x*v(2));  %conversion from meter to feet
sigma_xx= E_al*epsilon_xx;
plot(sigma_xx)