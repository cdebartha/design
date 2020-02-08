clear all;clc;
rho=0.964*0.062428;               %% density at 8000 ft FPS
AR=8;                             %% Aspect Ratio
e=0.7;
c=1.73*3.28084;                   %% Mean Aerodynamic Chord (feet)
S=22.33*10.7639;                  %% Surface Area of Wing   (feet^2)
W=1738.814*2.20462;               %% Weight of Aircraft     (pound)
clmax=1.54;                       %% Clmax 3D
clo3d=0.231;
clalpha2d=5.89;
clalpha3d=4.412;
va=sqrt((W*32.2*3.1)/(0.5*rho*clmax*S));  %% Max Stall Speed    (ft/s)
vb=1.1*va;                                %% Min Maneuvre Speed (ft/s)
vc=70*3.28084;                            %% Cruise Velocity    (ft/s) 
vd=1.4*vc;          
clmaxneg=clo3d+clalpha3d*(-14.25*pi/180);
v=linspace (0,322,323);
n_pos=zeros(323,1);
n_pos=0.5*rho.*v.*v*clmax*S/(W*32.2);
for i = 1 : 323
    if(n_pos(i)>3.1)
        n_pos(i)=3.1;
    end
end
n_neg=zeros(323,1);
n_neg=0.5*rho.*v.*v*clmaxneg*S/(W*32.2);
for i = 1 : 323
    if(n_neg(i)<=-1.25)
        n_neg(i)=-1.25;
    end
end
y=[-1.25,3.1];
x=[323,323];
vel=[vb vc vd].*0.592484;        %% in knots
v_gust=[66 50 25];               %% in FPS
mu_g=(2*((W/S)*32.2))/((rho*0.031081)*(c)*32.2*clalpha3d);    %%Unit Conversion
k_g=(0.88*mu_g)/(5.3+mu_g);
del_n=k_g.*v_gust.*vel*clalpha3d/(498*((W/S)*32.2));

ng_pos=1+del_n;
ng_neg=1-del_n;
n_cruise=[0,1];
n_gust=[0,1,1 ; vb,ng_pos(1),ng_neg(1) ; vc,ng_pos(2),ng_neg(2) ; vd,ng_pos(3),ng_neg(3) ]
plot(v,n_pos);
hold on;
plot(v,n_neg);
plot(x,y);
for i = 1 : 3
    y1=[vel(i)*1.68781,ng_pos(i)];
    y2=[vel(i)*1.68781,ng_neg(i)];
    plot([0,y1(1)],[1,y1(2)]);
    plot([0,y2(1)],[1,y2(2)]);

end
xlim([0 340]);

