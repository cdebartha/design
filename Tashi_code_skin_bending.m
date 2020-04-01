clear all; clc;

wingdata = importdata('n63415-il.dat');
%c = -0.16604338070307*x + 2.22;
x = linspace(0,6.685,10); %% 10 uniform section alonf span

E_ref = 73.1e9;
t_skin = 1.9e-3;
E_skin = 73.1e9; Y_skin = 855e6;           % Skin Material 2024T3 Aluminium

E_spar = 70e9; Y_spar = 254972900;         % Spar Material 6082T6 Aluminium

t_flange = 2e-3; l_flange = 10e-3;          %b_flange = span/2; 
t_web = 4e-3;  

for i=1:9
    chord = -0.16604338070307*((x(i)+x(i+1))/2) + 2.22;
    z_cordi = wingdata(:,1)*chord/100;
    y_cordi = wingdata(:,2)*chord/100;
    
    for j=1:41
        s = sqrt((y_cordi(i)-y_cordi(i+1))^2 + ((z_cordi(i)-z_cordi(i+1))^2));
        dA(j) = s * t_skin * (E_skin/E_ref);
        ydA(j) = dA(j)*0.5*(y_cordi(i)+y_cordi(i+1));
        zdA(j) = dA(j)*0.5*(z_cordi(i)+z_cordi(i+1));
        y2dA(j) = ydA(j)*0.5*(y_cordi(i)+y_cordi(i+1));
        z2dA(j) = zdA(j)*0.5*(z_cordi(i)+z_cordi(i+1));
        yzdA(j) = zdA(j)*0.5*(y_cordi(i)+y_cordi(i+1));
    end
    
    A_skin = sum(dA);
    ybar_skin = sum(ydA)/A_skin;
    zbar_skin = sum(zdA)/A_skin;
    Iz_skin = sum(y2dA);
    Iy_skin = sum(z2dA);
    Iyz_skin = sum(yzdA);
    
    h_web_fs = 12.843*chord/100;
    h_web_rs = 2.992*chord/100;
    A_fs = (2*l_flange*t_flange + h_web_fs*t_web)*E_spar/E_ref;
    A_rs = (2*l_flange*t_flange + h_web_rs*t_web)*E_spar/E_ref;
    
    ybar_fs = 0.3*chord;
    zbar_fs = 3.158*chord/100;
    ybar_rs = 0.7*chord;
    zbar_rs = 1.861*chord/100;
    
    Iy_fs = (l_flange*t_flange^3)/6 + 2*(l_flange*t_flange)*(t_flange/2 + h_web_fs - zbar_fs)^2 + (t_web*h_web_fs^3)/12;
    Iz_fs = (t_flange*l_flange^3)/6 +(t_web^3*h_web_fs)/12;
    Iyz_fs = 0;
    
    Iy_rs = (l_flange*t_flange^3)/6 + 2*(l_flange*t_flange)*(t_flange/2 + h_web_rs - zbar_rs)^2 + (t_web*h_web_rs^3)/12;
    Iz_rs = (t_flange*l_flange^3)/6 +(t_web^3*h_web_rs)/12;
    Iyz_fs = 0;
    
    ybar_final = (A_skin*ybar_skin + A_fs*ybar_fs + A_rs*ybar_rs)/(A_skin+A_fs+A_rs);
    zbar_final = (A_skin*zbar_skin + A_fs*zbar_fs + A_rs*zbar_rs)/(A_skin+A_fs+A_rs);
    
    Iy_final = Iy_fs + A_fs*(zbar_fs-zbar_final)^2 + Iy_rs + A_rs*(zbar_rs-zbar_final)^2 + Iy_skin + A_skin*(zbar_skin-zbar_final)^2;
    Iz_final = Iz_fs + A_fs*(ybar_fs-ybar_final)^2 + Iz_rs + A_rs*(ybar_rs-ybar_final)^2 + Iz_skin + A_skin*(ybar_skin-ybar_final)^2;
    Iyz_final = A_fs*(ybar_fs-ybar_final)*(zbar_fs-zbar_final) + A_rs*(zbar_rs-zbar_final)*(ybar_rs-ybar_final) + Iyz_skin + A_skin*(zbar_skin-zbar_final)*(ybar_skin-ybar_final);    
end