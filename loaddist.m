v_cruise = 75*3.28084; % in ft/s
g = 9.8*3.28084; % in ft/s
w = 1800*0.0685218; % in slug
w_cruise = 0.97*0.985*w;
wbys = 82.0114*0.00636588; % in slug/ft^2

%%
pressure = 10.9174 % in lb/ft^2
density = 18.685e-4 % in slug/ft^3
viscosity = 3.58e-7 % in lb.s/ft^2
%%
S = 225.7235;                                     % area of the wing in ft2
b = 42.4944;                                     % total span of the wing in ft
a = b/2 ;                                       % span of one wing in ft                          

AR = 8 ;                                        % aspect ratio of the wing
MAC = 5.3117 ;                                   % Mean Aerodynamic chord of the wing in ft
%%
n1 = 3.1 ;                                   % maximum load factor
W = 3968.3;                                       % gross take off weight in lbs

Lift = n1*W ;                                   % total lift 

L = Lift/2 ;                                    % lift loads bore by one wing
%%
x = linspace(0,a,1000);
c1 = sqrt(1 - (x/a)^2);                         %  elliptic chord distribution
L0 = L/(trapz(c1, 0, a));
L = L0*linspace(0,a,100);                       % elliptic lift distribution



