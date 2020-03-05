c = 1.619*3.28084;
coordinates = load('n63415-il.dat');
x = .01*c*coordinates(:,1);
y = .01*c*coordinates(:,2);
x1 = 0.3*c;
x2 = 0.7*c; 
tw = 0.004*3.28084;
tf = 0.012*3.28084;
y11=0.1517*3.28084;
y12=-0.08855*3.28084;
y21=0.09546*3.28084;
y22=-0.031975*3.28084;
t1 = 0.1*(y11-y12);
t2 = 0.1*(y21-y22);
plot(x,y,'Linewidth',1.1);
hold on
plot([x1-t1/2, x1+t1/2, x1+t1/2, x1+tw/2, ...
    x1+tw/2, x1+t1/2, x1+t1/2, x1-t1/2, x1-t1/2, x1-tw/2, ...
    x1-tw/2, x1-t1/2, x1-t1/2],[y11, y11, y11-tf, y11-tf, y12+tf, y12+tf, ...
    y12, y12, y12+tf, y12+tf, y11-tf, y11-tf, y11],'Linewidth', 1.3);
xlim([0 1]*c);
ylim([-0.5 0.5]*c);
hold on
plot([x2-t2/2, x2+t2/2, x2+t2/2, x2+tw/2, ...
    x2+tw/2, x2+t2/2, x2+t2/2, x2-t2/2, x2-t2/2, x2-tw/2, ...
    x2-tw/2, x2-t2/2, x2-t2/2],[y21, y21, y21-tf, y21-tf, y22+tf, y22+tf, ...
    y22, y22, y22+tf, y22+tf, y21-tf, y21-tf, y21], 'Linewidth', 1.3);
ax = gca;
ax.FontSize = 14;
ax.XAxis.LineWidth = 1.2;
ax.YAxis.LineWidth = 1.2;
hold on
z_stinger=[x(12),x(9),x(40),x(43)];
y_stinger=[y(12),y(9),y(40),y(43)];
% plot(z_stinger,y_stinger, '*k')