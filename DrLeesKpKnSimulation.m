clc; close all; clear all;
% Variable damping simulation (2018.05.23)

% Profiles are generated based on Minimum jerk trajectory
% xi (0 deg) to xf (5 deg) in t=d seconds, the minimum jerk trajectory is 
% x(t) = xi + (xf-xi)*(10(t/d)^3-15(t/d)^4+6(t/d)^5)

xi = 0;
% xf = 5/180*pi;  % rad
xf = 0.1; % meters
d = 0.3;        % second
Hz = 1000;      % measurement rate

t=[0:1/Hz:d]';
x= xi + (xf-xi)*(10*(t/d).^3-15*(t/d).^4+6*(t/d).^5);
v = diff(x)./diff(t);
tv = t(2:end);
a = diff(v)./diff(tv);
ta = tv(2:end);

va = v(2:end).*a;

va_max = max(va);
va_min = min(va);

%% State dependent damping (single function)
b_UB = 60;
b_LB = -20;
r = 0.8;       % At va_max, b_var will have the r*b_LB value!

k = -log((1-r)*b_LB/(r*b_LB-b_UB))/va_max;

b_var = -((b_UB-b_LB)./(1+exp(-k*va)))+b_UB;

%% State dependent damping (dual functions)
b_UB = 1;
b_LB = -0.5;
r = 0.95;       % At va_max, b_var will have the r*b_LB value!
b_var2 = zeros(length(va),1);

for i=1:length(va)
    if va(i,1)>=0
        kp = -log((1-r)/(1+r))/va_max;
        b_var2(i,1) = 2*b_LB/(1+exp(-kp*va(i,1)))-b_LB;
    else
        kn = -log((1+r)/(1-r))/va_min;
        b_var2(i,1) = -(2*b_UB/(1+exp(-kn*va(i,1)))-b_UB);
    end
end

fullscreen = get(0,'ScreenSize');
figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
set(gcf,'Color',[1,1,1]);

% subplot 521
subplot(5,1,1)
plot(t,x,'linewidth',2);
box off;
%title ('Displacement (rad)');
title('Position')
ylabel('$x \hspace{4pt} [m]$', 'Interpreter', 'latex', 'FontSize', 12)

% subplot 523
subplot(5,1,2)
plot(tv,v,'linewidth',2);
box off;
%title ('Velocity (rad/s)');
title('Velocity')
ylabel('$\dot{x} \hspace{4pt} [m/s]$', 'Interpreter', 'latex', 'FontSize', 12)

% subplot 525
subplot(5,1,3)
plot(ta,a,'linewidth',2);
box off;
%title ('Acceleration (rad/s^2)');
title('Acceleration') 
ylabel('$\ddot{x} \hspace{4pt} [m/s^2]$', 'Interpreter', 'latex', 'FontSize', 12)

% subplot 527
subplot(5,1,4)
plot(ta,va,'linewidth',2);
box off;
title('User Intent')
ylabel('$\dot{x}\ddot{x} \hspace{4pt} [m^2/s^3]$', 'Interpreter', 'latex', 'FontSize', 12)
%title ('Velocity * Acceleration (rad^2/s^3)');

% subplot 529
% subplot(5,1,4)
% plot(ta,b_var,'linewidth',2);
% text(max(ta)*0.1,max(b_var)*0.9,strcat('Movement time = ', num2str(d), 'sec.'), 'FontSize', 12);
% text(max(ta)*0.1,max(b_var)*0.7,strcat('Max v*a =', num2str(va_max), 'rad^2/s^3'), 'FontSize', 12);
% text(max(ta)*0.1,max(b_var)*0.5,strcat('k =', num2str(k), 'Nms/rad'), 'FontSize', 12);

% subplot (5,2,10)
subplot(5,1,5)
plot(ta,b_var2,'linewidth',2);
hold on
% plot(ta,b_var, 'linewidth',2);
% hold off
% legend('Dual Damping', 'Single Damping')
% legend('boxoff')
% legend('Location', 'northwest')

% text(max(ta)*0.1,max(b_var)*0.9,strcat('Movement time = ', num2str(d), 'sec.'), 'FontSize', 12);
% text(max(ta)*0.1,max(b_var)*0.7,strcat('Max v*a =', num2str(va_max), 'rad^2/s^3'), 'FontSize', 12);
% text(max(ta)*0.1,max(b_var)*0.5,strcat('Min v*a =', num2str(va_min), 'rad^2/s^3'), 'FontSize', 12);
% text(max(ta)*0.1,max(b_var)*0.3,strcat('kp =', num2str(kp), 'Nms/rad'), 'FontSize', 12);
% text(max(ta)*0.1,max(b_var)*0.1,strcat('kn =', num2str(kn), 'Nms/rad'), 'FontSize', 12);
title('Robotic Damping')
ylabel('$b_r \hspace{4pt} [Ns/m]$', 'Interpreter', 'latex', 'FontSize', 12)
xlabel('Time [s]', 'FontSize', 12)
box off;
%title ('State-dependent variable damping [Nms/rad]');

%% Simulation over a wide range of movement speeds (single function)
d1 = 0.3;       % second
d2 = 0.6;
d3 = 0.9;

t1=[0:1/Hz:d1]';
x1= xi + (xf-xi)*(10*(t1/d1).^3-15*(t1/d1).^4+6*(t1/d1).^5);
v1 = diff(x1)./diff(t1);
tv1 = t1(2:end);
a1 = diff(v1)./diff(tv1);
ta1 = tv1(2:end);
va1 = v1(2:end).*a1;
va1_max = max(va1);

k1 = -log((1-r)*b_LB/(r*b_LB-b_UB))/va1_max;
b1_var = -((b_UB-b_LB)./(1+exp(-k1*va1)))+b_UB;

t2=[0:1/Hz:d2]';
x2= xi + (xf-xi)*(10*(t2/d2).^3-15*(t2/d2).^4+6*(t2/d2).^5);
v2 = diff(x2)./diff(t2);
tv2 = t2(2:end);
a2 = diff(v2)./diff(tv2);
ta2 = tv2(2:end);
va2 = v2(2:end).*a2;
va2_max = max(va2);

k2 = -log((1-r)*b_LB/(r*b_LB-b_UB))/va2_max;
b2_var = -((b_UB-b_LB)./(1+exp(-k2*va2)))+b_UB;

t3=[0:1/Hz:d3]';
x3= xi + (xf-xi)*(10*(t3/d3).^3-15*(t3/d3).^4+6*(t3/d3).^5);
v3 = diff(x3)./diff(t3);
tv3 = t3(2:end);
a3 = diff(v3)./diff(tv3);
ta3 = tv3(2:end);
va3 = v3(2:end).*a3;
va3_max = max(va3);

k3 = -log((1-r)*b_LB/(r*b_LB-b_UB))/va3_max;
b3_var = -((b_UB-b_LB)./(1+exp(-k3*va3)))+b_UB;

fullscreen = get(0,'ScreenSize');
figure('Position',[0 0 fullscreen(3) fullscreen(4)]);
set(gcf,'Color',[1,1,1]);

subplot 321
plot(t1,x1,'linewidth',2,'color','r');
hold on;
plot(t2,x2,'linewidth',2,'color','g');
hold on;
plot(t3,x3,'linewidth',2,'color','b');
title ('Position profile (minimum jerk trajectory');
box off;

subplot 323
plot(ta1,b1_var,'linewidth',2,'color','r');
hold on;
plot(ta2,b2_var,'linewidth',2,'color','g');
hold on;
plot(ta3,b3_var,'linewidth',2,'color','b');
box off;
title ('Single function');

%% Simulation over a wide range of movement speeds (dual functions)va1_max = max(va1);
va1_min = min(va1);
va2_min = min(va2);
va3_min = min(va3);

b1_var2 = zeros(length(va1),1);

for i=1:length(va1)
    if va1(i,1)>=0
        kp = -log((1-r)/(1+r))/va1_max;
        b1_var2(i,1) = 2*b_LB/(1+exp(-kp*va1(i,1)))-b_LB;
    else
        kn = -log((1+r)/(1-r))/va1_min;
        b1_var2(i,1) = -(2*b_UB/(1+exp(-kn*va1(i,1)))-b_UB);
    end
end

b2_var2 = zeros(length(va2),1);

for i=1:length(va2)
    if va2(i,1)>=0
        kp = -log((1-r)/(1+r))/va2_max;
        b2_var2(i,1) = 2*b_LB/(1+exp(-kp*va2(i,1)))-b_LB;
    else
        kn = -log((1+r)/(1-r))/va2_min;
        b2_var2(i,1) = -(2*b_UB/(1+exp(-kn*va2(i,1)))-b_UB);
    end
end

b3_var2 = zeros(length(va3),1);

for i=1:length(va3)
    if va3(i,1)>=0
        kp = -log((1-r)/(1+r))/va3_max;
        b3_var2(i,1) = 2*b_LB/(1+exp(-kp*va3(i,1)))-b_LB;
    else
        kn = -log((1+r)/(1-r))/va3_min;
        b3_var2(i,1) = -(2*b_UB/(1+exp(-kn*va3(i,1)))-b_UB);
    end
end

subplot 325
plot(ta1,b1_var2,'linewidth',2,'color','r');
hold on;
plot(ta2,b2_var2,'linewidth',2,'color','g');
hold on;
plot(ta3,b3_var2,'linewidth',2,'color','b');
title ('Dual functions');
box off;