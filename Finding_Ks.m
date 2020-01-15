clc
clear
close all

data=load('Test2_zerodamp.txt');

alpha = [pi/2, -pi/2, -pi/2, pi/2, pi/2, -pi/2, 0];
a = [0, 0, 0, 0, 0, 0, 0];
d = [0.36, 0, 0.42, 0, 0.4, 0, 0.126];
theta=[0 0 0 0 0 0 -0.958709];
dh = [theta' d' a' alpha'];
h=0.001;
s=2000;

% BUILD ROBOT--------------------------------------------------------------
for i = 1:length(dh(:,1))
    L{i} = Link('d', dh(i,2), 'a', dh(i,3), 'alpha', dh(i,4));
end
rob = SerialLink([L{1} L{2} L{3} L{4} L{5} L{6} L{7}]);

qm= data(:,2:7)';

V= data(s+1:end,17);

A = data(s+1:end,18);

trial = data(s+1:end,19);

damp = data(s+1:end,20);

for i=1:length(qm)
    theta=[qm(:,i)' -0.958709]; %This should be in both ways
    [T,all] = rob.fkine(theta);
    phi_euler = atan2(T(2, 3), T(1, 3));
    theta_euler = atan2(sqrt(T(2, 3)^2 + T(1, 3)^2), T(3, 3));
    psi_euler = atan2(T(3, 2), -T(3, 1));
    xm1(:,i)=[T(1,4) T(2,4) T(3,4) phi_euler theta_euler psi_euler]';
    %[Ja,x(:,i)]=Jacobian(alpha,a,d,theta);
end
%xm1=data(:,3:8)';
xm=xm1(:,s+1:end);

d1 = designfilt('lowpassiir','FilterOrder',12,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',1000);

for i=1:size(xm,1)
    xm(i,:)=filtfilt(d1,xm(i,:));
end

ind_tn = find(trial == 1);

j=1;
n=1;
for i=1:length(ind_tn)
    if i < length(ind_tn)
        if ind_tn(i+1)-ind_tn(i) > 1
            clear new
            new=[xm(:,ind_tn(j):ind_tn(i)); V(ind_tn(j):ind_tn(i))'; A(ind_tn(j):ind_tn(i))'; damp(ind_tn(j):ind_tn(i))'];
            content{n}=new;
            n=n+1;
            j = i+1;
        end
    end
    if i==length(ind_tn)
        clear new
        new=[xm(:,ind_tn(j):ind_tn(i)); V(ind_tn(j):ind_tn(i))'; A(ind_tn(j):ind_tn(i))'; damp(ind_tn(j):ind_tn(i))'];
        content{n}=new;
    end
end

for i=2:size(content,2)-1
    length_vec(i-1)=size(content{i},2);
end
le=min(length_vec);

r=0.95;
n1=1;
n2=1;
%figure
for i=2:size(content,2)-1
    clear m X Vel Acc va va_max va_min
    m = content{i};
    X = m(1,:);
    direction=struct.dir;
    if dir== 1
    v=struct.xdot;
    a=struct.xdotdot;
    va = v*
%     Vel = diff(X)/h;
%     Acc = diff(Vel)/h;
    Vel=m(end-2,:);
    Acc=m(end-1,:);
    va = Vel(1:end).*Acc;
    %va=filtfilt(d1,va);
    
    va_max = max(va);
    va_min = min(va);
    
    kp(i-1,1) = -log((1-r)/(1+r))/va_max;
    kn(i-1,1) = -log((1+r)/(1-r))/va_min;
    
    if rem(i,2)==0
        X_plot1(n1,:)=X(1,1:le);
        va_plot1(n1,:)=va(1,1:le);
        %damp_plot1(n1,:)=m(end,1:le);
        damp_plot1=filtfilt(d1,m(end,1:le));
        n1=n1+1;
    else
         X_plot2(n2,:)=X(1,1:le);
         va_plot2(n2,:)=va(1,1:le);
         %damp_plot2(n2,:)=m(end,1:le);
         damp_plot2=filtfilt(d1,m(end,1:le));
         n2=n2+1;
    end
    %subplot(3,1,1)
%     hold on
%     plot(X,'color',[0.4 0.4 0.4]);
%     grid on
%     subplot(3,1,2)
%     hold on
%     plot(Vel,'color',[0.4 0.4 0.4]);
%     grid on
%     subplot(3,1,3)
%     hold on
%     plot(va,'color',[0.4 0.4 0.4]);
%     grid on

end

figure
subplot(3,1,1)
plot(mean(X_plot1,1),'k')
grid on
hold on
plot(mean(X_plot1,1)+std(X_plot1),'-.k')
hold on
plot(mean(X_plot1,1)-std(X_plot1),'-.k')
ylabel('$x (m)$','Interpreter','latex');
axis([0 3000 -0.1 0.02])
subplot(3,1,2)
plot(mean(va_plot1,1),'k')
grid on
hold on
plot(mean(va_plot1,1)+std(va_plot1),'-.k')
hold on
plot(mean(va_plot1,1)-std(va_plot1),'-.k')
ylabel('$\dot{x}\ddot{x}$','Interpreter','latex');
axis([0 3000 -0.7 0.6])
subplot(3,1,3)
plot(mean(damp_plot1,1),'k');
hold on
plot(mean(damp_plot1,1)+std(damp_plot1),'-.k')
hold on
plot(mean(damp_plot1,1)-std(damp_plot1),'-.k')
grid on
axis([0 3000 -30 20])
ylabel('$Damping (Ns/m)$','Interpreter','latex');
xlabel('Time (ms)');

figure
subplot(3,1,1)
plot(mean(X_plot2,1),'k')
grid on
hold on
plot(mean(X_plot2,1)+std(X_plot2),'-.k')
hold on
plot(mean(X_plot2,1)-std(X_plot2),'-.k')
ylabel('$x (m)$','Interpreter','latex');
axis([0 3000 -0.02 0.1])
subplot(3,1,2)
plot(mean(va_plot2,1),'k')
hold on
plot(mean(va_plot2,1)+std(va_plot2),'-.k')
hold on
plot(mean(va_plot2,1)-std(va_plot2),'-.k')
grid on
ylabel('$\dot{x}\ddot{x}$','Interpreter','latex');
axis([0 3000 -0.7 0.6])
subplot(3,1,3)
plot(mean(damp_plot2,1),'k');
hold on
plot(mean(damp_plot2,1)+std(damp_plot2),'-.k')
hold on
plot(mean(damp_plot2,1)-std(damp_plot2),'-.k')
grid on
axis([0 3000 -30 20])
ylabel('$Damping (Ns/m)$','Interpreter','latex');
xlabel('Time (ms)');

figure
subplot(3,1,1)
for j=1:size(X_plot1,1)
    hold on
    plot(X_plot1(j,:),'color',[0.7 0.7 0.7]);
end
hold on
plot(mean(X_plot1,1),'k')
grid on
ylabel('$x (m)$','Interpreter','latex');
axis([0 3000 -0.1 0.02])
subplot(3,1,2)
for j=1:size(va_plot1,1)
    hold on
    plot(va_plot1(j,:),'color',[0.7 0.7 0.7]);
end
hold on
plot(mean(va_plot1,1),'k')
grid on
ylabel('$\dot{x}\ddot{x}$','Interpreter','latex');
axis([0 3000 -0.7 0.6])
subplot(3,1,3)
for j=1:size(damp_plot1,1)
    hold on
    plot(damp_plot1(j,:),'color',[0.7 0.7 0.7]);
end
hold on
plot(mean(damp_plot1,1),'k');
grid on
axis([0 3000 -30 20])
ylabel('$Damping (Ns/m)$','Interpreter','latex');
xlabel('Time (ms)');

figure
subplot(3,1,1)
for j=1:size(X_plot2,1)
    hold on
    plot(X_plot2(j,:),'color',[0.7 0.7 0.7]);
end
hold on
plot(mean(X_plot2,1),'k')
grid on
ylabel('$x (m)$','Interpreter','latex');
axis([0 3000 -0.02 0.1])
subplot(3,1,2)
for j=1:size(va_plot2,1)
    hold on
    plot(va_plot2(j,:),'color',[0.7 0.7 0.7]);
end
hold on
plot(mean(va_plot2,1),'k')
grid on
ylabel('$\dot{x}\ddot{x}$','Interpreter','latex');
axis([0 3000 -0.7 0.7])
subplot(3,1,3)
for j=1:size(damp_plot2,1)
    hold on
    plot(damp_plot2(j,:),'color',[0.7 0.7 0.7]);
end
hold on
plot(mean(damp_plot2,1),'k');
grid on
axis([0 3000 -30 20])
ylabel('$Damping (Ns/m)$','Interpreter','latex');
xlabel('Time (ms)');

kp_ave = mean(kp);
kn_ave = mean(kn);



