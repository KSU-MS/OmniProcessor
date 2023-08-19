clear;
clc;
close all;
dirinfo = dir();
[filename,pathname] = uigetfile('*');
data = readtable(filename);
columns = {'Time (RTC)','Time( milli )','FL','FR','RL','RR','Steering','Gyro X','Gyro Y','Gyro Z','Accel X','Accel Y','Accel Z'};
%%
%%TimeRTC = table2array(data(:,1));
EpochTime = table2array(data(:,1));
FL = table2array(data(:,13));
FR = table2array(data(:,14));
RL = table2array(data(:,15));
RR = table2array(data(:,16));
Steering = table2array(data(:,3));
Gyro_x = table2array(data(:,10));
Gyro_y = table2array(data(:,11));
Gyro_z = table2array(data(:,12));
Accel_X = table2array(data(:,4));
Accel_Y = table2array(data(:,5));
Accel_Z = table2array(data(:,6));
ModifiedAccel_Z = Accel_Z-9.81;


dt_time= datetime(EpochTime, 'convertfrom','posixtime');
Time = timeofday(dt_time);

dt = 0.01;
n = length(Time);
Accel_X_fft = fft(Accel_X,n);
Accel_Y_fft = fft(Accel_Y,n);
Accel_Z_fft = fft(ModifiedAccel_Z,n);

FL_fft = fft(FL,n);
FR_fft = fft(FR,n);
RL_fft = fft(RL,n);
RR_fft = fft(RR,n);

PSD_Accel_X = Accel_X_fft.*conj(Accel_X_fft)/n;
PSD_Accel_Y = Accel_Y_fft.*conj(Accel_Y_fft)/n;
PSD_Accel_Z = Accel_Z_fft.*conj(Accel_Z_fft)/n;

PSD_FL = FL_fft.*conj(FL_fft)/n;
PSD_FR = FR_fft.*conj(FR_fft)/n;
PSD_RL = RL_fft.*conj(RL_fft)/n;
PSD_RR = RR_fft.*conj(RR_fft)/n;

freq = 1/(dt*n)*(0:n);
L = 1:floor(n/2);

indicies_accel_x = PSD_Accel_X > prctile(PSD_Accel_X,98.5);
PSD_Accel_X_CLEAN = PSD_Accel_X.*indicies_accel_x;
Accel_X_fft = indicies_accel_x.*Accel_X_fft;
ffilt_x = ifft(Accel_X_fft);

indicies_accel_y = PSD_Accel_Y > prctile(PSD_Accel_Y,98);
PSD_Accel_Y_CLEAN = PSD_Accel_Y.*indicies_accel_y;
Accel_Y_fft = indicies_accel_y.*Accel_Y_fft;
ffilt_y = ifft(Accel_Y_fft);

indicies_accel_z = PSD_Accel_Z > prctile(PSD_Accel_Z,98);
PSD_Accel_Z_CLEAN = PSD_Accel_Z.*indicies_accel_z;
Accel_Z_fft = indicies_accel_z.*Accel_Z_fft;
ffilt_z = ifft(Accel_Z_fft);

windowSize= 20;
b = (1/windowSize)*ones(1,windowSize);
a = 0.95;
Y_percent_error1 = abs((Accel_Y-ffilt_y)./((Accel_Y + ffilt_y)./2)).*100;
Y_percent_error2 = abs((Accel_Y-filter(b,a,Accel_Y)./((Accel_Y + filter(b,a,Accel_Y)./2)))).*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y ACCEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1);

subplot(3,1,1);
plot(Time/1000,Accel_Y,'b');
hold on
% plot(Time/1000,Accel_y_smooth,'k','LineWidth',3);
ylim([-3,3]);
grid on
title("Noisy Y Accelerometer Data")
xlabel("Time (s)")
ylabel("Acceleration (g)")
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

subplot(3,1,2);
plot(Time/1000,ffilt_y,'w');
grid on
ylim([-3 3]);
xlabel('Time (s)')
ylabel('Acceleration (g)')
title('Post-Fourier Transformed Data')
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

subplot(3,1,3);
plot(Time/1000,Accel_Y,'b');
hold on
plot(Time/1000,ffilt_y,'w','LineWidth',1.5);
grid on
title('Pre and Post Transformed Data')
xlabel('Time (s)')
ylabel('Acceleration (g)')
ylim([-3 3]);
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X ACCEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(2);
subplot(3,1,1);
plot(Time/1000,Accel_X,'r');
hold on
ylim([-3,3]);
grid on
title("Noisy X Accelerometer Data")
xlabel("Time (s)")
ylabel("Acceleration (g)")
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

% plot(freq(L),PSD_Accel_Y_CLEAN(L),'-','Color',[.5 .1 0],'LineWidth',2.5);

subplot(3,1,2);
plot(Time/1000,ffilt_x,'w');
grid on
ylim([-3 3]);
xlabel('Time (s)')
ylabel('Acceleration (g)')
title('Post-Fourier Transformed Data')
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

subplot(3,1,3);
plot(Time/1000,Accel_X,'r');
hold on
plot(Time/1000,ffilt_x,'w','LineWidth',2.5);
grid on
title('Pre and Post Transformed Data')
xlabel('Time (s)')
ylabel('Acceleration (g)')
ylim([-3 3]);
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z ACCEL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(3);
subplot(3,1,1);
plot(Time/1000,ModifiedAccel_Z,'g');
hold on
ylim([-3,3]);
grid on
title("Noisy Z Accelerometer Data")
xlabel("Time (s)")
ylabel("Acceleration (g)")
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

subplot(3,1,2);
plot(Time/1000,ffilt_z,'w');
grid on
ylim([-3 3]);
xlabel('Time (s)')
ylabel('Acceleration (g)')
title('Post-Fourier Transformed Data')
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

subplot(3,1,3);
plot(Time/1000,ModifiedAccel_Z,'g');
hold on
plot(Time/1000,ffilt_z,'w','LineWidth',1.25);
grid on
title('Pre and Post Transformed Data')
xlabel('Time (s)')
ylabel('Acceleration (g)')
ylim([-3 3]);
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
subplot(3,1,1);
plot(Time/1000,Accel_X,'r');
hold on
plot(Time/1000,ffilt_x,'w','LineWidth',1.25);
grid on
title('Pre and Post Transformed X Data')
xlabel('Time (s)')
ylabel('Acceleration (g)')
ylim([-3 3]);
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

subplot(3,1,2);
plot(Time/1000,Accel_Y,'b');
hold on
plot(Time/1000,ffilt_y,'w','LineWidth',1.25);
plot(Time/1000,filter(b,a,Accel_Y),'k','LineWidth',1.25)
grid on
title('Pre and Post Transformed Y Data')
xlabel('Time (s)')
ylabel('Acceleration (g)')
ylim([-3 3]);
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

subplot(3,1,3);
plot(Time/1000,ModifiedAccel_Z,'g');
hold on
plot(Time/1000,ffilt_z,'w','LineWidth',1.25);
grid on
title('Pre and Post Transformed Z Data')
xlabel('Time (s)')
ylabel('Acceleration (g)')
ylim([-3 3]);
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5);
subplot(3,1,1);
plot(freq(L),PSD_Accel_X(L),'r','Linewidth',1.25);
grid on
hold on
plot(freq(L),PSD_Accel_X_CLEAN(L),'w','Linewidth',1.25);
xlabel('Frequency (Hz) (X)')
ylabel ('Power/Frequency (dB/Hz) (X)')
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);


subplot(3,1,2);
plot(freq(L),PSD_Accel_Y(L),'b','Linewidth',1.25);
grid on
hold on
plot(freq(L),PSD_Accel_Y_CLEAN(L),'w','Linewidth',1.25);
xlabel('Frequency (Hz) (Y)')
ylabel ('Power/Frequency (dB/Hz) (Y)')
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

subplot(3,1,3);
plot(freq(L),PSD_Accel_Z(L),'g','Linewidth',1.25);
grid on
hold on
plot(freq(L),PSD_Accel_Z_CLEAN(L),'w','Linewidth',1.25);
xlabel('Frequency (Hz) (Z)')
ylabel ('Power/Frequency (dB/Hz) (Z)')
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indicies_FL = PSD_FL > 1;
PSD_FL_CLEAN = PSD_FL.*indicies_FL;
FL_fft = indicies_FL.*FL_fft;
ffilt_FL = ifft(FL_fft);

indicies_FR = PSD_FL > 2;
PSD_FR_CLEAN = PSD_FR.*indicies_FR;
FR_fft = indicies_FR.*FR_fft;
ffilt_FR = ifft(FR_fft);

indicies_RL = PSD_RL > 1;
PSD_RL_CLEAN = PSD_RL.*indicies_RL;
RL_fft = indicies_RL.*RL_fft;
ffilt_RL = ifft(RL_fft);

indicies_RR = PSD_RR > 1;
PSD_RR_CLEAN = PSD_RR.*indicies_RR;
RR_fft = indicies_RR.*RR_fft;
ffilt_RR = ifft(RR_fft);



figure(9);
plot(Time/1000, FL); %FL
hold on
plot(Time/1000, FR); %FR
plot(Time/1000, RL); %RL 
plot(Time/1000, RR); %RR
legend('Col 3','Col 4','Col 5','Col 6');
grid on

figure(12);
plot(Time/1000, Gyro_z,'g','LineWidth',2); % Gyro Z
hold on
grid on
title('Yaw Rate and Lateral Acceleration vs. Time')
xlabel('Time (s)')
ylabel('Yaw Rate (rad/s)','FontSize',12,'FontWeight','bold','Color','g');
yyaxis right
ylim([-3 3]);
ylabel('Acceleration (g)','FontSize',12,'FontWeight','bold','Color','b')
plot(Time/1000,ffilt_y,'b','LineWidth',1.25);
legend('Gyro Z','Y Accel');
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);

% Accel_Y_sgolay = sgolayfilt(Accel_Y,3,81);
% figure(13);
% plot(Time/1000,ffilt_y,'k');
% hold on
% plot(Time/1000,Accel_Y_sgolay,'b');
% 

% figure(13);
% plot(Time/1000,Accel_Y,'b');
% hold on
% plot(Time/1000,ffilt_y,'k','LineWidth',1.25);
% grid on
% plot(Time/1000,filter(b,a,Accel_Y),'r','LineWidth',1.25);
% 
% yyaxis right
% ylim([0,1000])
% plot(Time/1000,Y_percent_error2,'y')
% 
% figure(14);
% plot(Time/1000,Accel_Z,'g');
% hold on
% plot(Time/1000,ffilt_z,'k','LineWidth',1.25);
% grid on
% plot(Time/1000,filter(b,a,Accel_Z),'r','LineWidth',1.25);

figure(15);
plot(Time/1000,Steering,'w');
hold on
plot(Time/1000, FL,'r'); %FL
plot(Time/1000, FR,'b'); %FR
plot(Time/1000, RL,'g'); %RL 
plot(Time/1000, RR,'y'); %RR
% plot(Time/1000,Gyro_z,'c'); %GyroZ
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);
grid on
legend('\color{white}Steering','\color{white}FL','\color{white}FR','\color{white}RL','\color{white}RR')

figure(16);
scatter(ffilt_y,ffilt_x,10,'filled','w');
hold on
scatter(Accel_Y,Accel_X,10,'filled','r');
set(gca,'color',[.25 .25 .25])
set(gca,'GridColor',[1 1 1]);
grid on
xlim([-2 2])
ylim([-2 2])
title('GG Circle')
xlabel('Lateral Acceleration [G]')
ylabel('Longitudinal Acceleration [G]')