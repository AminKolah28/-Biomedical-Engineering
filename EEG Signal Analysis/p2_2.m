close all;
clear;
clc;
d0=load('matlab.mat');
Data0=d0.stage0.Data(:,2);
Fs=d0.stage0.Fs(2);
N=length(Data0);
frequencies0=Fs*(-N/2:N/2-1)/N;
t0 = linspace(0, 10, length(Data0)); 
figure;
plot(t0, Data0);
xlabel('t');
ylabel('Amplitude');
title('stage 0');


Data1=d0.stage1.Data(:,2);
Fs=d0.stage1.Fs(2);
N=length(Data1);
frequencies1=Fs*(-N/2:N/2-1)/N;
t1 = linspace(0, 10, length(Data1)); 
figure;
plot(t1, Data1);
xlabel('t');
ylabel('Amplitude');
title('stage 1');

Data2=d0.stage2.Data(:,2);
Fs=d0.stage2.Fs(2);
N=length(Data2);
frequencies2=Fs*(-N/2:N/2-1)/N;
t2 = linspace(0, 10, length(Data2)); 
figure;
plot(t2, Data2);
xlabel('t');
ylabel('Amplitude');
title('stage 2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_f0=fftshift(fft(Data0));

%gamma
fmin_gamma = 35;
fmax_gamma = 800;

gamma_keep = find(frequencies0 >= fmin_gamma & frequencies0 <= fmax_gamma);

Gamma=data_f0(gamma_keep);
gamma=ifft(ifftshift(Gamma))
t_gamma=(0:length(Gamma)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies0(gamma_keep), abs(Gamma));
title('gamma 0');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_gamma, real(gamma));
title('gamma 0');
xlabel('t');
ylabel('Magnitude');

power_gamma_0 = trapz(t_gamma, abs(gamma).^2) / (t_gamma(end) - t_gamma(1));
disp(['power_gamma: ', num2str(power_gamma_0)]);

%alpha
fmin_alpha = 8;
fmax_alpha = 13;

alpha_keep = find(frequencies0 >= fmin_alpha & frequencies0 <= fmax_alpha);

Alpha=data_f0(alpha_keep);
alpha=ifft(ifftshift(Alpha))
t_alpha=(0:length(Alpha)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies0(alpha_keep), abs(Alpha));
title('alpha 0');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_alpha, real(alpha));
title('alpha 0');
xlabel('t');
ylabel('Magnitude');

power_alpha_0 = trapz(t_alpha, abs(alpha).^2) / (t_alpha(end) - t_alpha(1));
disp(['power_alpha: ', num2str(power_alpha_0)]);

%beta
fmin_beta = 13;
fmax_beta = 35;

beta_keep = find(frequencies0 >= fmin_beta & frequencies0 <= fmax_beta);

Beta=data_f0(beta_keep);
beta=ifft(ifftshift(Beta))
t_beta=(0:length(Beta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies0(beta_keep), abs(Beta));
title('beta 0');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_beta, real(beta));
title('beta 0');
xlabel('t');
ylabel('Magnitude');

power_beta_0 = trapz(t_beta, abs(beta).^2) / (t_beta(end) - t_beta(1));
disp(['power_beta: ', num2str(power_beta_0)]);

%theta
fmin_theta = 4;
fmax_theta = 8;

theta_keep = find(frequencies0 >= fmin_theta & frequencies0 <= fmax_theta);

Theta=data_f0(theta_keep);
theta=ifft(ifftshift(Theta))
t_theta=(0:length(Theta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies0(theta_keep), abs(Theta));
title('theta 0');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_theta, real(theta));
title('theta 0');
xlabel('t');
ylabel('Magnitude');

power_theta_0 = trapz(t_theta, abs(theta).^2) / (t_theta(end) - t_theta(1));
disp(['power_theta: ', num2str(power_theta_0)]);

%delta
fmin_delta = 0.5;
fmax_delta = 4;

delta_keep = find(frequencies0 >= fmin_delta & frequencies0 <= fmax_delta);

Delta=data_f0(delta_keep);
delta=ifft(ifftshift(Delta))
t_delta=(0:length(Delta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies0(delta_keep), abs(Delta));
title('delta 0');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_delta, real(delta));
title('delta 0');
xlabel('t');
ylabel('Magnitude');

power_delta_0 = trapz(t_delta, abs(delta).^2) / (t_delta(end) - t_delta(1));
disp(['power_delta: ', num2str(power_delta_0)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_f1=fftshift(fft(Data1));

%gamma
fmin_gamma = 35;
fmax_gamma = 800;

gamma_keep = find(frequencies1 >= fmin_gamma & frequencies1 <= fmax_gamma);

Gamma=data_f1(gamma_keep);
gamma=ifft(ifftshift(Gamma))
t_gamma=(0:length(Gamma)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies1(gamma_keep), abs(Gamma));
title('gamma 1');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_gamma, real(gamma));
title('gamma 1');
xlabel('t');
ylabel('Magnitude');

power_gamma_1 = trapz(t_gamma, abs(gamma).^2) / (t_gamma(end) - t_gamma(1));
disp(['power_gamma: ', num2str(power_gamma_1)]);

%alpha
fmin_alpha = 8;
fmax_alpha = 13;

alpha_keep = find(frequencies1 >= fmin_alpha & frequencies1 <= fmax_alpha);

Alpha=data_f1(alpha_keep);
alpha=ifft(ifftshift(Alpha))
t_alpha=(0:length(Alpha)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies1(alpha_keep), abs(Alpha));
title('alpha 1');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_alpha, real(alpha));
title('alpha 1');
xlabel('t');
ylabel('Magnitude');

power_alpha_1 = trapz(t_alpha, abs(alpha).^2) / (t_alpha(end) - t_alpha(1));
disp(['power_alpha: ', num2str(power_alpha_1)]);

%beta
fmin_beta = 13;
fmax_beta = 35;

beta_keep = find(frequencies1 >= fmin_beta & frequencies1 <= fmax_beta);

Beta=data_f1(beta_keep);
beta=ifft(ifftshift(Beta))
t_beta=(0:length(Beta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies1(beta_keep), abs(Beta));
title('beta 1');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_beta, real(beta));
title('beta 1');
xlabel('t');
ylabel('Magnitude');

power_beta_1 = trapz(t_beta, abs(beta).^2) / (t_beta(end) - t_beta(1));
disp(['power_beta: ', num2str(power_beta_1)]);


%theta
fmin_theta = 4;
fmax_theta = 8;

theta_keep = find(frequencies1 >= fmin_theta & frequencies1 <= fmax_theta);

Theta=data_f1(theta_keep);
theta=ifft(ifftshift(Theta))
t_theta=(0:length(Theta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies1(theta_keep), abs(Theta));
title('theta 1');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_theta, real(theta));
title('theta 1');
xlabel('t');
ylabel('Magnitude');

power_theta_1 = trapz(t_theta, abs(theta).^2) / (t_theta(end) - t_theta(1));
disp(['power_theta: ', num2str(power_theta_1)]);

%delta
fmin_delta = 0.5;
fmax_delta = 4;

delta_keep = find(frequencies1 >= fmin_delta & frequencies1 <= fmax_delta);

Delta=data_f1(delta_keep);
delta=ifft(ifftshift(Delta))
t_delta=(0:length(Delta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies1(delta_keep), abs(Delta));
title('delta 1');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_delta, real(delta));
title('delta 1');
xlabel('t');
ylabel('Magnitude');

power_delta_1 = trapz(t_delta, abs(delta).^2) / (t_delta(end) - t_delta(1));
disp(['power_delta: ', num2str(power_delta_1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_f2=fftshift(fft(Data2));

%gamma
fmin_gamma = 35;
fmax_gamma = 800;

gamma_keep = find(frequencies2 >= fmin_gamma & frequencies2 <= fmax_gamma);

Gamma=data_f2(gamma_keep);
gamma=ifft(ifftshift(Gamma))
t_gamma=(0:length(Gamma)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies2(gamma_keep), abs(Gamma));
title('gamma 2');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_gamma, real(gamma));
title('gamma 2');
xlabel('t');
ylabel('Magnitude');

power_gamma_2 = trapz(t_gamma, abs(gamma).^2) / (t_gamma(end) - t_gamma(1));
disp(['power_gamma: ', num2str(power_gamma_2)]);

%alpha
fmin_alpha = 8;
fmax_alpha = 13;

alpha_keep = find(frequencies2 >= fmin_alpha & frequencies2 <= fmax_alpha);

Alpha=data_f2(alpha_keep);
alpha=ifft(ifftshift(Alpha))
t_alpha=(0:length(Alpha)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies2(alpha_keep), abs(Alpha));
title('alpha 2');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_alpha, real(alpha));
title('alpha 2');
xlabel('t');
ylabel('Magnitude');

power_alpha_2 = trapz(t_alpha, abs(alpha).^2) / (t_alpha(end) - t_alpha(1));
disp(['power_alpha: ', num2str(power_alpha_2)]);

%beta
fmin_beta = 13;
fmax_beta = 35;

beta_keep = find(frequencies2 >= fmin_beta & frequencies2 <= fmax_beta);

Beta=data_f2(beta_keep);
beta=ifft(ifftshift(Beta))
t_beta=(0:length(Beta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies2(beta_keep), abs(Beta));
title('beta 2');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_beta, real(beta));
title('beta 2');
xlabel('t');
ylabel('Magnitude');

power_beta_2 = trapz(t_beta, abs(beta).^2) / (t_beta(end) - t_beta(1));
disp(['power_beta: ', num2str(power_beta_2)]);

%theta
fmin_theta = 4;
fmax_theta = 8;

theta_keep = find(frequencies2 >= fmin_theta & frequencies2 <= fmax_theta);

Theta=data_f2(theta_keep);
theta=ifft(ifftshift(Theta))
t_theta=(0:length(Theta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies2(theta_keep), abs(Theta));
title('theta 2');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_theta, real(theta));
title('theta 2');
xlabel('t');
ylabel('Magnitude');

power_theta_2 = trapz(t_theta, abs(theta).^2) / (t_theta(end) - t_theta(1));
disp(['power_theta: ', num2str(power_theta_2)]);

%delta
fmin_delta = 0.5;
fmax_delta = 4;

delta_keep = find(frequencies2 >= fmin_delta & frequencies2 <= fmax_delta);

Delta=data_f2(delta_keep);
delta=ifft(ifftshift(Delta))
t_delta=(0:length(Delta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies2(delta_keep), abs(Delta));
title('delta 2');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_delta, real(delta));
title('delta 2');
xlabel('t');
ylabel('Magnitude');

power_delta_2 = trapz(t_delta, abs(delta).^2) / (t_delta(end) - t_delta(1));
disp(['power_delta: ', num2str(power_delta_2)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
disp(['power_gamma_0: ', num2str(power_gamma_0)]);
disp(['power_alpha_0: ', num2str(power_alpha_0)]);
disp(['power_beta_0: ', num2str(power_beta_0)]);
disp(['power_theta_0: ', num2str(power_theta_0)]);
disp(['power_delta_0: ', num2str(power_delta_0)]);

disp(['power_gamma_1: ', num2str(power_gamma_1)]);
disp(['power_alpha_1: ', num2str(power_alpha_1)]);
disp(['power_beta_1: ', num2str(power_beta_1)]);
disp(['power_theta_1: ', num2str(power_theta_1)]);
disp(['power_delta_1: ', num2str(power_delta_1)]);

disp(['power_gamma_2: ', num2str(power_gamma_2)]);
disp(['power_alpha_2: ', num2str(power_alpha_2)]);
disp(['power_beta_2: ', num2str(power_beta_2)]);
disp(['power_theta_2: ', num2str(power_theta_2)]);
disp(['power_delta_2: ', num2str(power_delta_2)]);