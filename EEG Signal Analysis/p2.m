close all;
clear;
clc;
d1=load('matlab.mat');
Data=d1.stage0.Data(:,2);
Fs=d1.stage0.Fs(2);
N=length(Data);
frequencies=Fs*(-N/2:N/2-1)/N;
t1 = linspace(0, 10, length(Data)); 
figure;
plot(t1, Data);
xlabel('t');
ylabel('Amplitude');
title('stage 0');

%%%%%%%%%%
data_f=fftshift(fft(Data));

%gamma
fmin_gamma = 35;
fmax_gamma = 800;

gamma_keep = find(frequencies >= fmin_gamma & frequencies <= fmax_gamma);

Gamma=data_f(gamma_keep);
gamma=ifft(ifftshift(Gamma))
t_gamma=(0:length(Gamma)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies(gamma_keep), abs(Gamma));
title('gamma 0');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_gamma, real(gamma));
title('gamma 0');
xlabel('t');
ylabel('Magnitude');

power_gamma = trapz(t_gamma, abs(gamma).^2) / (t_gamma(end) - t_gamma(1));
disp(['power_gamma: ', num2str(power_gamma)]);

%alpha
fmin_alpha = 8;
fmax_alpha = 13;

alpha_keep = find(frequencies >= fmin_alpha & frequencies <= fmax_alpha);

Alpha=data_f(alpha_keep);
alpha=ifft(ifftshift(Alpha))
t_alpha=(0:length(Alpha)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies(alpha_keep), abs(Alpha));
title('alpha 0');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_alpha, real(alpha));
title('alpha 0');
xlabel('t');
ylabel('Magnitude');

power_alpha = trapz(t_alpha, abs(alpha).^2) / (t_alpha(end) - t_alpha(1));
disp(['power_alpha: ', num2str(power_alpha)]);

%beta
fmin_beta = 13;
fmax_beta = 35;

beta_keep = find(frequencies >= fmin_beta & frequencies <= fmax_beta);

Beta=data_f(beta_keep);
beta=ifft(ifftshift(Beta))
t_beta=(0:length(Beta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies(beta_keep), abs(Beta));
title('beta 0');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_beta, real(beta));
title('beta 0');
xlabel('t');
ylabel('Magnitude');

power_beta = trapz(t_beta, abs(beta).^2) / (t_beta(end) - t_beta(1));
disp(['power_beta: ', num2str(power_beta)]);

%theta
fmin_theta = 4;
fmax_theta = 8;

theta_keep = find(frequencies >= fmin_theta & frequencies <= fmax_theta);

Theta=data_f(theta_keep);
theta=ifft(ifftshift(Theta))
t_theta=(0:length(Theta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies(theta_keep), abs(Theta));
title('theta 0');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_theta, real(theta));
title('theta 0');
xlabel('t');
ylabel('Magnitude');

power_theta = trapz(t_theta, abs(theta).^2) / (t_theta(end) - t_theta(1));
disp(['power_theta: ', num2str(power_theta)]);

%delta
fmin_delta = 0.5;
fmax_delta = 4;

delta_keep = find(frequencies >= fmin_delta & frequencies <= fmax_delta);

Delta=data_f(delta_keep);
delta=ifft(ifftshift(Delta))
t_delta=(0:length(Delta)-1)/Fs;

figure;
subplot(2,1,1);
plot(frequencies(delta_keep), abs(Delta));
title('delta 0');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_delta, real(delta));
title('delta 0');
xlabel('t');
ylabel('Magnitude');

power_delta = trapz(t_delta, abs(delta).^2) / (t_delta(end) - t_delta(1));
disp(['power_delta: ', num2str(power_delta)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
disp(['power_gamma: ', num2str(power_gamma)]);
disp(['power_alpha: ', num2str(power_alpha)]);
disp(['power_beta: ', num2str(power_beta)]);
disp(['power_theta: ', num2str(power_theta)]);
disp(['power_delta: ', num2str(power_delta)]);

