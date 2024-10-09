close all;
clear;
clc;
d1=load('v1.mat');
data1 = d1.val ;
selectedLine1 = data1(20, :);
Fs = 1599; 
N = length(selectedLine1);
frequencies = Fs*(-N/2:N/2-1)/N;
data_f=fftshift(fft(selectedLine1));
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
title('gamma');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_gamma, real(gamma));
title('gamma');
xlabel('t');
ylabel('Magnitude');

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
title('alpha');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_alpha, real(alpha));
title('alpha');
xlabel('t');
ylabel('Magnitude');

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
title('beta');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_beta, real(beta));
title('beta');
xlabel('t');
ylabel('Magnitude');

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
title('theta');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_theta, real(theta));
title('theta');
xlabel('t');
ylabel('Magnitude');

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
title('delta');
xlabel('f');
ylabel('Magnitude');

subplot(2,1,2);
plot(t_delta, real(delta));
title('delta');
xlabel('t');
ylabel('Magnitude');
