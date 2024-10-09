close all;
clear;
clc;
%1st signal
d1=load('v1.mat');
data1 = d1.val ;
selectedLine1 = data1(20, :);
t1 = linspace(0, 10, length(selectedLine1)); 
%timeInSeconds = -0:1/length(selectedLine):10;
figure;
plot(t1, selectedLine1);
xlabel('t');
ylabel('Amplitude');
title('Signal from Line 20');


data1_f=fftshift(fft(selectedLine1));
phase_data1 = angle(data1_f);
magnitude_data1 = abs(data1_f); 

figure;
subplot(2,1,1);
plot(phase_data1);
title('phase');
xlabel('f');
ylabel('phase');

subplot(2,1,2);
plot(magnitude_data1);
title('magnitude ');
xlabel('f');
ylabel('magnitude');
%2ed signal
d2=load('v2.mat');
data2 = d2.val ;
selectedLine2 = data2(20, :);
t2 = linspace(0, 10, length(selectedLine2)); 

figure;
plot(t2, selectedLine2);
xlabel('t');
ylabel('Amplitude');
title('Signal from Line 20');


data2_f=fftshift(fft(selectedLine2));
phase_data2 = angle(data2_f);
magnitude_data2 = abs(data2_f); 

figure;
subplot(2,1,1);
plot(phase_data2);
title('phase');
xlabel('f');
ylabel('phase');

subplot(2,1,2);
plot(magnitude_data2);
title('magnitude ');
xlabel('f');
ylabel('magnitude');
%3th signal
d3=load('v3.mat');
data3 = d3.val ;
selectedLine3 = data3(20, :);
t3 = linspace(0, 10, length(selectedLine3)); 

figure;
plot(t3, selectedLine3);
xlabel('t');
ylabel('Amplitude');
title('Signal from Line 20');


data3_f=fftshift(fft(selectedLine3));
phase_data3 = angle(data3_f);
magnitude_data3 = abs(data3_f); 

figure;
subplot(2,1,1);
plot(phase_data3);
title('phase');
xlabel('f');
ylabel('phase');

subplot(2,1,2);
plot(magnitude_data3);
title('magnitude ');
xlabel('f');
ylabel('magnitude');