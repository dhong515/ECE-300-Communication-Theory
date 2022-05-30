%% Danny Hong  ECE 300  HW 2
clc
clear
close all

load handel.mat %loading audio file
audiowrite('handel.wav', y, Fs) %writing audio file
[m, fs] = audioread('handel.wav'); %reading in the audio file

fs = 50*fs; %Upsampling sampling rate
m = m./max(abs(m(:))); %Normalizing the message sample
m = m(:,1).'; %Only one column is required
N = 50*(size(m, 2)); %Amount of times the message signal is being sampled
m = interp(m, 50); %Upsampling the message

%% Questions 1-6

%Question 1
t = linspace(1, N, N); %Declaring time values.

%Plotting m(t)
figure;
subplot(3, 1, 1); 
plot(t, m);
title("m(t) plot");
xlabel("time (s)");
ylabel("m(t)");

M = fftshift(fft(m)); %M(w) is the fourier transform of m(w)
wd = linspace(-pi, pi, N); %declaring angular frequency values
f = wd*fs/(2*pi); %declaring frequency values

%Plotting Amplitude plot for M(w)
subplot(3, 1, 2);
plot(f, abs(M));
title("Amplitude plot for M(\omega)");
xlabel("frequency (Hz)");
ylabel("|M(\omega)|");

%Plotting Phase plot for M(w)
subplot(3, 1, 3);
plot(f, unwrap(angle(M)));
title("Phase plot for M(\omega)");
xlabel("frequency (Hz)");
ylabel("\angleM(\omega)");

%The bandwidth is about 4.5 kHz as observed from the amplitude plot for M(w)

%Question 2 
%Variable definitions remain the same from question 1

Ac = 1; %Chose carrier amplitude to be 1 
fc = 550000; %Chose carrier frequency to be 550,000 Hz (AM radio frequency typically ranges from (550 kHz - 1700 kHz).
wc = 2*pi*fc; %Angular carrier frequency value (in radians)
c = Ac*cos(wc*t); %Carrier signal function c(t)

%Question 3
%Variable definitions remain the same from questions 1 and 2

u = m.*c; %u(t) is the DSB-SC modulated signal
U = fftshift(fft(u)); %U(w) is the Fourier Transform of DSBSC

%Plotting DSB-SC modulated signal
figure;
subplot(3, 1, 1);
plot(t, u);
title("DSB-SC modulated message plot");
xlabel("time (seconds)");
ylabel("u(t)");

%Plotting Amplitude plot of DSB-SC modulated signal
subplot(3, 1, 2);
plot(f, abs(U));
title("Fourier Amplitude plot of the DSB-SC modulated message");
xlabel("Frequency (Hz)");
ylabel("|U(\omega)|");

%Plotting Phase plot of DSB-SC modulated signal
subplot(3, 1, 3);
plot(f, unwrap(angle(U)));
title("Fourier Phase plot of the DSB-SC modulated message");
xlabel("Frequency (Hz)");
ylabel("\angleU(\omega)");

Ap = 1; %declaring pilot term coefficient to be 1
p = Ap*cos(wc*t); %pilot term p(t)

x = u + p; %x(t) is the DSB modulated signal which equals u(t) + p(t)
X = fftshift(fft(x)); %X(w) is the fourier transform of DSB

%Plotting DSB modulated signal
figure;
subplot(3, 1, 1);
plot(t, x);
title("DSB modulated message plot");
xlabel("time (seconds)");
ylabel("x(t)");

%Plotting Fourier Transform of DSB modulated signal
subplot(3, 1, 2);
plot(f, abs(X));
title("Fourier Amplitude plot of the DSB modulated message");
xlabel("Frequency (Hz)");
ylabel("|X(\omega)|");

%Plotting Fourier Transform of DSB modulated signal
subplot(3, 1, 3);
plot(f, unwrap(angle(X)));
title("Fourier Phase plot of the DSB modulated message");
xlabel("Frequency (Hz)");
ylabel("\angleX(\omega)");

%For DSB-SC modulation, the output ranges from -1 to 1 while for DSB modulation, the output ranges from 0 to 2.

%Question 4
%Variable definitions remain the same from questions 1-3

[mt_hilbert, mf_hilbert] = Hilbert(m, f); %Reusing Hilbert.m from HW1 to find the Hilbert Transform of m(t) and M(w)
ulssb = u + Ac.*mt_hilbert.*sin(wc*t); %ulssb(t) is the lower ssb AM signal
ulssbFT = fftshift(fft(ulssb)); %ulssbFT(w) is the Fourier Transform of ulssb(t)
uussb = u - Ac.*mt_hilbert.*sin(wc*t); %uussb(t) is the upper ssb AM signal
uussbFT = fftshift(fft(uussb));%uussbFT(w) is the Fourier Transform of uussb(t)

%Plotting ulssb(t) 
figure;
subplot(3, 1, 1);
plot(t, abs(ulssb));
title("Lower SSB AM modulated message plot");
xlabel("time (seconds)");
ylabel("ulssb(t)");

%Plotting Fourier Amplitude Plot of ulssb(t) 
subplot(3, 1, 2);
plot(f, abs(ulssbFT));
title("Fourier Amplitdue Plot of Lower SSB AM modulated message");
xlabel("Frequency (Hz)");
ylabel("|ulssbFT(\omega)|");

%Plotting Fourier Phase Plot of ulssb(t) 
subplot(3, 1, 3);
plot(f, unwrap(angle(ulssbFT)));
title("Fourier Phase Plot of Lower SSB AM modulated message");
xlabel("Frequency (Hz)");
ylabel("\angleulssbFT(\omega)");

%Plotting uussb(t)
figure;
subplot(3, 1, 1);
plot(f, abs(uussb));
title("Upper SSB AM modulated message plot");
xlabel("time (seconds)");
ylabel("uussb(t)");

%Plotting Fourier Amplitude Plot of uussb(t)
subplot(3, 1, 2);
plot(f, abs(uussbFT));
title("Fourier Amplitude Plot of Lower SSB AM modulated message");
xlabel("Frequency (Hz)");
ylabel("|uussbFT(\omega)|");

%Plotting Fourier Phase Plot of uussb(t)
subplot(3, 1, 3);
plot(f, unwrap(angle(ulssbFT)));
title("Fourier Phase Plot of Lower SSB AM modulated message");
xlabel("Frequency (Hz)");
ylabel("\angleuussbFT(\omega)");

%Question 5
%Variable definitions remain the same from questions 1-4

convAM = (1 + m).*c; %convAM(t) is the conventional AM signal
convAMFT = fftshift(fft(convAM)); %convAMFT(w) is the Fourier Transform of convAM(t)

%Plotting convAM(t)
figure;
subplot(3, 1, 1);
plot(t, convAM);
title("Conventional AM signal plot");
xlabel("time (seconds)");
ylabel("convAM(t)");

%Plotting Fourier Amplitude Plot of Conventional AM modulated signal
subplot(3, 1, 2);
plot(f, abs(convAMFT));
title("Fourier Amplitude Plot of Conventional AM modulated message");
xlabel("Frequency (Hz)");
ylabel("|convAMFT(\omega)|");

%Plotting Fourier Phase Plot of Conventional AM modulated signal
subplot(3, 1, 3);
plot(f, unwrap(angle(convAMFT)));
title("Fourier Phase Plot of Conventional AM modulated message");
xlabel("Frequency (Hz)");
ylabel("\angleconvAMFT(\omega)");

%Question 6
%Variable definitions remain the same from questions 1-5

convAM(convAM < 0) = 0; %Rectifying the Conventional AM signal
convAMFT = fftshift(fft(convAM)); %Fourier Transform of the rectified signal

%Plotting Rectified Conventional AM message signal
figure;
subplot(3, 1, 1);
plot(t, convAM);
title("Rectified Conventional AM message plot");
xlabel("time (seconds)");
ylabel("convAM(t)");

%Plotting Fourier Amplitude Plot of Rectified Conventional AM message signal
subplot(3, 1, 2);
plot(f, abs(convAMFT));
title("Fourier Amplitude Plot of Rectified Conventional AM message");
xlabel("Frequency (Hz)");
ylabel("|convAMFT(\omega)|");

%Plotting Fourier Phase Plot of Rectified Conventional AM message signal
subplot(3, 1, 3);
plot(f, unwrap(angle(convAMFT)));
title("Fourier Phase Plot of Rectified Conventional AM message");
xlabel("Frequency (Hz)");
ylabel("\angleconvAMFT(\omega)");

demodm = lowpass(convAM, 5000, fs); %Applying a lowpass filter to m(t)
demodM = fftshift(fft(demodm)); %Fourier Transform of the lowpass filtered signal

%Plotting the Demodulated message m(t)
figure;
subplot(3, 1, 1);
plot(t, demodm);
title("Demodulated message plot");
xlabel("time (seconds)");
ylabel("m(t)");

%Plotting Amplitude Plot of Demodulated message m(t)
subplot(3, 1, 2);
plot(f, abs(demodM));
title("Fourier Amplitude Plot of Demodulated message");
xlabel("Frequency (Hz)");
ylabel("|M(\omega)|");

%Plotting Phase Plot of Demodulated message m(t)
subplot(3, 1, 3);
plot(f, unwrap(angle(demodM)));
title("Fourier Phase Plot of Demodulated message");
xlabel("Frequency (Hz)");
ylabel("\angleM(\omega)");

demodm = downsample(demodm, 50); %Downsampling the demodulated message signal
sound(demodm, fs/50); %Playing the sound after downsampling the sampling frequency