%% Hilbert function for Question 8

%The function Hilbert takes in an input time domain signal, performs a
%Hilbert Transform on that signal, and then outputs it.
function [xt, xf] = Hilbert(x, f)

X = fftshift(fft(x)); %Shifted Fast fourier transform of a function in the time domain.
xf = X.*(-1j).*sign(f); %Calculates the function's Hilbert Transform
xt = ifft(fftshift(xf)); %Inverse Fast Fourier Transform to return to time domain.

end
