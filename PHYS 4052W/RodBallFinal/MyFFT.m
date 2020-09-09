function [y, f, Fs] = MyFFT(x, t)
L = length(x);
y = fft(x);
% Compute the two-sided spectrum P2. 
% Then compute the single-sided spectrum P1 
% based on P2 and the even-valued signal length L.

P2 = abs(y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

%Define the frequency domain f and calculate its range
P1(1)=0; %zero out DC offset
y = P1;
Fs = 1/(t(2)-t(1));
f = Fs*(0:(L/2))/L;
end
