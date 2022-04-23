c = 3 * 10^8; %speed of light in m/s


%RADAR parameters given
f0 = 77 * 10^9; %starting frequency, given as 77 GHz, in Hz
fs = 10 * 10^6; %starting frequency, given as 10 MHz, in Hz
N = 256; %number of samples in one chirp transmitted
L = 128; %number of chirps transmitted
meu = 25 * 10^12; %chirp-rate, given as 25 MHz/us, in Hz/s
Tcri = 50 * 10^-6; %chirp repetition interval, given as 50 us, in s
P = 16; %number of receiver sensors 

%RADAR parameters assumed
Pt = 1; %transmit power, in W 
G = 10^(30 / 10); %gain, given as 30 dB
Ls = 1; %system loss
La_R = 1; %atmospheric loss

%RADAR parameters determined
lambda = c / f0; %wavelength, in m
Tc = N * (1 / fs); %chirp duration, in s
d = 0.5 * lambda; %distance between two receiver sensors

%target parameters given
R1 = [20, 31, 40] ; %range, in m
R2 = [20, 31, 22] ; %range, in m
R3 = [20, 31, 21] ; %range, in m
v1 = [15, -6, 12.5]; %velocity, in m/s
v2 = [15, -6, 12.5]; %velocity, in m/s
v3 = [15, -6, 15]; %velocity, in m/s
angle = [20, -20, 30]; %direction of arrival 
K = length(R1); %number of targets

%target parameters assumed
sigma = [2, 2, 2]; %RADAR cross-section, in m^2 

t = 0:1 / fs:Tc - 1 / fs; %time to transmit one chirp, in s
N = length(t); %number of samples in one chirp transmitted

%generating the down-converted received signal for the 1 st case
Rx_downconverted_p_l = zeros(P, L, N);
for p = 1:P
    for l = 1:L
        for k = 1:K
            R_k_l = R1(k) + v1(k) * (l - 1) * Tcri; %distance travelled by lth chirp transmitted before it is received, in m; assuming positive velocity corresponds to the target moving away from the RADAR
            Pr_k_l = (Pt * lambda^2 * G^2 * sigma(k)) / ((4 * pi)^3 * R_k_l^4 * Ls * La_R); %power of the lth chirp received, due to the kth target, in W
            tau_k_l = 2 * R_k_l / c; %delay in receiving the lth chirp transmitted, due to the kth target, in s    
            Rx_downconverted_p_l(p, l, 1:N) = Rx_downconverted_p_l(p, l, 1:N) + reshape((sqrt(Pr_k_l) / 2) * exp(1j * (2 * pi * meu * tau_k_l * t + 2 * pi * f0 * tau_k_l + 2 * pi * (p - 1) * (d / lambda) * sind(angle(k)))), 1, 1, N); %ignoring the noise component 
        end
    end
    disp("In the 1 st case, amplitude of the down-converted chirp-" + l + ", received by the receiver-" + p + ", due to target-" + k + ", is: " + sqrt(Pr_k_l) / 2); %amplitude of the down-converted received signal
    disp("In the 1 st case, amplitude of the down-converted chirp-" + l + ", received by the receiver-" + p + ", due to target-" + k + ", after taking FFT once is: " + N * sqrt(Pr_k_l) / 2); %amplitude of the down-converted received signal after taking FFT once 
    disp("In the 1 st case, amplitude of the down-converted chirp-" + l + ", received by the receiver-" + p + ", due to target-" + k + ", after taking FFT twice is: " + L * N * sqrt(Pr_k_l) / 2); %amplitude of the down-converted received signal after taking FFT twice
    disp("In the 1 st case, amplitude of the down-converted chirp-" + l + ", received by the receiver-" + p + ", due to target-" + k + ", after taking FFT thrice is: " + P * L * N * sqrt(Pr_k_l) / 2); %amplitude of the down-converted received signal after taking FFT thrice
end

fR = 0:fs / N:fs - fs/N;
range = (c / (2 * meu)) * fR;

Rx_downconverted_1_l(1:L, 1:N) = Rx_downconverted_p_l(1, 1:L, 1:N);
Rx_downconverted_1D_FFT = fft(Rx_downconverted_1_l, N, 2);
        
%plotting the DFT for the 1 st chirp received by the 1 st receiver in the 1 st case
figure; 
subplot(2, 1, 1);
plot(range, abs(Rx_downconverted_1D_FFT(1, :)),'r');
grid on;
xlabel('range (m)');
ylabel('amplitude');
title('DFT plot-1 for the 1 st chirp received by the 1 st receiver');

%plotting the DFT for the last chirp received by the 1 st receiver in the 1 st case
subplot(2, 1, 2);
plot(range, abs(Rx_downconverted_1D_FFT(L, :)),'r');
grid on;
xlabel('range (m)');
ylabel('amplitude');
grid on;
title('DFT plot-1 for the last chirp received by the 1 st receiver');

%plotting the range-velocity map for 1 st receiver in the 1 st case
plotRDmap(Rx_downconverted_1_l, fs, Tcri, f0, meu, '-1');

fD = -0.5 / Tcri:1 / (L * Tcri):0.5 / Tcri - 1 / (L * Tcri);
velocity = (c / 2) * (fD / f0);

ftheta = -(d / lambda): (d / lambda) / (P / 2): (d / lambda) - ((d / lambda) / (P / 2));
theta = asind((lambda / d) * ftheta);

Rx_downconverted_1D_FFT = fft(Rx_downconverted_p_l, N, 3);
Rx_downconverted_2D_FFT = fftshift(fft(Rx_downconverted_1D_FFT, L, 2), 2);

[~, Rindex] = min(abs(range - R1(1)));
[~, vindex] = min(abs(velocity - v1(1)));

Rx_downconverted_3D_FFT = fftshift(fft(Rx_downconverted_2D_FFT(:, vindex, Rindex), P, 1), 1);

%plotting the angular spectrum for the 1 st chirp of 1 st target in the 1 st case
figure; 
subplot(3, 1, 1);
plot(theta, abs(Rx_downconverted_3D_FFT),'r');
grid on;
xlabel('theta (degree)');
ylabel('amplitude');
title('Angular plot-1 for the 1 st target');

[~, Rindex] = min(abs(range - R1(2)));
[~, vindex] = min(abs(velocity - v1(2)));

Rx_downconverted_3D_FFT = fftshift(fft(Rx_downconverted_2D_FFT(:, vindex, Rindex), P, 1), 1);

%plotting the angular spectrum for the 1 st chirp of 2 nd target in the 1 st case
subplot(3, 1, 2);
grid on;
plot(theta, abs(Rx_downconverted_3D_FFT),'r');
grid on;
xlabel('theta (degree)');
ylabel('amplitude');
title('Angular plot-1 for the 2 nd target');

[~, Rindex] = min(abs(range - R1(3)));
[~, vindex] = min(abs(velocity - v1(3)));

Rx_downconverted_3D_FFT = fftshift(fft(Rx_downconverted_2D_FFT(:, vindex, Rindex), P, 1), 1);

%plotting the angular spectrum for the 1 st chirp of 3 rd target in the 1 st case
subplot(3, 1, 3);
grid on;
plot(theta, abs(Rx_downconverted_3D_FFT),'r');
grid on;
xlabel('theta (degree)');
ylabel('amplitude');
title('Angular plot-1 for the 3 rd target');

%generating the down-converted received signal for the 2 nd case
Rx_downconverted_p_l = zeros(P, L, N);
for p = 1:P
    for l = 1:L
        for k = 1:K
            R_k_l = R2(k) + v2(k) * (l - 1) * Tcri; %distance travelled by lth chirp transmitted before it is received, in m; assuming positive velocity corresponds to the target moving away from the RADAR
            Pr_k_l = (Pt * lambda^2 * G^2 * sigma(k)) / ((4 * pi)^3 * R_k_l^4 * Ls * La_R); %power of the lth chirp received, due to the kth target, in W
            tau_k_l = 2 * R_k_l / c; %delay in receiving the lth chirp transmitted, due to the kth target, in s    
            Rx_downconverted_p_l(p, l, 1:N) = Rx_downconverted_p_l(p, l, 1:N) + reshape((sqrt(Pr_k_l) / 2) * exp(1j * (2 * pi * meu * tau_k_l * t + 2 * pi * f0 * tau_k_l + 2 * pi * (p - 1) * (d / lambda) * sind(angle(k)))), 1, 1, N); %ignoring the noise component 
        end
    end
end

fR = 0:fs / N:fs - fs/N;
range = (c / (2 * meu)) * fR;

Rx_downconverted_1_l(1:L, 1:N) = Rx_downconverted_p_l(1, 1:L, 1:N);
Rx_downconverted_1D_FFT = fft(Rx_downconverted_1_l, N, 2);
        
%plotting the DFT for the 1 st chirp received by the 1 st receiver in the 2 nd case
figure; 
subplot(2, 1, 1);
grid on;
plot(range, abs(Rx_downconverted_1D_FFT(1, :)),'r');
grid on;
xlabel('range (m)');
ylabel('amplitude');
title('DFT plot-2 for the 1 st chirp received by the 1 st receiver');

%plotting the DFT for the last chirp received by the 1 st receiver in the 2 nd case
subplot(2, 1, 2);
grid on;
plot(range, abs(Rx_downconverted_1D_FFT(L, :)),'r');
grid on;
xlabel('range (m)');
ylabel('amplitude');
title('DFT plot-2 for the last chirp received by the 1 st receiver');

%plotting the range-velocity map for 1 st receiver in the 2 nd case
plotRDmap(Rx_downconverted_1_l, fs, Tcri, f0, meu, '-2');

fD = -0.5 / Tcri:1 / (L * Tcri):0.5 / Tcri - 1 / (L * Tcri);
velocity = (c / 2) * (fD / f0);

ftheta = -(d / lambda): (d / lambda) / (P / 2): (d / lambda) - ((d / lambda) / (P / 2));
theta = asind((lambda / d) * ftheta);

Rx_downconverted_1D_FFT = fft(Rx_downconverted_p_l, N, 3);
Rx_downconverted_2D_FFT = fftshift(fft(Rx_downconverted_1D_FFT, L, 2), 2);

[~, Rindex] = min(abs(range - R2(1)));
[~, vindex] = min(abs(velocity - v2(1)));

Rx_downconverted_3D_FFT = fftshift(fft(Rx_downconverted_2D_FFT(:, vindex, Rindex), P, 1), 1);

%plotting the angular spectrum for the 1 st chirp of 1 st target in the 2 nd case
figure; 
subplot(3, 1, 1);
grid on;
plot(theta, abs(Rx_downconverted_3D_FFT),'r');
grid on;
xlabel('theta (degree)');
ylabel('amplitude');
title('Angular plot-2 for the 1 st target');

[~, Rindex] = min(abs(range - R2(2)));
[~, vindex] = min(abs(velocity - v2(2)));

Rx_downconverted_3D_FFT = fftshift(fft(Rx_downconverted_2D_FFT(:, vindex, Rindex), P, 1), 1);

%plotting the angular spectrum for the 1 st chirp of 2 nd target in the 2 nd case
subplot(3, 1, 2);
grid on;
plot(theta, abs(Rx_downconverted_3D_FFT),'r');
grid on;
xlabel('theta (degree)');
ylabel('amplitude');
title('Angular plot-2 for the 2 nd target');

[~, Rindex] = min(abs(range - R2(3)));
[~, vindex] = min(abs(velocity - v2(3)));

Rx_downconverted_3D_FFT = fftshift(fft(Rx_downconverted_2D_FFT(:, vindex, Rindex), P, 1), 1);

%plotting the angular spectrum for the 1 st chirp of 3 rd target in the 2 nd case
subplot(3, 1, 3);
grid on;
plot(theta, abs(Rx_downconverted_3D_FFT),'r');
grid on;
xlabel('theta (degree)');
ylabel('amplitude');
title('Angular plot-2 for the 3 rd target');

%generating the down-converted received signal for the 3 rd case
Rx_downconverted_p_l = zeros(P, L, N);
for p = 1:P
    for l = 1:L
        for k = 1:K
            R_k_l = R3(k) + v3(k) * (l - 1) * Tcri; %distance travelled by lth chirp transmitted before it is received, in m; assuming positive velocity corresponds to the target moving away from the RADAR
            Pr_k_l = (Pt * lambda^2 * G^2 * sigma(k)) / ((4 * pi)^3 * R_k_l^4 * Ls * La_R); %power of the lth chirp received, due to the kth target, in W
            tau_k_l = 2 * R_k_l / c; %delay in receiving the lth chirp transmitted, due to the kth target, in s    
            Rx_downconverted_p_l(p, l, 1:N) = Rx_downconverted_p_l(p, l, 1:N) + reshape((sqrt(Pr_k_l) / 2) * exp(1j * (2 * pi * meu * tau_k_l * t + 2 * pi * f0 * tau_k_l + 2 * pi * (p - 1) * (d / lambda) * sind(angle(k)))), 1, 1, N); %ignoring the noise component 
        end
    end
end

fR = 0:fs / N:fs - fs/N;
range = (c / (2 * meu)) * fR;

Rx_downconverted_1_l(1:L, 1:N) = Rx_downconverted_p_l(1, 1:L, 1:N);
Rx_downconverted_1D_FFT = fft(Rx_downconverted_1_l, N, 2);
        
%plotting the DFT for the 1 st chirp received by the 1 st receiver in the 3 rd case
figure; 
subplot(2, 1, 1);
grid on;
plot(range, abs(Rx_downconverted_1D_FFT(1, :)),'r');
grid on;
xlabel('range (m)');
ylabel('amplitude');
title('DFT plot-3 for the 1 st chirp received by the 1 st receiver');

%plotting the DFT for the last chirp received by the 1 st receiver in the 3 rd case
subplot(2, 1, 2);
grid on;
plot(range, abs(Rx_downconverted_1D_FFT(L, :)),'r');
grid on;
xlabel('range (m)');
ylabel('amplitude');
title('DFT plot-3 for the last chirp received by the 1 st receiver');

%plotting the range-velocity map for 1 st receiver in the 3 rd case
plotRDmap(Rx_downconverted_1_l, fs, Tcri, f0, meu, '-3');

fD = -0.5 / Tcri:1 / (L * Tcri):0.5 / Tcri - 1 / (L * Tcri);
velocity = (c / 2) * (fD / f0);

ftheta = -(d / lambda): (d / lambda) / (P / 2): (d / lambda) - ((d / lambda) / (P / 2));
theta = asind((lambda / d) * ftheta);

Rx_downconverted_1D_FFT = fft(Rx_downconverted_p_l, N, 3);
Rx_downconverted_2D_FFT = fftshift(fft(Rx_downconverted_1D_FFT, L, 2), 2);

[~, Rindex] = min(abs(range - R3(1)));
[~, vindex] = min(abs(velocity - v3(1)));

Rx_downconverted_3D_FFT = fftshift(fft(Rx_downconverted_2D_FFT(:, vindex, Rindex), P, 1), 1);

%plotting the angular spectrum for the 1 st chirp of 1 st target in the 3 rd case
figure; 
subplot(3, 1, 1);
grid on;
plot(theta, abs(Rx_downconverted_3D_FFT),'r');
grid on;
xlabel('theta (degree)');
ylabel('amplitude');
title('Angular plot-3 for the 1 st target');

[~, Rindex] = min(abs(range - R3(2)));
[~, vindex] = min(abs(velocity - v3(2)));

Rx_downconverted_3D_FFT = fftshift(fft(Rx_downconverted_2D_FFT(:, vindex, Rindex), P, 1), 1);

%plotting the angular spectrum for the 1 st chirp of 2 nd target in the 3 rd case
subplot(3, 1, 2);
grid on;
plot(theta, abs(Rx_downconverted_3D_FFT),'r');
grid on;
xlabel('theta (degree)');
ylabel('amplitude');
title('Angular plot-3 for the 2 nd target');

[~, Rindex] = min(abs(range - R3(3)));
[~, vindex] = min(abs(velocity - v3(3)));

Rx_downconverted_3D_FFT = fftshift(fft(Rx_downconverted_2D_FFT(:, vindex, Rindex), P, 1), 1);

%plotting the angular spectrum for the 1 st chirp of 3 rd target in the 3 rd case
subplot(3, 1, 3);
grid on;
plot(theta, abs(Rx_downconverted_3D_FFT),'r');
grid on;
xlabel('theta (degree)');
ylabel('amplitude');
title('Angular plot-3 for the 3 rd target');

function plotRDmap(z_downconverted, fs, Tcri, f0, meu, plot)
    c = 3 * 10^8;
    [L, N] = size(z_downconverted);
    fR = 0:fs / N:fs - fs/N;
    R = (c / (2 * meu)) * fR;
    z_downconverted_1D_FFT = fft(z_downconverted, N, 2);
    fD = -0.5 / Tcri:1 / (L * Tcri):0.5 / Tcri - 1 / (L * Tcri);
    v = (c / 2) * (fD / f0);
    z_downconverted_2D_FFT = fftshift(fft(z_downconverted_1D_FFT, L, 1), 1);
    
    figure; 
    imagesc(R, v, 10 * log10(abs(z_downconverted_2D_FFT) + eps));
    colorbar;
    xlabel('range (m)');
    ylabel('velocity (m/s)');
    title(['Range-velocity map', plot]);
end