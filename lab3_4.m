%% Given parameters, basic assumptions and derivations
f0 = 77e9;
fs = 10e6;
nChirps = 128;
N = 256;
mu = 25e6/1e-6;
c = 3e8;
Tcri = 50e-6;
Tc = nChirps/fs;
BW = mu*Tc;
dR = c/(2*BW);

Pt = 0.1;
G = 30;
lam = c/f0;
Ls = 1;
La_R = 1; 

t = 0:1/fs:Tcri-1/fs; %Time duration
Nsamps = length(t);

%% Target Parameters
Rk = [20 , 21, 20];        %meters
Vk = [30 , -12, 25];    %m/s
rcs = [1, 1, 1];                %radar cross section
nTgts = length(Rk);             %number of targets

%% zIF Signal
zIF = zeros(nChirps, Nsamps);
for kc = 1:nChirps
    for kt = 1:nTgts
        dist = Rk(kt)+Vk(kt)*(kc-1)*Tcri;
        tk = 2*dist/c;
        Pr = Pt*G^2*lam^2*rcs(kt)/((4*pi)^3*(dist)^4*Ls*La_R);
        
        zIF(kc,:) = zIF(kc,:) + (sqrt(Pr)/2)*exp(1i*(2*pi*f0*tk + 2*pi*mu*t*tk - pi*mu*(tk.^2)));
    end
end

%% Plotting
fRk         = (0:N-1)/N*fs;
Rng         = c/2/mu*fRk;
zIF_1Dfft   = fft(zIF,N,2);
zIF_1Dfft_1   = fft(zIF(1,:),N,2);
zIF_1Dfft_2   = fft(zIF(128,:),N,2);

fDk         = ((0:nChirps-1)/nChirps - 0.5)/Tcri;
Vel         = c/2/f0*fDk;
zIF_2Dfft  = fftshift(fft(zIF_1Dfft(:,2),nChirps,1),1);

figure;
plot(Rng,(abs(zIF_1Dfft_1)))
hold on
plot(Rng,(abs(zIF_1Dfft_2)))
xlabel('Range (m)')
ylabel('Amplitude')
grid ON
title('DFT1 map case-4')
legend('chirp 1','chirp 128')
hold off

figure;
plot(Vel,(abs(zIF_2Dfft)))
xlabel('Velocity (m/s)')
ylabel('Amplitude')
title('DFT2 map case-4')
grid ON

plotRDmap(zIF,fs,Tcri,c,f0,mu)

function plotRDmap(zIF,fs,Tcri,c,f0,mu)
    [nChirps,N] = size(zIF);
    fRk         = (0:N-1)/N*fs;
    Rng         = c/2/mu*fRk;
    zIF_1Dfft   = fft(zIF,N,2);
    fDk         = ((0:nChirps-1)/nChirps - 0.5)/Tcri;
    Vel         = c/2/f0*fDk;
    zIF_2Dfft   = fftshift(fft(zIF_1Dfft,nChirps,1),1);
    
    figure; imagesc(Rng,Vel,10*log10(abs(zIF_2Dfft)+eps)); colorbar;
    xlabel('Range (m)')
    ylabel('Velocity (m/s)')
    title('Range-Velocity Map case-4')
end