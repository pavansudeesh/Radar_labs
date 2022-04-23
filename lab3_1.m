clear all;
% Radar parameters
f0          = 77e9;
c           = 3e8;
lam         = c/f0;
mu          = 25e12;

fs          = 10e6;
N           = 256;
Tc          = N/fs;
Tcri        = 50e-6;
N           = round(Tc*fs);
idleTime    = Tcri - Tc;
Pt          = 0.1;
G           = 10^(30/10);
Ls          = 1;
La_R        = 1;

nChirps     = 128;
Msens       = 8;
d           = lam/2;
MsensPosn   = d*(0:Msens-1).';
zIF         = zeros(Msens,nChirps,N);

% Target parameters
Rk          = [20, 31, 40];
vk          = [15, -6, 12.5];
Dirk        = [20 -22 30];
rcs         = [3,4, 7];
nTgts       = length(Rk);

t     = (0:N-1)/fs;
for kc = 1:nChirps
    for kt  = 1:nTgts
        tk  = 2/c*(Rk(kt) + vk(kt)*(kc-1)*Tcri);
        Pr  = Pt*G^2*lam^2*rcs(kt)/((4*pi)^3*(Rk(kt)+vk(kt)*(kc-1)*Tcri)^4*Ls*La_R);
        
        for kp = 1:Msens
            ph = exp(1j*2*pi/lam*MsensPosn(kp)*sind(Dirk(kt)));
            sens_sig = sqrt(Pr)*rectpuls(t - Tc/2 - tk,Tc - tk).*exp(1j*(2*pi*f0*tk + 2*pi*mu*tk*t))*ph;
            zIF(kp,kc,:) = zIF(kp,kc,:) + reshape(sens_sig,[1,1,N]);
        end
%         spv = exp(1j*2*pi/lam*MsensPosn*sind(Dirk(kt)));
%         zIF(:,kc,:) = zIF(:,kc,:) + reshape(sqrt(Pr)spv(rectpuls(t - Tc/2 - tk,Tc - tk).exp(1j(2*pi*f0*tk + 2*pi*mu*tk*t))),[Msens , 1, N]);
    end
end

zIF = awgn(zIF,10,10*log10(var(zIF(:))));

% First DFT/FFT
k_ind       = (0:N-1);
fR_ind      = k_ind/N*fs;
R_ind       = c*fR_ind/(2*mu);
zIF_1F      = fft(zIF,N,3);
figure; 
plot(R_ind,abs(squeeze(zIF_1F(1,1,:)))-max(abs(zIF_1F(1,1,:))))
xlim([15 45])
xlabel('Range (m)')
% 
% Second DFT/FFT
zIF_2F      = fftshift(fft(zIF_1F,nChirps,2),2);
fD_idcs     = ((0:nChirps-1)/nChirps - 0.5)/Tcri;
v_idcs      = c*fD_idcs/(2*f0);
[aa, idx1]    = min(abs(R_ind-Rk(1)));

figure; plot(v_idcs, abs(squeeze(zIF_2F(1,:,idx1))));
xlabel('Velocity (m/s)')

plotRDmap(squeeze(zIF(2,:,:)),fs,Tcri,c,f0,mu)

ftta        = ((0:Msens-1)/Msens) - 0.5;
Tta_idcs    = asind(lam*ftta/d);
Tta_Interp  = linspace(-80, 80, 41);

zIF_3F = zeros(Msens,nTgts);
zIF_3F_interp = zeros(length(Tta_Interp),nTgts);
figure;

for kt = 1:nTgts
    [aa, ridx1]    = min(abs(R_ind-Rk(kt)));
    Ridx{kt}       = ridx1;
    
    [aa, vidx1]    = min(abs(v_idcs-vk(kt)));
    Vidx{kt}       = vidx1;
    
    zIF_3F(:,kt) = fftshift(fft(zIF_2F(:,Vidx{kt},Ridx{kt})));
    zIF_3F_interp(:,kt) = interp1(Tta_idcs,abs(zIF_3F(:,kt)),Tta_Interp,'v5cubic');
    subplot(nTgts,1,kt); hold on;
    plot(Tta_idcs,abs(zIF_3F(:,kt)))
    plot(Tta_Interp,zIF_3F_interp(:,kt),'r')
end