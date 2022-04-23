% Radar Parameters

c = 3e8;
fo = 1.5e9;
Pt = 5;
GdB = 30;
G = 10^(GdB/10);
lambda = c/fo;
Ls = 10^(1.5/10);
La_R = 1;
f_if = 20e6;
B = 20e6;
Tp = 1e-6;
Rmax = 10e3;
Tpri = (2/c)*Rmax + Tp;
mu = B/Tp;
fs = 4*(f_if + B);
L = 256;


%targer Parameters

R_tgt = [2.5e3 8e3 3.5e3];
v_tgt = [25 60 -15]*1000/3600;
rcs = 10.^([10 20 15]/10);


t = 0:1/fs:Tpri - Tp-1/fs;
y = zeros(L,length(t));

for kl = 1:L
    for kt = 1:length(R_tgt)
        tau_kl = 2/c*(R_tgt(kt) + (kl-1)*Tpri*v_tgt(kt));
        Pr_k = Pt*lambda^2*G^2*rcs(kt)/((4*pi)^3*R_tgt(kt)^4*Ls*La_R);
        bk = sqrt(Pr_k);
        y(kl,:) = y(kl,:) + bk/2*rectpuls(t-tau_kl-Tp/2,Tp).*exp(1j*(2*pi*f_if*t-2*pi*fo*tau_kl+pi*mu*(t-tau_kl).^2));
    end
end

figure;
plot(c*t/2,real(y(1,:)));
xlabel('Range(m)');
ylabel( 'Amplitude');