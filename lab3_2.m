clear;
%%fmcw transmitted churp
%given 
Fo = 77*10^6;
Fs = 150*10^6;
mu = 25 *10^13; % MHz/us
BW = 10^9;
Tc = BW/mu;
Tcri = 10*Tc;
nchirps = 70;
G =30;
%implications
t = 1/Fs : 1/Fs : Tc;

%%received intermediate frequency(=>In harware)
%Case1
figure(1);
nTargets = 1;
R_targets = [30];
vel_targets = [100];
rcs = [3];
y = received_if_signal_generator(nTargets,R_targets,vel_targets,rcs,nchirps,Fo,mu,Fs,Tcri,G);
plotRDmap(y,Fs,Tcri,3*10^8,Fo,mu);

%Case2  
figure(2)
nTargets = 3;
R_targets = [30 42 51];
vel_targets = [10 13 21];
rcs = [3 6 0.1];
y = received_if_signal_generator(nTargets,R_targets,vel_targets,rcs,nchirps,Fo,mu,Fs,Tcri,G);
plotRDmap(y,Fs,Tcri,3*10^8,Fo,mu);


function y = received_if_signal_generator(ntargets,R_tgts,vel_tgts,rcs,nchirps,fo,mu,Fs,Tcri,G)
    c = 3*10^8;
    Rmax = 250;
    r_tgts = R_tgts;
    y = [];
    t = 0 : 1/Fs : Tcri/10 - 1/Fs; 
    for i = 1:nchirps
        x = 0;
        for j = 1:ntargets
            tk = 2*r_tgts(j)/c;
            lambda = c/fo;
            bk = G^2*lambda^2*rcs(j)/((4*pi)^3*(r_tgts(j) + vel_tgts(j)*(nchirps-1)*Tcri)^4);
            x = x + sqrt(bk) * exp(1j*(2*pi*fo*tk + 2*pi*mu*tk*t - pi*mu*(tk^2)));
        end
        y = [y;x];
        for k = ntargets
            r_tgts(k) = r_tgts(k) + vel_tgts(k)*Tcri;
        end
    end
end

function y = integration(t,fs)
    y = [];
    previous = 0;
    for i = t
        y = [y previous+i*(1/fs)];
        previous = previous + i*(1/fs);
    end
end
