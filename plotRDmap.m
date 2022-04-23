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