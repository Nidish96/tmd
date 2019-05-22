% Demonstrate multisine excitation
% it is seen that interpolating a multisine destroys the frequency content

close all
clearvars
clc

srcpath = '~/src/matlab/PNLSS_v1_0';
addpath(genpath(srcpath));
figpath = '../fig/';

%%

fs = 1e3;
N = 2^13; % 8192    % Number of samples
M = round(0.2*N/2); % Last excited line. 20 % of Nyquist freq

% freq resolution
f0 = fs/N;
% fmax = 100;
% M = fmax/f0;
fmax=  M*f0;  % lines * freq_reso
fprintf('fmax: %0.3f\n', fmax);

% predictable random number. With obligatory references...
% https://xkcd.com/221/, https://dilbert.com/strip/2001-10-25
rng(10);
RMSu = 0.05; % Root mean square value for the input signal
R = 4;       % Number of phase realizations (one for validation and one for performance testing)
P = 5;       % Number of periods
kind = 'Full';           % 'Full','Odd','SpecialOdd', or 'RandomOdd': kind of multisine
[u,lines,non_exc_odd,non_exc_even] = fMultisine(N, kind, M, R); % Multisine signal, excited and detection lines
u = u/rms(u(:,1))*RMSu; % Scale multisine to the correct rms level
u = repmat(u,[1 1 P]);  % N x R x P
u = permute(u,[1 3 2]); % N x P x R
u_ext = u(:); % Concatenate the data: N*P*R x m. m is the number of inputs

% Normalized frequency vector
freq = (0:N-1)/N*fs;
time = (0:N*P*R-1)/fs;
% 0:1/fs:N*P/fs - 1/fs;

%% plot power spectrum, etc
figure; 
plot(time,u_ext)
% indicate periods
h1 = vline(time((1:R*P)*N),'--g');
% indicate realizations
h2 = vline(time((1:R)*N*P),'--k');
set([h1 h2],'LineWidth',0.5)
xlabel('time (s)')
ylabel('magnitude')
title(['Multisine: ' num2str(R) ' realizations of ' num2str(P) ' periods of ' num2str(N) ' samples per period'])
% export_fig([figpath, 'multisine_time.pdf'])

% plot only half spectrum(single sided spectrum)
dft = fft(u(1:end,1,1));
nt = length(dft)/2+1;
xdft = dft(1:nt);
% freq = freq(1:end);

figure; subplot(2,1,1)
plot(freq(1:nt),db(xdft),'-*')
xlabel('frequency (Hz)')
ylabel('magnitude (dB)')
title('FFT of one period of the multisine realizations')

subplot(2,1,2)
plot(freq(1:nt),angle(xdft),'-*')
xlabel('frequency (Hz)')
ylabel('phase (rad)')
title('FFT of one period of the multisine realizations')
% export_fig([figpath, 'multisine_freq.pdf'])

%% test of interpolation/freq content

% interpolate multisine. Do we destroy freq content? Yes!
u_interp1 = @(t) interp1(time, u_ext, t, 'nearest');
Nt = 3*N;
ttest = linspace(0,time(N),Nt);
freqtest = (0:Nt-1)/Nt*fs;
utest = u_interp1(ttest);
% plot only half spectrum(single sided spectrum)
dft = fft(utest);
nt = Nt/2+1;
xdft = dft(1:nt);

figure; subplot(2,1,1)
plot(freqtest(1:nt),db(xdft),'-*')
xlabel('frequency (Hz)')
ylabel('magnitude (dB)')
title(['FFT of one period of the interpolated multisine. N: ' num2str(Nt)])

subplot(2,1,2)
plot(freqtest(1:nt),angle(xdft),'-*')
xlabel('frequency (Hz)')
ylabel('phase (rad)')
title('FFT of one period of the interpolated multisine')
% export_fig([figpath, 'multisine_freq_interp1.pdf'])

%% test of downsampled multisine
% downsamp should match M. Ie. if downsamp is 4, then M should not be more
% than 25% of Nyquist freq.
downsamp = 4;  % 8192/4=2048
u1 = downsample(u(:,1,1),downsamp);

% plot only half spectrum(single sided spectrum)
dft = fft(u1);
nt = length(dft)/2+1;
xdft = dft(1:nt);

figure; subplot(2,1,1)
plot(freq(1:nt),db(xdft),'-*')
xlabel('frequency (Hz)')
ylabel('magnitude (dB)')
title('FFT of one period of the downsampled multisine realizations')

subplot(2,1,2)
plot(freqtest(1:nt),angle(xdft),'-*')
xlabel('frequency (Hz)')
ylabel('phase (rad)')
title('FFT of one period of the downsampled multisine')
% export_fig([figpath, 'multisine_freq_downs.pdf'])