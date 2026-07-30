%% PSD of White, Pink, and Brown Noise
clear; close all; clc;

%% Parameters
Fs = 2;          % Sampling frequency (Hz)
N  = 2^18;          % Number of samples

%% Create noise generators
whiteGen = dsp.ColoredNoise( ...
    'Color','white', ...
    'SamplesPerFrame',N, ...
    'NumChannels',1);

pinkGen = dsp.ColoredNoise( ...
    'Color','pink', ...
    'SamplesPerFrame',N, ...
    'NumChannels',1);

brownGen = dsp.ColoredNoise( ...
    'Color','brown', ...
    'SamplesPerFrame',N, ...
    'NumChannels',1);

%% Generate signals
white = whiteGen();
pink  = pinkGen();
brown = brownGen();

%% Estimate PSD
window   = hamming(4096);
noverlap = 2048;
nfft     = 4096;

[Pw,f] = pwelch(white,window,noverlap,nfft,Fs);
[Pp,~] = pwelch(pink, window,noverlap,nfft,Fs);
[Pb,~] = pwelch(brown,window,noverlap,nfft,Fs);

%% Plot PSD
fig=figure('WindowState','maximized','Theme','Light');
loglog(f,Pw,'k','LineWidth',1.5); hold on;
loglog(f,Pp,'r','LineWidth',1.5);
loglog(f,Pb,'b','LineWidth',1.5);

grid off; axis tight;
xlim([0 Fs/2]);

xlabel('Frequency (Hz)');
ylabel('PSD (Power/Hz)');

legend('White (1/f^0)', ...
    'Pink (1/f)', ...
    'Brown (1/f^2)', ...
    'Location','southwest','box','off');
ax=gca;
ax.LineWidth=2;
ax.FontSize=18;

% Fs = 2;                      % Sampling frequency (Hz)
% 
% % Frequency vector from 10^-3 Hz to Nyquist (1 Hz)
% f = logspace(-3, log10(Fs/2), 1000);
% 
% % Normalized theoretical PSDs
% Swhite = ones(size(f));
% Spink  = 1./f;
% Sbrown = 1./f.^2;
% 
% figure;
% loglog(f,Swhite,'k','LineWidth',2); hold on;
% loglog(f,Spink,'r','LineWidth',2);
% loglog(f,Sbrown,'b','LineWidth',2);
% 
% grid on;
% box on;
% 
% xlabel('Frequency (Hz)');
% ylabel('Normalized PSD');
% 
% legend('White ($1$)', ...
%     'Pink ($1/f$)', ...
%     'Brown ($1/f^2$)', ...
%     'Interpreter','latex', ...
%     'Location','southwest');
% 
% xlim([1e-3 1]);
% ylim([1e0 1e6]);