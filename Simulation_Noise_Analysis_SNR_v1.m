%% Simulation_Noise_Analysis_SNR.m
clear; close all; clc;
addpath('functions');

%% ============================================================
%  PARAMETERS  (network identical to Simulation.m)
% ============================================================
fs    = 2;     % sampling frequency (Hz)
nfft  = 512;    % FFT resolution (scalars are nfft-robust -> modest nfft)
nobs  = 10000;    % time samples
n_osc_cycles = 2.5;  % number of coupling on/off cycles over the whole record
DC    = 50;      % square-wave duty cycle (%)

% coupling oscillation frequency, adapted to nobs and fs so the number of
% cycles across the record stays fixed regardless of nobs/fs:
%   n_osc_cycles = f_osc * (nobs/fs)  ->  f_osc = n_osc_cycles * fs / nobs
f_osc = n_osc_cycles * fs / nobs;   % (Hz)
c1max = 0.8;       % coupling amplitude
M     = 7;       % number of VAR processes
c_ff  = 0.98;    % RLS forgetting factor
nReal = 1;     % number of realizations
base_seed = 4000;% realization seeds: base_seed + r

snr_list = [100, 10, 5, 1];  
nSNR     = numel(snr_list);

% --- output folder for saved figures ---
outRoot        = fullfile(pwd, 'Noise_Analysis_Figures_v2_4000_fs2');
closeAfterSave = true;        % close each figure after saving
saveRes        = 600;         % PNG resolution (dpi)
if ~exist(outRoot, 'dir'), mkdir(outRoot); end

k = 1:nobs;
t = k / fs;      % time axis (s)

% Time-varying coupling
c1 = -c1max * square(2*pi*f_osc*t, DC);
c1(c1 == -c1max) = 0;
c2 = 1 - c1;    

t_in = 7;        % first valid spectrogram column (skip RLS transient)

% Measure display labels
meas_lbl = {'$F_{Y;X_1}$', '$F_{Y;X_2}$', '$F_{Y;X_1,X_2}$', '$I_{Y;X_1;X_2}$'};

% Custom colormap (shared with Simulation.m) for the spectrogram
cmap = load('colmap.mat', 'VRVmap');

% Conditions to test (clean baseline + three noise colors)
noiseTypes = {'clean', 'white', 'pink', 'brown'};
noiseName  = {'Clean (no added noise)', 'White noise (flat)', ...
              'Pink noise (1/f)', 'Brown noise (1/f^2)'};
nN = numel(noiseTypes);

fprintf('=== Simulation_Noise_Analysis_SNR.m ===\n');
fprintf('  fs=%d Hz | nobs=%d | M=%d | ff=%.2f | %d realizations\n', ...
        fs, nobs, M, c_ff, nReal);
fprintf('  SNR sweep (linear ratio): %s\n', mat2str(snr_list));
fprintf('  figures -> %s\n\n', outRoot);

%% ============================================================
%  1) TRUE time-varying VAR  (noise-independent -> build once)
% ============================================================
par.poles{1} = [0.9 0.35];
par.poles{4} = [0.9 0.10];
par.poles{6} = [0.9 0.10];
par.poles{7} = [];

AmT = []; SuT = [];
for n = 1:length(t)
    C = c1(n);
    par.coup = [1 2 1 0.8;
                1 3 1 C;
                4 5 2 C;
                6 5 1 0.8;
                5 7 2 C;
                2 7 2 0.8];
    par.Su = ones(1, M);
    [Am1, Su1] = theoreticalVAR(M, par);

    Am1 = Am1';
    Am1 = reshape(Am1, [M, M, size(Am1,2)/M]);

    AmT(:,:,:,n) = Am1;
    SuT(:,:,n)   = Su1;
end
p = size(AmT, 3);                         % true (and estimation) model order
for n = 1:nobs
    SuT(:,:,n) = (SuT(:,:,n) + SuT(:,:,n)')/2 + 1e-12*eye(M);
end

%% ============================================================
%  1b) THEORETICAL measures from the TRUE TV-VAR  (ground truth)
%      Noise-independent and estimation-free -> compute once.
%      These exact curves replace the old coupling reference and
%      are drawn in BLACK on every trace panel.
% ============================================================
Fy1_th  = zeros(1,nobs);  Fy2_th = zeros(1,nobs);
Fy12_th = zeros(1,nobs);  Iy_th  = zeros(1,nobs);
fy1_th  = zeros(nfft,nobs);  fy2_th = zeros(nfft,nobs);
fy12_th = zeros(nfft,nobs);  iy_th  = zeros(nfft,nobs);

wsth = warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
for n = 1:nobs
    Am_n = reshape(AmT(:,:,:,n), M, M*p);   % [A(1) ... A(p)] for sir_VARspectra
    S_n  = sir_VARspectra(Am_n, SuT(:,:,n), nfft, fs);
    out  = iispectral(S_n, 7, [1 2 3], [4 5 6]);
    Fy1_th(n)  = out.Fy_x1;    Fy2_th(n)  = out.Fy_x2;
    Fy12_th(n) = out.Fy_x1x2;  Iy_th(n)   = out.Iy_x1_x2;
    fy1_th(:,n)  = out.fy_x1;   fy2_th(:,n)  = out.fy_x2;
    fy12_th(:,n) = out.fy_x1x2; iy_th(:,n)   = out.iy_x1_x2;
end
warning(wsth);

% sanitize (near-singular spectra can yield complex/Inf)
TR_th = {Fy1_th, Fy2_th, Fy12_th, Iy_th};
for m = 1:4
    v = real(TR_th{m});  v(~isfinite(v)) = NaN;  TR_th{m} = v;
end

%% ============================================================
%  2) COMPUTE — all SNR x noise colors x realizations
% ============================================================
% reduced storage, indexed [SNR, color]
TRmed = cell(nSNR,nN); TRq1 = cell(nSNR,nN); TRq3 = cell(nSNR,nN); SPmed = cell(nSNR,nN);
f_vec = [];

tStart = tic;
for si = 1:nSNR
    snr_lin = snr_list(si);
    fprintf('############  SNR = %g (%.1f dB)  ############\n', snr_lin, 10*log10(snr_lin));

    for e = 1:nN
        fprintf('--- %s  (%d realizations) ---\n', noiseName{e}, nReal);

        trace = {zeros(nReal,nobs), zeros(nReal,nobs), ...
                 zeros(nReal,nobs), zeros(nReal,nobs)};         % {Fy1,Fy2,Fy12,Iy}
        SP    = {zeros(nReal,nfft,nobs), zeros(nReal,nfft,nobs), ...
                 zeros(nReal,nfft,nobs), zeros(nReal,nfft,nobs)};

        for r = 1:nReal
            rng(base_seed + r);            % same clean realization across colours/SNR

            Xgen = var_to_tsdata_nonstat(AmT, SuT, 1);
            Xsim = Xgen(:,:,1);                                  % M x nobs (clean)
            Xn   = addColoredNoise(Xsim, noiseTypes{e}, snr_lin);% corrupt BEFORE est.

            X = permute(Xn, [1 3 2]);
            X = zscore(X, 0, 3);           % normalize data

            [F1,F2,F12,Iy, s1,s2,s12,siy, f_vec] = ...
                runPipeline(X, c_ff, p, nfft, fs, nobs);

            Fm = {F1,F2,F12,Iy};
            sp = {s1,s2,s12,siy};
            for m = 1:4
                trace{m}(r,:) = Fm{m};
                SP{m}(r,:,:)  = sp{m};
            end

            if mod(r,20) == 0 || r == 1
                fprintf('  realization %3d/%d  [elapsed %.0fs]\n', r, nReal, toc(tStart));
            end
        end

        % --- sanitize: near-singular spectra can make fy complex/Inf ---
        for m = 1:4
            v = real(trace{m});  v(~isfinite(v)) = NaN;  trace{m} = v;
            s = real(SP{m});     s(~isfinite(s)) = NaN;  SP{m}    = s;
        end

        % --- reduce over realizations ---
        TRmed{si,e} = cell(1,4); TRq1{si,e} = cell(1,4);
        TRq3{si,e}  = cell(1,4); SPmed{si,e} = cell(1,4);
        for m = 1:4
            TRmed{si,e}{m} = median(trace{m}, 1, 'omitnan');        % 1 x nobs
            TRq1{si,e}{m}  = prctile(trace{m}, 25, 1);
            TRq3{si,e}{m}  = prctile(trace{m}, 75, 1);
            SPmed{si,e}{m} = squeeze(median(SP{m}, 1, 'omitnan'));  % nfft x nobs
        end
        fprintf('\n');
    end
end

%% ============================================================
%  FIGURE STYLING
% ============================================================
col4 = [0.00 0.20 0.90;    % Fy1   blue
        0.00 0.75 0.40;    % Fy2   green
        0.90 0.50 0.00;    % Fy12  orange
        0.85 0.05 0.05];   % Iy    red

% SINGLE spectrogram color scale shared by ALL measures (incl. Iy) so every
% panel is directly comparable on one colorbar.
d = [];
for m = 1:4
    for si = 1:nSNR
        for e = 1:nN
            d = [d; SPmed{si,e}{m}(:)]; %#ok<AGROW>
        end
    end
end
d = d(isfinite(d));
cl = prctile(abs(d), 99);                 % symmetric limit from the data
clim_all = [-cl cl];
clims_sp = {clim_all, clim_all, clim_all, clim_all};

% SINGLE y-limit shared by ALL measures, SNR & colors -> every trace panel
% (Fy1, Fy2, Fy12, Iy) is drawn on one common scale for direct comparison.
d = [];
for m = 1:4
    d = [d; TR_th{m}(:)];                  % include theoretical (black) curve
    for si = 1:nSNR
        for e = 1:nN
            d = [d; TRq1{si,e}{m}(:); TRq3{si,e}{m}(:)]; %#ok<AGROW>
        end
    end
end
d  = d(isfinite(d));
lo = prctile(d, 1);  hi = prctile(d, 99);
mg = max(0.1*(hi - lo), 0.05);
ylim_all = [lo - mg, hi + mg];
ylims_sc = repmat({ylim_all}, 1, 4);      % same scale for every measure

% SNR color ramp for the overlay (high SNR = light, low SNR = saturated red)
base = col4(4,:);
frac = linspace(0.30, 1.0, nSNR)';
cols_snr = 1 - frac .* (1 - base);

%% ============================================================
%  (a) PER-COLOUR TRACE GRID   (rows = SNR, cols = measure)
% ============================================================
for e = 1:nN
    fh = figure('Units','inches','Position',[0.5 0.5 13 9.5], ...
                'Name',sprintf('Traces | %s', noiseName{e}), 'Theme','Light','WindowState','maximized');
    tl = tiledlayout(nSNR, 4, 'TileSpacing','compact', 'Padding','compact');
    hB = []; hM = []; hR = [];
    for si = 1:nSNR
        for m = 1:4
            ax = nexttile; hold(ax,'on');
            md = TRmed{si,e}{m};  q1 = TRq1{si,e}{m};  q3 = TRq3{si,e}{m};
            bad = ~isfinite(q1) | ~isfinite(q3);
            q1(bad) = md(bad);  q3(bad) = md(bad);
            q1(~isfinite(q1)) = 0;  q3(~isfinite(q3)) = 0;
            tint = 1 - 0.65*(1 - col4(m,:));
            hb = fill([t fliplr(t)], [q3 fliplr(q1)], tint, ...
                      'FaceAlpha',0.55, 'EdgeColor','none');
            hm = plot(t, md, 'Color', col4(m,:), 'LineWidth', 1.3);
            hr = plot(t, TR_th{m}, 'k-', 'LineWidth', 1.4);   % theoretical (ground truth)
            hold(ax,'off'); grid(ax,'off'); box(ax,'on');
            ylim(ax, ylims_sc{m}); xlim(ax, [t(1) t(end)]);
            if si == 1, title(ax, meas_lbl{m}, 'Interpreter','latex', 'FontSize', 13); end
            if m  == 1
                ylabel(ax, sprintf('SNR = %g', snr_list(si)), ...
                       'FontWeight','bold', 'FontSize', 11);
            end
            if si < nSNR, set(ax,'XTickLabel',[]); else, xlabel(ax,'Time (s)'); end
            if si == 1 && m == 1, hB = hb; hM = hm; hR = hr; end
            ax=gca;
            ax.FontSize=15;
        end
    end
    lg = legend([hB hM hR], {'IQR (25-75%)','Median','Theoretical'}, ...
                'Orientation','horizontal');
    lg.Layout.Tile = 'north';
    title(tl, sprintf('%s — scalar measures vs SNR (median + IQR over %d realizations)', ...
          noiseName{e}, nReal), 'FontSize', 14, 'FontWeight','bold');

    fname = fullfile(outRoot, sprintf('grid_traces_%s.eps', noiseTypes{e}));
    drawnow;
    exportgraphics(fh, fname, 'Resolution', saveRes, 'ContentType','vector');   % EPS, 600 dpi
    savefig(fh, fullfile(outRoot, sprintf('grid_traces_%s.fig', noiseTypes{e}))); % MATLAB .fig
    fprintf('  saved %s\n', fname);
    if closeAfterSave, close(fh); end
end

%% ============================================================
%  (b) PER-COLOUR SPECTROGRAM GRID   (rows = SNR, cols = measure)
% ============================================================
for e = 1:nN
    fh = figure('Units','inches','Position',[0.5 0.5 13 9.5], ...
                'Name',sprintf('Spectrograms | %s', noiseName{e}), 'Theme','Light','WindowState','maximized');
    colormap(fh, cmap.VRVmap);
    tl = tiledlayout(nSNR, 4, 'TileSpacing','compact', 'Padding','compact');
    for si = 1:nSNR
        for m = 1:4
            ax = nexttile;
            data_sp = SPmed{si,e}{m}(:, t_in:end);
            imagesc(ax, t(t_in:end), f_vec/fs, data_sp, clims_sp{m});   % normalized freq (f/fs)
            set(ax,'YDir','normal'); axis(ax,'tight'); ylim(ax,[0 0.5]);% 0 .. Nyquist = 0.5
            if si == 1, title(ax, meas_lbl{m}, 'Interpreter','latex', 'FontSize', 13); end
            if m  == 1
                ylabel(ax, sprintf('SNR = %g\nNorm. freq (f/f_s)', snr_list(si)), ...
                       'FontWeight','bold', 'FontSize', 11);
            end
            if si < nSNR, set(ax,'XTickLabel',[]); else, xlabel(ax,'Time (s)'); end
            % one colorbar per column (shared clim down the column)
            if si == nSNR
                cb = colorbar(ax,'Location','southoutside');
                cb.Label.String = 'nats';
            end
            ax=gca;
            ax.FontSize=15;
        end
    end
    title(tl, sprintf('%s — median spectrograms vs SNR (over %d realizations)', ...
          noiseName{e}, nReal), 'FontSize', 14, 'FontWeight','bold');

    fname = fullfile(outRoot, sprintf('grid_spectro_%s.eps', noiseTypes{e}));
    drawnow;
    exportgraphics(fh, fname, 'Resolution', saveRes, 'ContentType','vector');    % EPS, 600 dpi
    savefig(fh, fullfile(outRoot, sprintf('grid_spectro_%s.fig', noiseTypes{e}))); % MATLAB .fig
    fprintf('  saved %s\n', fname);
    if closeAfterSave, close(fh); end
end

%% ============================================================
%  (c) Iy OVERLAY SUMMARY   (2x2, one panel per noise color)
% ============================================================
fh = figure('Units','inches','Position',[0.5 0.5 12 8.5], ...
            'Name','Iy overlay vs SNR', 'Theme','Light','WindowState','maximized');
tl = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');
hSNR = gobjects(1,nSNR); hRef = [];
for e = 1:nN
    ax = nexttile; hold(ax,'on');
    for si = 1:nSNR
        h = plot(t, TRmed{si,e}{4}, 'Color', cols_snr(si,:), 'LineWidth', 1.6);
        if e == 1, hSNR(si) = h; end
    end
    hr = plot(t, TR_th{4}, 'k-', 'LineWidth', 1.4);   % theoretical (ground truth)
    if e == 1, hRef = hr; end
    hold(ax,'off'); grid(ax,'off'); box(ax,'on');
    ylim(ax, ylims_sc{4}); xlim(ax, [t(1) t(end)]);
    xlabel(ax,'Time (s)'); ylabel(ax,'$I_{Y;X_1;X_2}$', 'Interpreter','latex', 'FontSize', 13);
    title(ax, noiseName{e}, 'FontSize', 12);
    ax=gca;
    ax.FontSize=15;
end
lg = legend([hSNR hRef], [arrayfun(@(s) sprintf('SNR = %g', s), snr_list, ...
            'UniformOutput', false), {'Theoretical'}], 'Orientation','horizontal');
lg.Layout.Tile = 'north';
title(tl, sprintf('Interaction information I_{Y;X_1;X_2} — SNR degradation (over %d realizations)', ...
      nReal), 'FontSize', 14, 'FontWeight','bold');

fname = fullfile(outRoot, 'summary_Iy_overlay.eps');
drawnow;
exportgraphics(fh, fname, 'Resolution', saveRes, 'ContentType','vector');   % EPS, 600 dpi
savefig(fh, fullfile(outRoot, 'summary_Iy_overlay.fig'));                    % MATLAB .fig
fprintf('  saved %s\n', fname);
if closeAfterSave, close(fh); end

%% ============================================================
%  (d) COUPLING PARAMETERS + NOISE SPECTRA
%      Left : the time-varying coupling a(t_n) (and 1 - a(t_n))
%      Right: power spectra of the colored noises (flat / 1/f / 1/f^2)
% ============================================================
spec_types = {'white', 'pink', 'brown'};
spec_name  = {'White (flat)', 'Pink (1/f)', 'Brown (1/f^2)'};
spec_cols  = [0.35 0.35 0.35;    % white -> grey
              0.85 0.20 0.55;    % pink
              0.55 0.30 0.10];   % brown
Npsd = 4096;     % FFT length per segment for the PSD estimate
nAvg = 200;      % segments averaged (variance reduction)

fh = figure('Units','inches','Position',[0.5 0.5 13 5], ...
            'Name','Coupling parameters & noise spectra', 'Theme','Light','WindowState','maximized');
tl = tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');

% --- left: time-varying coupling parameters ---
ax = nexttile; hold(ax,'on');
plot(ax, t, c1, 'Color', col4(1,:), 'LineWidth', 1.8);
plot(ax, t, c2, '--', 'Color', col4(3,:), 'LineWidth', 1.4);
hold(ax,'off'); grid(ax,'off'); box(ax,'on');
xlim(ax,[t(1) t(end)]); ylim(ax,[-0.15 1.15]);
xlabel(ax,'Time (s)'); ylabel(ax,'Coupling strength');
legend(ax, {'$a(t_n)$','$1-a(t_n)$'}, 'Interpreter','latex', ...
       'Location','east', 'FontSize', 12);
title(ax, sprintf('Time-varying coupling ( %.2f Hz square wave, DC = %d%% )', ...
      f_osc, DC), 'FontSize', 12);
ax=gca;
ax.FontSize=15;

% --- right: power spectra of the colored noises ---
rng(base_seed);                              % reproducible spectra
ax = nexttile; hold(ax,'on');
hS = gobjects(1, numel(spec_types));
for q = 1:numel(spec_types)
    [Pxx, fpsd] = noiseSpectrum(spec_types{q}, Npsd, nAvg, fs);
    idx = 2:numel(fpsd);                     % skip DC (zeroed for pink/brown)
    hS(q) = loglog(ax, fpsd(idx), Pxx(idx), 'Color', spec_cols(q,:), 'LineWidth', 1.8);
end
set(ax, 'XScale','log', 'YScale','log');
hold(ax,'off'); grid(ax,'off'); box(ax,'on');
xlim(ax, [fpsd(2) fs/2]);
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'PSD (a.u.)');
legend(ax, hS, spec_name, 'Location','southwest', 'FontSize', 11);
ax=gca;
ax.FontSize=15;
title(ax, sprintf('Noise power spectra (unit variance, %d-avg periodogram)', nAvg), ...
      'FontSize', 12);

title(tl, 'Coupling parameters and colored-noise spectra', ...
      'FontSize', 14, 'FontWeight','bold');

fname = fullfile(outRoot, 'summary_coupling_noise_spectra.eps');
drawnow;
exportgraphics(fh, fname, 'Resolution', saveRes, 'ContentType','vector');   % EPS, 600 dpi
savefig(fh, fullfile(outRoot, 'summary_coupling_noise_spectra.fig'));       % MATLAB .fig
fprintf('  saved %s\n', fname);
if closeAfterSave, close(fh); end

fprintf('\nTotal time %.0fs. Done — %d trace grids + %d spectro grids + 1 overlay + 1 coupling/noise figure saved under %s\n', ...
        toc(tStart), nN, nN, outRoot);

%% ============================================================
%  (e) SINGLE-SNR SPECTROGRAM STRIP ACROSS NOISE COLOURS
%      Worst-case SNR, noise colours (rows) x measures (cols).
%      Shows whether pink (1/f) and brown (1/f^2) erode the LF band
%      relative to white -- the spectral evidence the time-domain
%      overlay cannot provide.
% ============================================================
snr_pick  = 1;            % SNR (linear) to display; 1 = 0 dB (worst case)
cols_pick = [2 3 4];      % rows: noise colour idx 1=clean 2=white 3=pink 4=brown
meas_pick = 1:4;          % cols: 1=Fy1 2=Fy2 3=Fy12 4=Iy

% --- frequency axis: NORMALIZED FREQUENCY (f/fs) -------------------------
% fs is arbitrary in this simulation, so the axis is shown as normalized
% frequency f/fs (cycles/sample), running 0 .. 0.5 (Nyquist). Dividing
% f_vec by fs makes the peaks land at the pole definitions directly:
% node 1 pole 0.35 and nodes 4/6 pole 0.10. Relabels the axis only; the
% data are untouched.
FAXIS = f_vec / fs;
ymax  = 0.5;
% -------------------------------------------------------------------------

si = find(snr_list == snr_pick, 1);
if isempty(si), error('snr_pick=%g not found in snr_list', snr_pick); end

nR = numel(cols_pick);  nC = numel(meas_pick);   % rows = noise colours, cols = measures
fh = figure('Units','inches','Position',[0.5 0.5 3.2*nC 2.4*nR], ...
    'Name',sprintf('Spectro strip | SNR=%g', snr_pick), 'Theme','Light','WindowState','maximized');
colormap(fh, cmap.VRVmap);
tl = tiledlayout(nR, nC, 'TileSpacing','compact', 'Padding','compact');

for ri = 1:nR
    e = cols_pick(ri);
    for ci = 1:nC
        m  = meas_pick(ci);
        ax = nexttile;
        data_sp = SPmed{si,e}{m}(:, t_in:end);
        % common clim (clims_sp is one shared scale for all measures) -> the
        % single colorbar below is directly comparable across every panel.
        imagesc(ax, t(t_in:end), FAXIS, data_sp, clims_sp{m});
        set(ax,'YDir','normal'); axis(ax,'tight'); ylim(ax,[0 ymax]);
        if ri == 1, title(ax, meas_lbl{m}, 'Interpreter','latex', 'FontSize', 13); end
        if ci == 1
            ylabel(ax, sprintf('%s\nNorm. freq (f/f_s)', noiseName{e}), 'FontSize', 11);
        end
        if ri < nR, set(ax,'XTickLabel',[]); else, xlabel(ax,'Time (s)'); end
        ax=gca;
        ax.FontSize=15;
        
    end
end
% single shared colorbar: valid ONLY because clims_sp is one common scale for
% all measures, so measures can be compared directly. If you switch to
% per-measure clims, add one colorbar per column instead.
cb = colorbar('Location','eastoutside');  cb.Label.String = 'nats';
cb.Layout.Tile = 'east';
title(tl, sprintf('Median spectrograms at SNR = %g (%.0f dB) across noise colours', ...
    snr_pick, 10*log10(snr_pick)), 'FontSize', 13, 'FontWeight','bold');

fname = fullfile(outRoot, sprintf('strip_spectro_SNR%g.eps', snr_pick));
drawnow;
exportgraphics(fh, fname, 'Resolution', saveRes, 'ContentType','vector');   % EPS, 600 dpi
savefig(fh, fullfile(outRoot, sprintf('strip_spectro_SNR%g.fig', snr_pick))); % MATLAB .fig
fprintf('  saved %s\n', fname);

%% ============================================================
%  LOCAL FUNCTIONS
% ============================================================
function Xn = addColoredNoise(X, type, snr_lin)
% snr_lin : LINEAR SNR power ratio (signal power / noise power)
if strcmpi(type, 'clean')
    Xn = X;  return;            % baseline: no noise added
end
[Mc, T] = size(X);
Xn = zeros(Mc, T);
for ch = 1:Mc
    s   = X(ch, :);
    eta = coloredNoise(type, T);          % unit-variance colored noise
    psig = var(s, 1);
    pn   = psig / snr_lin;                 % target noise power for this SNR ratio
    Xn(ch, :) = s + sqrt(pn) * eta;
end
end


function [Pxx, fpsd] = noiseSpectrum(type, N, nAvg, fs)
% Averaged one-sided periodogram of unit-variance colored noise.
%   type : 'white' | 'pink' | 'brown'
%   N    : FFT length per segment
%   nAvg : number of independent segments averaged
%   fs   : sampling frequency (for the frequency axis)
acc = zeros(1, N);
for i = 1:nAvg
    y   = coloredNoise(type, N);          % unit-variance segment
    Y   = fft(y);
    acc = acc + (abs(Y).^2) / N;          % periodogram
end
acc  = acc / nAvg;
half = 1:floor(N/2);
Pxx  = acc(half);
fpsd = (0:floor(N/2)-1) * fs / N;
end


function y = coloredNoise(type, N)
w = randn(1, N);
switch lower(type)
    case 'white'
        y = w;
    otherwise
        switch lower(type)
            case 'pink',  alpha = 1;
            case 'brown', alpha = 2;
            otherwise, error('unknown noise type "%s"', type);
        end
        Xf = fft(w);
        kk = 0:N-1;
        fk = min(kk, N - kk);              % folded frequency (symmetric, DC=0)
        fk(1) = 1;                          % avoid 1/0 (DC zeroed next line)
        scale = fk .^ (-alpha/2);           % amplitude ~ 1/f^(alpha/2)
        scale(1) = 0;                       % remove DC
        y = real(ifft(Xf .* scale));
end
y = y - mean(y);
sd = std(y, 1);
if sd > 0, y = y / sd; end                  % unit variance
end


function [Fy1,Fy2,Fy12,Iy, fy1,fy2,fy12,iy, f] = ...
    runPipeline(X, c_ff, p, nfft, fs, nobs)
[A_i,  Su_i]     = idMVAR(squeeze(X), p, 0);
[Am_rls, Su_rls] = tvID_VAR_RLS_IC(X, p, c_ff, A_i, Su_i);

[~, ~, f] = sir_VARspectra(Am_rls(:,:,1), Su_rls(:,:,1), nfft, fs);

ws = warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

Fy1 = zeros(1,nobs);  Fy2 = zeros(1,nobs);
Fy12 = zeros(1,nobs); Iy  = zeros(1,nobs);
fy1 = zeros(nfft,nobs);  fy2 = zeros(nfft,nobs);
fy12 = zeros(nfft,nobs); iy  = zeros(nfft,nobs);

for n = 1:nobs
    S_n = sir_VARspectra(Am_rls(:,:,n), Su_rls(:,:,n), nfft, fs);
    out = iispectral(S_n, 7, [1 2 3], [4 5 6]);
    Fy1(n)  = out.Fy_x1;    Fy2(n) = out.Fy_x2;
    Fy12(n) = out.Fy_x1x2;  Iy(n)  = out.Iy_x1_x2;
    fy1(:,n)  = out.fy_x1;    fy2(:,n) = out.fy_x2;
    fy12(:,n) = out.fy_x1x2;  iy(:,n)  = out.iy_x1_x2;
end

warning(ws);
end
