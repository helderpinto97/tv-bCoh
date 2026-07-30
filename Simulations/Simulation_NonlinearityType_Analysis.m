%% Simulation_NonlinearityType_Analysis_Realizations.m
clear; close all; clc;
addpath('../functions');

%% ---------------- parameters ----------------
fs    = 2;     nfft = 1024;    nobs = 4000;
n_osc_cycles = 2.5;                         % number of coupling on/off cycles over the record
DC    = 50;      c1max = 0.8;   M = 7;          c_ff = 0.98;
p     = 2;                                  % true model order

% coupling oscillation frequency, adapted to nobs and fs so the number of
% cycles across the record stays fixed regardless of nobs/fs:
%   n_osc_cycles = f_osc * (nobs/fs)  ->  f_osc = n_osc_cycles * fs / nobs
f_osc = n_osc_cycles * fs / nobs;           % (Hz)
nReal     = 100;                             % Monte-Carlo realizations per type
base_seed = 4000;                           % per-realization seed: base_seed + r
                                            % (shared across types -> paired comparison)

t = (1:nobs) / fs;
c1_true = -c1max * square(2*pi*f_osc*t, DC);
c1_true(c1_true == -c1max) = 0;
f_LF_nom = 0.10;  f_HF_nom = 0.35;

meas_lbl = {'$F_{Y;X_1}$','$F_{Y;X_2}$','$F_{Y;X_1,X_2}$','$I_{Y;X_1;X_2}$'};

%% ---------------- define & normalize the nonlinearity shapes ----------------
% Raw shapes phi(u). They are standardized below so that, for U~N(0,1),
% E[phi_hat] = 0 and Var[phi_hat] = 1.  Identity -> exactly linear coupling.
raw = { @(u) u, ...                 % 1 Linear            
        @(u) tanh(1.5*u), ...       % 2 Saturating tanh   
        @(u) u.^3, ...              % 3 Cubic             
        @(u) u.^2, ...              % 4 Quadratic         
        @(u) abs(u), ...            % 5 Absolute value    
        @(u) max(0,u) };            % 6 ReLU / half-wave  
type_name = {'Linear','Saturating (tanh)','Cubic','Quadratic','Abs |x|','ReLU'};
nT = numel(raw);

% Standardize each raw shape so that, for U~N(0,1), E[phihat]=0 and Var[phihat]=1.
rng(7);  Zref = randn(1, 20000);      % fixed reference sample for normalization
phihat = cell(1, nT);  rho = zeros(1, nT);
for kk = 1:nT
    v  = raw{kk}(Zref);
    m  = mean(v);   sd = std(v);
    phihat{kk} = @(u) (raw{kk}(u) - m) / sd;         % zero-mean, unit-var
    rho(kk)    = corr(phihat{kk}(Zref).', Zref.');   % linear projection (column vectors)
end


%% ---------------- reference per-process std (from the linear run) ----------
% Used to feed each shape a standardized source, so g operates in the same
% regime regardless of type.
rng(42);
Xlin = genNL(M, c1_true, f_LF_nom, f_HF_nom, nobs, @(u) u, ones(1,M));
sref = std(Xlin, 0, 2)';            % 1 x M

%% ============================================================
%  THEORETICAL measures from the TRUE (linear) TV-VAR  (ground truth)
%  Estimation-free closed form -> compute once.
% ============================================================
par.poles{1} = [0.9 f_HF_nom];      % HF autonomous oscillator (process 1)
par.poles{4} = [0.9 f_LF_nom];      % LF autonomous oscillator (process 4)
par.poles{6} = [0.9 f_LF_nom];      % LF autonomous oscillator (process 6)
par.poles{7} = [];                  % driven process, no self-pole

AmT = []; SuT = [];
for n = 1:nobs
    C = c1_true(n);
    par.coup = [1 2 1 0.8;          % 1 -> 2  lag 1   0.8
                1 3 1 C;            % 1 -> 3  lag 1   a(t)
                4 5 2 C;            % 4 -> 5  lag 2   a(t)
                6 5 1 0.8;          % 6 -> 5  lag 1   0.8
                5 7 2 C;            % 5 -> 7  lag 2   a(t)
                2 7 2 0.8];         % 2 -> 7  lag 2   0.8
    par.Su = ones(1, M);
    [Am1, Su1] = theoreticalVAR(M, par);
    Am1 = Am1';
    Am1 = reshape(Am1, [M, M, size(Am1,2)/M]);
    AmT(:,:,:,n) = Am1;             
    SuT(:,:,n)   = Su1;             
end
p_th = size(AmT, 3);                % = p
for n = 1:nobs
    SuT(:,:,n) = (SuT(:,:,n) + SuT(:,:,n)')/2 + 1e-12*eye(M);
end

Fy1_th  = zeros(1,nobs);  Fy2_th = zeros(1,nobs);
Fy12_th = zeros(1,nobs);  Iy_th  = zeros(1,nobs);
fy2_th  = zeros(nfft,nobs);  iy_th = zeros(nfft,nobs);
wsth = warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
for n = 1:nobs
    Am_n = reshape(AmT(:,:,:,n), M, M*p_th);   % [A(1) ... A(p)]
    S_n  = sir_VARspectra(Am_n, SuT(:,:,n), nfft, fs);
    out  = iispectral(S_n, 7, [1 2 3], [4 5 6]);
    Fy1_th(n)  = out.Fy_x1;    Fy2_th(n)  = out.Fy_x2;
    Fy12_th(n) = out.Fy_x1x2;  Iy_th(n)   = out.Iy_x1_x2;
    fy2_th(:,n) = out.fy_x2;   iy_th(:,n) = out.iy_x1_x2;
end
warning(wsth);
TR_th = {Fy1_th, Fy2_th, Fy12_th, Iy_th};
for m = 1:4                                   % sanitize
    v = real(TR_th{m});  v(~isfinite(v)) = NaN;  TR_th{m} = v;
end

%% ============================================================
%  COMPUTE — all nonlinearity types x realizations
% ============================================================
TRmed = cell(1,nT); TRq1 = cell(1,nT); TRq3 = cell(1,nT);   % {type}{measure} 1 x nobs
SPmed = cell(1,nT);                                          % {type}{measure} nfft x nobs
R2med = nan(4, nT);  SNRmed = nan(4, nT);                    % median quality vs c1_true
f_vec = [];

fprintf('=== Simulation_NonlinearityType_Analysis_Realizations.m ===\n');
fprintf('  fs=%d Hz | nobs=%d | M=%d | ff=%.2f | %d realizations/type\n\n', ...
        fs, nobs, M, c_ff, nReal);
tStart = tic;
for kk = 1:nT
    trace = {zeros(nReal,nobs), zeros(nReal,nobs), ...
             zeros(nReal,nobs), zeros(nReal,nobs)};          % {Fy1,Fy2,Fy12,Iy}
    SP    = {zeros(nReal,nfft,nobs), zeros(nReal,nfft,nobs), ...
             zeros(nReal,nfft,nobs), zeros(nReal,nfft,nobs)};% {fy1,fy2,fy12,iy}
    R2r   = nan(4,nReal);  SNRr = nan(4,nReal);

    for r = 1:nReal
        rng(base_seed + r);             % same innovations across types (paired)
        X  = genNL(M, c1_true, f_LF_nom, f_HF_nom, nobs, phihat{kk}, sref);
        Xz = permute(X, [1 3 2]);
        Xz = zscore(Xz, 0, 3);

        [F1,F2,F12,Iy, s1,s2,s12,siy, f_vec] = ...
            runScalarPipeline(Xz, c_ff, M, p, nfft, fs, nobs);

        Fm = {F1,F2,F12,Iy};  sp = {s1,s2,s12,siy};
        for m = 1:4
            trace{m}(r,:) = Fm{m};
            SP{m}(r,:,:)  = sp{m};
            [R2r(m,r), SNRr(m,r)] = snrMetric(Fm{m}, c1_true);
        end
    end

    % --- sanitize: near-singular spectra can make values complex/Inf ---
    for m = 1:4
        v = real(trace{m});  v(~isfinite(v)) = NaN;  trace{m} = v;
        s = real(SP{m});     s(~isfinite(s)) = NaN;  SP{m}    = s;
    end

    % --- reduce over realizations ---
    TRmed{kk} = cell(1,4); TRq1{kk} = cell(1,4); TRq3{kk} = cell(1,4); SPmed{kk} = cell(1,4);
    for m = 1:4
        TRmed{kk}{m} = median(trace{m}, 1, 'omitnan');        % 1 x nobs
        TRq1{kk}{m}  = prctile(trace{m}, 25, 1);
        TRq3{kk}{m}  = prctile(trace{m}, 75, 1);
        SPmed{kk}{m} = squeeze(median(SP{m}, 1, 'omitnan'));  % nfft x nobs
        R2med(m,kk)  = median(R2r(m,:),  'omitnan');
        SNRmed(m,kk) = median(SNRr(m,:), 'omitnan');
    end

    fprintf('  %-18s  rho=%+.2f | median R^2: Fy1=%.2f Fy2=%.2f Fy12=%.2f Iy=%.2f  [%.0fs]\n', ...
            type_name{kk}, rho(kk), R2med(1,kk), R2med(2,kk), R2med(3,kk), R2med(4,kk), toc(tStart));
end

%% ============================================================
%  PLOT — 6 x 4 grid  (rows = nonlinearity type, cols = measure)
%  Coloured median + IQR band per panel. The BLACK square wave is the
%  coupling GATE (ON/OFF).
% ============================================================
outRoot = fullfile(pwd, 'NonlinearityType_Figures_NEW');
if ~exist(outRoot,'dir'), mkdir(outRoot); end

% per-row (per nonlinearity) colours
rowCol = [0.00 0.00 0.85;    % Identity/Linear   dark blue
          0.15 0.60 0.90;    % Tanh              light blue
          0.10 0.70 0.40;    % Cubic             green
          0.85 0.10 0.10;    % Quadratic         red
          0.95 0.55 0.10;    % Abs               orange
          0.55 0.25 0.75];   % ReLu              purple
rowLbl = {'Identity','Tanh','Cubic','Quadratic','Abs','ReLu'};

% per-column (per measure) titles, y-limits and y-ticks
colLbl  = {'$F_{X_1;Y}(t_n)$','$F_{X_2;Y}(t_n)$', ...
           '$F_{X_1,X_2;Y}(t_n)$','$I_{X_1;X_2;Y}(t_n)$'};
ylimC   = {[0 3],[0 3],[0 3],[-1 1]};
ytickC  = {[0 3],[0 3],[0 3],[-1 0 1]};

gate    = double(c1_true > 0.5*c1max);     % 1 = coupling ON, 0 = OFF (used for TF transitions)

fh = figure('Units','inches','Position',[0.4 0.4 14 9], ...
            'Name','Nonlinearities — time domain', 'Theme','Light', ...
            'WindowState','maximized');
tl = tiledlayout(nT, 4, 'TileSpacing','compact', 'Padding','compact');

for kk = 1:nT                               % rows = nonlinearity type
    for m = 1:4                             % cols = measure
        ax = nexttile; hold(ax,'on');
        yl = ylimC{m};

        % theoretical value of this measure for the TRUE (linear) TV-VAR,
        % drawn in black as the ground-truth reference (replaces the old
        % coupling-parameter square wave).
        plot(ax, t, TR_th{m}, 'k-', 'LineWidth', 1.4);

        % median + IQR band
        [a,b] = bandfix(TRmed{kk}{m}, TRq1{kk}{m}, TRq3{kk}{m});
        tint  = 1 - 0.45*(1-rowCol(kk,:));
        fill(ax, [t fliplr(t)], [b fliplr(a)], tint, 'FaceAlpha',0.5, 'EdgeColor','none');
        plot(ax, t, TRmed{kk}{m}, 'Color', rowCol(kk,:), 'LineWidth', 1.8);

        hold(ax,'off'); box(ax,'on'); grid(ax,'off');
        xlim(ax,[t(1) t(end)]); ylim(ax, yl);
        set(ax,'YTick',ytickC{m}, 'FontSize',13, 'LineWidth',1.4, 'Layer','top');

        if kk == 1                                    % column titles on top row
            title(ax, colLbl{m}, 'Interpreter','latex', 'FontSize',15);
        end
        if m == 1                                     % bold row label on left column
            ylabel(ax, rowLbl{kk}, 'FontWeight','bold', 'FontSize',14);
        end
        if m ~= 1 && m ~= 4, set(ax,'YTickLabel',[]); end   % keep y-ticks only on col 1 & 4
        if kk == nT
            xlabel(ax,'Time [s]', 'FontSize',13);
        else
            set(ax,'XTickLabel',[]);
        end
    end
end

fname = fullfile(outRoot, 'Nonlinearities_TimeDomain.eps');
drawnow; exportgraphics(fh, fname, 'Resolution', 600);
savefig(fh, fullfile(outRoot, 'Nonlinearities_TimeDomain.fig'));   % editable MATLAB figure
fprintf('  figure ready: %s\n', fname);

%% ============================================================
%  PLOT — 6 x 4 grid, TIME-FREQUENCY measures (spectrograms)
%  rows = nonlinearity type, cols = measure.  The coupling GATE is marked
%  by BLACK dashed vertical lines at every ON/OFF transition.
% ============================================================
cmapTF = parula(256);
try, cm = load('colmap.mat','VRVmap'); cmapTF = cm.VRVmap;
catch, try, cm = load('../colmap.mat','VRVmap'); cmapTF = cm.VRVmap; end
end

% gate transition instants (ON/OFF edges)
dd   = diff([0 gate 0]);
onSt = t(dd == 1);  onEn = t(find(dd == -1) - 1);
trEdges = sort([onSt onEn]);

% per-column symmetric colour limits (shared down each column, over all types)
climC = cell(1,4);
for m = 1:4
    d = [];
    for kk = 1:nT, d = [d; SPmed{kk}{m}(:)]; end %#ok<AGROW>
    d  = d(isfinite(d));
    cl = prctile(abs(d), 99);
    if isempty(cl) || cl == 0, cl = 1; end
    climC{m} = [-cl cl];
end

fh2 = figure('Units','inches','Position',[0.4 0.4 14 9], ...
             'Name','Nonlinearities — time-frequency', 'Theme','Light', ...
             'WindowState','maximized');
colormap(fh2, cmapTF);
tl2 = tiledlayout(nT, 4, 'TileSpacing','compact', 'Padding','compact');

for kk = 1:nT                               % rows = nonlinearity type
    for m = 1:4                             % cols = measure
        ax = nexttile;
        imagesc(ax, t, f_vec, SPmed{kk}{m}, climC{m});
        set(ax,'YDir','normal'); ylim(ax,[0 fs/2]);
        hold(ax,'on');
        for q = 1:numel(trEdges)
            xline(ax, trEdges(q), 'k--', 'LineWidth', 1.0);
        end
        hold(ax,'off');
        xlim(ax,[t(1) t(end)]);
        set(ax,'FontSize',13,'LineWidth',1.4,'Layer','top');

        if kk == 1                                    % column titles on top row
            title(ax, colLbl{m}, 'Interpreter','latex', 'FontSize',15);
        end
        if m == 1                                     % bold row label on left column
            ylabel(ax, sprintf('%s\nFreq [Hz]', rowLbl{kk}), 'FontWeight','bold', 'FontSize',12);
        else
            set(ax,'YTickLabel',[]);
        end
        if kk == nT
            xlabel(ax,'Time [s]', 'FontSize',13);
            cb = colorbar(ax,'Location','southoutside'); cb.Label.String = 'nats';
        else
            set(ax,'XTickLabel',[]);
        end
    end
end

fname2 = fullfile(outRoot, 'Nonlinearities_TimeFreq.eps');
drawnow; exportgraphics(fh2, strrep(fname2,'.eps','.png'), 'Resolution', 600,'ContentType','vector');
% exportgraphics(fh2, fname2, 'ContentType','vector','Resolution',600);
savefig(fh2, fullfile(outRoot, 'Nonlinearities_TimeFreq.fig'));   % editable MATLAB figure
fprintf('  figure ready: %s\n', fname2);

%% ============================================================
%  LOCAL FUNCTIONS
% ============================================================

function [a,b] = bandfix(med,q1,q3)
% replace non-finite band edges so fill() is well-defined
a = q1; b = q3; bad = ~isfinite(a)|~isfinite(b);
a(bad) = med(bad); b(bad) = med(bad);
a(~isfinite(a)) = 0; b(~isfinite(b)) = 0;
end

function X = genNL(M, c1_true, f_LF, f_HF, nobs, phihat, sref)
% 7-node network of Simulation.m with a configurable coupling nonlinearity.
%   Each coupling x_j += coef * g_i(x_i),  where for source process i
%   g_i(x) = sref(i) * phihat( x / sref(i) )  (shape applied in standardised
%   units, rescaled to the source's natural scale). phihat = @(u) u reproduces
%   the original LINEAR network exactly.
r   = 0.9;
aHF = [2*r*cos(2*pi*f_HF), -r^2];
aLF = [2*r*cos(2*pi*f_LF), -r^2];

gi = @(x, idx) sref(idx) * phihat( x / sref(idx) );   % per-source nonlinearity

U = randn(M, nobs);
X = zeros(M, nobs);
pp = 2;
for n = 1:nobs
    C = c1_true(n);
    if n <= pp
        X(:,n) = U(:,n);
        continue;
    end
    % autonomous oscillators
    X(1,n) = aHF(1)*X(1,n-1) + aHF(2)*X(1,n-2) + U(1,n);
    X(4,n) = aLF(1)*X(4,n-1) + aLF(2)*X(4,n-2) + U(4,n);
    X(6,n) = aLF(1)*X(6,n-1) + aLF(2)*X(6,n-2) + U(6,n);
    % driven processes (couplings passed through the nonlinearity)
    X(2,n) = 0.8*gi(X(1,n-1),1)                    + U(2,n);
    X(3,n) = C  *gi(X(1,n-1),1)                    + U(3,n);
    X(5,n) = C  *gi(X(4,n-2),4) + 0.8*gi(X(6,n-1),6) + U(5,n);
    X(7,n) = C  *gi(X(5,n-2),5) + 0.8*gi(X(2,n-2),2) + U(7,n);
end
end


function [Fy1,Fy2,Fy12,Iy, fy1,fy2,fy12,iy, f] = ...
    runScalarPipeline(X, c_ff, ~, p, nfft, fs, nobs)
% Linear TV-VAR identification of order p, then IIR decomposition.
% Returns both the scalar (frequency-integrated) measures and the full
% time-frequency spectral measures (nfft x nobs).
[A_i, Su_i]      = idMVAR(squeeze(X), p, 0);
[Am_rls, Su_rls] = tvID_VAR_RLS_IC(X, p, c_ff, A_i, Su_i);
[~,~,f] = sir_VARspectra(Am_rls(:,:,1), Su_rls(:,:,1), nfft, fs);
ws = warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
Fy1 = zeros(1,nobs); Fy2 = zeros(1,nobs);
Fy12 = zeros(1,nobs); Iy = zeros(1,nobs);
fy1 = zeros(nfft,nobs); fy2 = zeros(nfft,nobs);
fy12 = zeros(nfft,nobs); iy = zeros(nfft,nobs);
for n = 1:nobs
    S_n = sir_VARspectra(Am_rls(:,:,n), Su_rls(:,:,n), nfft, fs);
    out = iispectral(S_n, 7, [1 2 3], [4 5 6]);
    Fy1(n)  = out.Fy_x1;   Fy2(n) = out.Fy_x2;
    Fy12(n) = out.Fy_x1x2; Iy(n)  = out.Iy_x1_x2;
    fy1(:,n)  = out.fy_x1;   fy2(:,n) = out.fy_x2;
    fy12(:,n) = out.fy_x1x2; iy(:,n)  = out.iy_x1_x2;
end
warning(ws);
end


function [r2, snr_est_dB] = snrMetric(meas_est, c1_true)
% R^2 and steady-state d-prime SNR of one scalar measure vs the clean coupling.
c_full = c1_true(:);  y_full = meas_est(:);  N = numel(c_full);
n_settle = 75;
phase_starts = [1; find(diff(c_full) ~= 0) + 1];
transient = false(N,1);
for ii = 1:numel(phase_starts)
    i1 = phase_starts(ii);  i2 = min(i1 + n_settle - 1, N);
    transient(i1:i2) = true;
end
valid = ~isnan(y_full);
c_r2 = c_full(valid);  y_r2 = y_full(valid);
if numel(y_r2) < 10, r2 = NaN; snr_est_dB = NaN; return; end
coef   = [c_r2, ones(numel(c_r2),1)] \ y_r2;
y_fit  = coef(1)*c_r2 + coef(2);
ss_res = sum((y_r2 - y_fit).^2);
ss_tot = sum((y_r2 - mean(y_r2)).^2);
r2     = max(0, 1 - ss_res / max(ss_tot, 1e-12));
steady = valid & ~transient;
c_ss = c_full(steady);  y_ss = y_full(steady);
y_on = y_ss(c_ss > 0.5);  y_off = y_ss(c_ss < 0.1);
if numel(y_on) < 2 || numel(y_off) < 2, snr_est_dB = NaN; return; end
separation   = abs(mean(y_on) - mean(y_off));
sigma_pooled = sqrt((var(y_on) + var(y_off)) / 2 + 1e-12);
snr_est_dB   = 20 * log10(max(separation, 1e-12) / sigma_pooled);
end
