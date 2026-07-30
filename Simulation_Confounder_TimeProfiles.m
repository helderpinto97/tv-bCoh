%% Simulation_Confounder_TimeProfiles.m
clear all; close all; clc;
addpath('functions');

%% ---------------- parameters ----------------
fs    = 2;      % sampling frequency (Hz)
nfft  = 1024;      % FFT resolution
nobs  = 4000;     % time samples
n_gate_cycles = 2.5;   % number of confounder on/off cycles over the whole record
DC    = 50;       % duty cycle (%)

% on/off switching frequency of the confounder, adapted to nobs and fs so the
% number of cycles across the record stays fixed regardless of nobs/fs:
%   n_gate_cycles = f_osc * (nobs/fs)  ->  f_osc = n_gate_cycles * fs / nobs
f_osc = n_gate_cycles * fs / nobs;   % (Hz)

c     = 0.9;      % genuine coupling X1->Y , X2->Y (strong -> visible synergy)
gON   = 1.0;      % common-noise level when ON (relative to each channel std)

nReal = 2;       % number of realizations
p_est = 2;        % estimation order (common noise -> needs > true VAR order)
c_ff  = 0.98;     % RLS forgetting factor
base_seed = 9000;

iY = 3; iX1 = 1; iX2 = 2;
t  = (1:nobs)/fs;

% confounder gate: 1 = ON, 0 = OFF  (square wave)
gate = double(square(2*pi*f_osc*t, DC) < 0);

fprintf('=== Confounder impact: common Gaussian noise (with vs without) ===\n');
fprintf('  nobs=%d | p_est=%d | ff=%.2f | gON=%.2f | %d realizations\n\n', ...
        nobs, p_est, c_ff, gON, nReal);

%% ---------------- base 3-node VAR (stationary, built once) ----------------
parB.poles = {[0.90 0.35], [0.90 0.35], [0.90 0.35]};  % X1, X2, Y all resonant at 35 Hz (0.35*fs)
parB.coup  = [1 3 1 c; 2 3 1 c];         % X1->Y , X2->Y
parB.Su    = [1 1 1];                   % small target self-drive -> sources dominate Y
[AmB, SuB] = theoreticalVAR(3, parB);
AmB2d = AmB';                            % 3 x 3p
pB    = size(AmB2d,2)/3;
AmT   = repmat(reshape(AmB2d,[3,3,pB]), [1 1 1 nobs]);   % constant in time
SuT   = repmat(SuB, [1 1 nobs]);

%% ---------------- theoretical (ground-truth) profiles ----------------
% clean spectrum (constant) + per-channel std from it
[Sclean,~,~] = sir_VARspectra(AmB2d, SuB, nfft, fs);
sd = zeros(3,1);
for ch = 1:3, sd(ch) = sqrt(mean(real(squeeze(Sclean(ch,ch,:))))); end

oWO  = iispectral(Sclean, iY, iX1, iX2);
THWO = repmat([real(oWO.Fy_x1);real(oWO.Fy_x2);real(oWO.Fy_x1x2);real(oWO.Iy_x1_x2)], 1, nobs);

THW = zeros(4,nobs);                                  % common white noise adds
RR  = sd*sd';                                         % rank-1 flat term (sd_i*sd_j)
for n = 1:nobs
    Sn = Sclean + (gON*gate(n))^2 * RR;               % add to every frequency
    o  = iispectral(Sn, iY, iX1, iX2);
    THW(:,n) = [real(o.Fy_x1);real(o.Fy_x2);real(o.Fy_x1x2);real(o.Iy_x1_x2)];
end

%% ---------------- estimation over realizations (paired seeds) ----------------
trW  = zeros(nReal,4,nobs);        % WITH common noise  (scalar measures)
trWO = zeros(nReal,4,nobs);        % WITHOUT (clean)
SPW  = zeros(nReal,4,nfft,nobs);   % WITH  (spectral measures)
SPWO = zeros(nReal,4,nfft,nobs);   % WITHOUT
ws = warning('off','MATLAB:singularMatrix'); warning('off','MATLAB:nearlySingularMatrix');
for r = 1:nReal
    [trWO(r,:,:), spwo, ~] = runReal(AmT, SuT, 0,   gate, p_est, c_ff, nfft, fs, nobs, base_seed+r, iY,iX1,iX2);
    [trW(r,:,:),  spw,  fvec] = runReal(AmT, SuT, gON, gate, p_est, c_ff, nfft, fs, nobs, base_seed+r, iY,iX1,iX2);
    SPWO(r,:,:,:) = spwo;  SPW(r,:,:,:) = spw;
    fprintf('  realization %2d/%d\n', r, nReal);
end
warning(ws);

[medW, q1W, q3W]  = reduceTrace(trW);
[medWO,q1WO,q3WO] = reduceTrace(trWO);
medSPW  = squeeze(median(SPW,  1, 'omitnan'));   % 4 x nfft x nobs
medSPWO = squeeze(median(SPWO, 1, 'omitnan'));

%% ---------------- plot (2x2, one measure per panel) ----------------
colWO = [0.00 0.45 0.85];  tintWO = 1 - 0.55*(1-colWO);   % no confounder (blue)
colW  = [0.85 0.10 0.10];  tintW  = 1 - 0.55*(1-colW);    % with confounder (red)
lbl   = {'$F_{Y;X_1}$','$F_{Y;X_2}$','$F_{Y;X_1,X_2}$','$I_{Y;X_1;X_2}$'};

fh = figure('Units','inches','Position',[0.5 0.5 12 7.5], ...
            'Name','Confounder impact','Theme','Light','WindowState','maximized');
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
H = gobjects(1,5);
% uniform y-limits across all four panels
yall = [];
for mm = 1:4
    [a0,b0] = bandfix(medWO(mm,:),q1WO(mm,:),q3WO(mm,:));
    [a1,b1] = bandfix(medW(mm,:), q1W(mm,:), q3W(mm,:));
    yall = [yall a0 b0 a1 b1 THWO(mm,:)]; %#ok<AGROW>
end
yl = [min(yall)-0.1, max(yall)+0.1];
for m = 1:4
    ax = nexttile; hold(ax,'on');
    [aWO,bWO] = bandfix(medWO(m,:),q1WO(m,:),q3WO(m,:));
    [aW ,bW ] = bandfix(medW(m,:), q1W(m,:), q3W(m,:));
    shadeON(ax, t, gate, yl);                                  % grey = confounder active (shared yl)
    hbWO = fill(ax,[t fliplr(t)],[bWO fliplr(aWO)],tintWO,'FaceAlpha',0.5,'EdgeColor','none');
    hbW  = fill(ax,[t fliplr(t)],[bW  fliplr(aW )],tintW, 'FaceAlpha',0.5,'EdgeColor','none');
    hmWO = plot(ax,t,medWO(m,:),'Color',colWO,'LineWidth',1.8);
    hmW  = plot(ax,t,medW(m,:), 'Color',colW, 'LineWidth',1.8);
    ht   = plot(ax,t,THWO(m,:),'-','Color','k','LineWidth',1.3);   % theoretical (clean)
    % yline(ax,0,':','Color',[0.4 0.4 0.4]);
    hold(ax,'off'); grid(ax,'off'); box(ax,'on');
    xlim(ax,[t(1) t(end)]); ylim(ax,yl);
    xlabel(ax,'Time (s)'); ylabel(ax,'nats');
    title(ax, lbl{m}, 'Interpreter','latex', 'FontSize', 13);
    if m==1, H = [hbWO hmWO hbW hmW ht]; end
    axis square;
    ax=gca;
    ax.FontSize=16;
    ax.LineWidth=2;
end
lg = legend(H, {'No confounder — IQR','No confounder — median', ...
                'With confounder — IQR','With confounder — median','Theoretical'}, ...
            'Orientation','horizontal','NumColumns',5);
lg.Layout.Tile = 'north';
outRoot = fullfile(pwd,'Confounder_Analysis_Figures');
if ~exist(outRoot,'dir'), mkdir(outRoot); end
fname = fullfile(outRoot,'confounder_impact_common_noise.png');
% drawnow; exportgraphics(fh, fname, 'Resolution', 300);
fprintf('\nsaved %s\n', fname);

%% ---------------- time-frequency measures (spectrogram) ----------------
% rows: no confounder / with confounder ; cols: the four spectral measures
cmapTF = parula(256);
try 
    cm = load('colmap.mat','VRVmap'); cmapTF = cm.VRVmap; 
end 

% uniform symmetric color limit across all panels
d = [reshape(medSPWO,1,[]), reshape(medSPW,1,[])];
d = d(isfinite(d));
clm = prctile(abs(d),99);
if isempty(clm) || clm==0, clm = 1; end
dd = diff([0 gate 0]);                              % confounder on/off edges
onSt = t(dd==1); onEn = t(find(dd==-1)-1);

fh2 = figure('Units','inches','Position',[0.5 0.5 13 6.5], ...
             'Name','Confounder impact (time-frequency)','Theme','Light','WindowState','maximized');
colormap(fh2, cmapTF);
tl2 = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
rowsSP = {medSPWO, medSPW}; rowname = {'No confounder','With confounder'};
for ri = 1:2
    for m = 1:4
        ax = nexttile;
        data = squeeze(rowsSP{ri}(m,:,:));          % nfft x nobs
        imagesc(ax, t, fvec, data, [-clm clm]);
        set(ax,'YDir','normal'); ylim(ax,[0 fs/2]);
        hold(ax,'on');
        for q = 1:numel(onSt)
            xline(ax, onSt(q), 'w--', 'LineWidth',0.8);
            xline(ax, onEn(q), 'w--', 'LineWidth',0.8);
        end
        hold(ax,'off');
        if ri==1, title(ax, lbl{m}, 'Interpreter','latex', 'FontSize',13); end
        if m==1,  ylabel(ax, sprintf('%s\nFreq (Hz)', rowname{ri}), 'FontWeight','bold'); end
        if ri==2, xlabel(ax,'Time (s)'); else, set(ax,'XTickLabel',[]); end
        cb = colorbar(ax); cb.Label.String = 'nats';
        axis square;
        ax=gca;
        ax.FontSize=16;
        ax.LineWidth=2;
    end
end
% fname2 = fullfile(outRoot,'confounder_impact_timefreq.png');
% drawnow; exportgraphics(fh2, fname2, 'Resolution', 300);
% fprintf('saved %s\n', fname2);

%% ================= local functions =================
function [tr, sp, f] = runReal(AmT, SuT, gconf, gate, p_est, c_ff, nfft, fs, nobs, seed, iY,iX1,iX2)
% one realization: generate clean data, optionally add a COMMON gated
% Gaussian noise to all channels, then estimate the measure profiles.
%   tr : 4 x nobs        scalar (time-domain) measures
%   sp : 4 x nfft x nobs spectral (time-frequency) measures
%   f  : 1 x nfft        frequency axis
rng(seed);
Xclean = var_to_tsdata_nonstat(AmT, SuT, 1);   % 3 x nobs
w      = randn(1, nobs);                        % the common confounder noise
sdc    = std(Xclean, 0, 2);                     % per-channel std (3 x 1)
Xobs   = Xclean + (gconf*gate) .* (sdc .* w);   % add same w to all 3 channels (gated)
X      = zscore(permute(Xobs,[1 3 2]), 0, 3);   % 3 x 1 x nobs
[A_i, Su_i]      = idMVAR(squeeze(X), p_est, 0);
[Am_rls, Su_rls] = tvID_VAR_RLS_IC(X, p_est, c_ff, A_i, Su_i);
tr = zeros(4,nobs);
sp = zeros(4,nfft,nobs);
for n = 1:nobs
    [S,~,f] = sir_VARspectra(Am_rls(:,:,n), Su_rls(:,:,n), nfft, fs);
    o = iispectral(S, iY, iX1, iX2);
    tr(:,n)   = [real(o.Fy_x1); real(o.Fy_x2); real(o.Fy_x1x2); real(o.Iy_x1_x2)];
    sp(1,:,n) = real(o.fy_x1);   sp(2,:,n) = real(o.fy_x2);
    sp(3,:,n) = real(o.fy_x1x2); sp(4,:,n) = real(o.iy_x1_x2);
end
sp(~isfinite(sp)) = NaN;
end

function [med,q1,q3] = reduceTrace(tr)
% tr: nReal x 4 x nobs  ->  median / quartiles : 4 x nobs
nobs = size(tr,3); med = zeros(4,nobs); q1 = med; q3 = med;
for m = 1:4
    v = squeeze(tr(:,m,:)); v(~isfinite(v)) = NaN;   % nReal x nobs
    med(m,:) = median(v,1,'omitnan');
    q1(m,:)  = prctile(v,25,1);
    q3(m,:)  = prctile(v,75,1);
end
end

function [a,b] = bandfix(med,q1,q3)
% replace non-finite band edges so fill() is well-defined
a = q1; b = q3; bad = ~isfinite(a)|~isfinite(b);
a(bad) = med(bad); b(bad) = med(bad);
a(~isfinite(a)) = 0; b(~isfinite(b)) = 0;
end

function shadeON(ax, t, gate, yl)
% Shade the time intervals where gate==1 (confounder active).
d  = diff([0 gate 0]);
st = find(d==1);  en = find(d==-1)-1;
for k = 1:numel(st)
    patch(ax, t([st(k) en(k) en(k) st(k)]), [yl(1) yl(1) yl(2) yl(2)], ...
          [0.88 0.88 0.88], 'EdgeColor','none', 'HandleVisibility','off');
end
end
