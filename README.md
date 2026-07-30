# Time-Varying Block Coherence (`tv-bCoh`)

A MATLAB toolbox for estimating **time-varying block coherence** and related
linear interdependence measures in the **time–frequency domain**. The method
extends Geweke's spectral measures of linear dependence to groups (blocks) of
signals and tracks them over time through **time-varying vector autoregressive
(TV-VAR)** modeling estimated with recursive least squares (RLS).

It targets multichannel biomedical signals whose dependency structure changes
across time and frequency — e.g. cerebrovascular–autonomic interactions.

## What it computes

Given a target block `Y` and two source blocks `X1`, `X2`, the core routine
(`iispectral.m`) returns, at every time point and frequency:

- `f_{Y;X1}`, `f_{Y;X2}` — spectral linear association between `Y` and each source block
- `f_{Y;X1,X2}` — joint spectral association between `Y` and both source blocks
- `i_{Y;X1;X2} = f_{Y;X1} + f_{Y;X2} − f_{Y;X1,X2}` — the interaction (source-interaction) measure

The interaction term separates **redundant** from **synergistic** contributions
of the two source blocks to the target (positive net → redundancy,
negative → synergy, following the interaction-information convention). Each
spectral quantity is also integrated over frequency to give its time-domain
counterpart (`F_{...}`, `I_{...}`).

All block coherences are formed from determinant ratios of the (sub)blocks of
the TV-VAR cross-spectral matrix; the measures are therefore **linear /
second-order** by construction and assume the dependence of interest is captured
by the estimated VAR model.

## Repository structure

```
tv-bCoh/
├── functions/          Core and auxiliary MATLAB functions
├── Simulations/        Scripts reproducing the simulation study
├── colmap.mat          Colormap used by the plotting scripts
└── README.md
```

### `functions/` — core pipeline

These implement the estimation → spectra → measures pipeline used by the
simulation scripts:

| File | Purpose |
|------|---------|
| `theoreticalVAR.m` | Build theoretical VAR coefficients from poles/couplings (ground-truth generator) |
| `var_to_tsdata_nonstat.m` | Generate non-stationary (time-varying) Gaussian VAR realizations |
| `idMVAR.m` | Stationary MVAR identification (least squares) — used to initialize the RLS estimator |
| `tvID_VAR_RLS_IC.m` | TV-VAR identification via RLS **with initial conditions** (the estimator used in the demos) |
| `sir_VARspectra.m` | Parametric spectral matrix / transfer function from VAR coefficients |
| `iispectral.m` | Block coherences and the association / interaction measures (the method's core) |
| `subplot_tight.m` | Tight-margin subplot helper for the figures |

### `functions/` — auxiliary

Provided for completeness; **not called** by any script in `Simulations/`:

| File | Purpose | Note |
|------|---------|------|
| `tvID_VAR_RLS.m` | RLS TV-VAR identification without initial conditions | — |
| `mos_id_tv_VAR_RLS.m` | RLS-based model-order selection | requires external `knee_pt` (not bundled) |
| `lrp_LinReg.m` | Linear regression of processes from VAR parameters | — |
| `lrp_Yule.m` | Yule–Walker / Lyapunov solution for VAR covariances | — |
| `lrp_oir_mir.m` | MIR / OIR between two blocks | requires `lrp_MIR` and `oir_subindexes` (not bundled) |

> **Note:** `mos_id_tv_VAR_RLS.m` and `lrp_oir_mir.m` depend on functions that
> are not included in this repository (`knee_pt.m`, `lrp_MIR.m`,
> `oir_subindexes.m`). Add these before calling those routines. The main
> simulation pipeline does not use them.

### `Simulations/`

| Script | What it does |
|--------|--------------|
| `Simulation.m` | Minimal end-to-end demo: 7-variate TV-VAR → RLS estimation → time-varying spectra → block measures for the full network and single blocks |
| `Simulation_Confounder.m` | Effect of a common (confounding) Gaussian noise on the interaction measure, with vs. without the confounder, over Monte-Carlo realizations |
| `Simulation_Noise_Analysis_SNR.m` | Robustness of the measures across SNR levels; saves figures to an output folder |
| `Simulation_NonlinearityType_Analysis.m` | Behavior of the (linear) measures when the true coupling passes through different nonlinearities (tanh, cubic, quadratic, abs, ReLU) |
| `Noise_Plot.m` | PSD of white / pink / brown noise (illustrative; needs the DSP System Toolbox) |

## Requirements

- MATLAB (no specific version pinned)
- **Signal Processing Toolbox** — `square`
- **Statistics and Machine Learning Toolbox** — `zscore`, `corr`
- **DSP System Toolbox** — only for `Simulations/Noise_Plot.m` (`dsp.ColoredNoise`)

## Quick start

From the `Simulations/` folder (scripts add `../functions` to the path):

```matlab
cd Simulations
Simulation      % runs the full pipeline and plots the time–frequency measures
```

The typical workflow inside a script is:

```matlab
% 1) TV-VAR ground truth and realization
[Am,Su]  = theoreticalVAR(M, par);
Xgen     = var_to_tsdata_nonstat(AmT, SuT, 1);

% 2) TV-VAR estimation (RLS, initialized with a stationary fit)
[A_i,Su_i]        = idMVAR(squeeze(X), p, 0);
[Am_rls,Su_rls]   = tvID_VAR_RLS_IC(X, p, c_forget_factor, A_i, Su_i);

% 3) Time-varying spectra
[S,~,f] = sir_VARspectra(Am_rls(:,:,n), Su_rls(:,:,n), nfft, fs);

% 4) Block association / interaction measures
out = iispectral(S, iY, iX1, iX2);   % iY, iX1, iX2 = index sets of the three blocks
```

`iispectral` requires `length(iY) + length(iX1) + length(iX2)` to equal the
number of processes in the spectral matrix.

## Reference

If you use this toolbox, please cite:

> Pinto, H., Dias, C., Vergara, V. R., Barà, C., Pernice, R., Rocha, A. P.,
> Faes, L., and Antonacci, Y. (2025). *Time-Frequency Linear Interdependence
> Measures Reveal Multivariate Patterns of Cerebrovascular-Autonomic
> Interactions in Traumatic Brain Injury.* (Manuscript under review — DOI to be
> added upon publication.)

## License

No license file is currently included. Until one is added, all rights are
reserved by default; contact the author regarding reuse.
