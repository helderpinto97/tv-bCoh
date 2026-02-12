# Time-Varying Block Coherence
This repository provides a toolbox for estimating **time-varying block coherence** by extending Geweke’s linear dependence measures into the **time–frequency domain**. The method combines **block coherence** with **time-varying vector autoregressive (TV-VAR)** modeling to quantify dynamic multivariate interactions.

The toolbox is intended for researchers and practitioners analyzing multichannel signals whose dependency structure changes over time and across frequencies.

## Overview

The implemented framework enables:

- Time–frequency characterization of multivariate linear dependence  
- Block-level interaction analysis between groups of signals  
- Detection of redundant and synergistic effects across variable sets  
- Model-based estimation using TV-VAR dynamics  
- Distribution-free dependence measures at the metric level (no specific signal distribution assumed beyond model estimation)


## Structure
- **Simulation script**: runs the simulations used to demonstrate and validate the proposed method  
- **Functions folder**: contains the core functions implementing the time–frequency and TV-VAR analysis

## Reference
If you use this toolbox, please cite the associated study describing the method:

**Pinto, H.**, Dias, C., Vergara, V. R., Barà, C., Pernice, R., Rocha, A. P., Faes, L., and Antonacci, Y. (2025). Time-Frequency Linear Interdependence Measures Reveal Multivariate Patterns of Cerebrovascular-Autonomic Interactions in Traumatic Brain Injury. *Under Submission*.

