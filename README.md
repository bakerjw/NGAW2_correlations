# NGA-West2 Intensity Measure Correlations

This MATLAB codebase analyzes intensity measure correlations from the NGA-West2 ground motion database and reproduces the results presented in:

**Baker, J. W., and Bradley, B. A. (2017). "Intensity measure correlations observed in the NGA-West2 database, and dependence of correlations on rupture and site parameters." *Earthquake Spectra*, 33(1), 145–156.**

## Quick Start
=======
The data were updated in December 2021 to correct errors in the Z1 and region values associated with the ground motion database, and the code now produces figures that are slightly updated relative to the originally published results. Thank you to Eduardo Miranda and Alan Poulos for identifying the errors and helping resolve them.
>>>>>>> b0a2c3c81e085576252448d429a96ba882f02252

1. **Run Analysis**: Execute the main analysis pipeline:
   ```matlab
   mainScript
   ```
   ⚠️ **Note**: The `getResiduals.m` script is computationally intensive and may take significant time to run.

2. **View Results**: Generated figures are saved in `results/figures/` as PDF files.

## Project Structure

```
NGAW2_correlations/
├── README.md                 # This file
├── LICENSE                   # License information
├── mainScript.m             # Main analysis pipeline
├── setup.m                  # MATLAB path configuration
│
├── data/                    # Data files
│   ├── raw/                 # Original NGA-West2 data
│   │   ├── NGA_W2_corr_meta_data.mat
│   │   ├── Z1_Z25_updated.mat
│   │   └── durObs.mat
│   └── processed/           # Generated intermediate data
│       ├── allIMsResids.mat
│       ├── mixedEffectsResids.mat
│       └── rho*.mat
│
├── src/                     # Source code
│   ├── gmpe/               # Ground Motion Prediction Equations
│   ├── residuals/          # Residual computation functions
│   ├── correlations/       # Correlation analysis functions
│   ├── models/             # Correlation prediction models
│   ├── plotting/           # Visualization utilities
│   └── utils/              # General utility functions
│
└── results/                # Generated outputs
    ├── figures/            # PDF figures
    └── supplement/         # Electronic supplement data
```

## Key Components

### Analysis Pipeline
1. **Residual Computation** (`src/residuals/`): Computes within- and between-event residuals using mixed-effects approach
2. **Correlation Analysis** (`src/correlations/`): Estimates correlations from residuals using the approach of Jayaram & Baker (2009)
3. **Model Comparison** (`src/models/`): Compares different correlation prediction models
4. **Visualization** (`src/plotting/`): Generates publication-quality figures

### Ground Motion Models
- **Spectral Acceleration**: BSSA 2014, CY 2014 (modified NGA-West2 implementations)
- **Duration**: Afshari & stewart (2016), Boore et al. (2009)
- **Ratios**: Shahi & Baker (2014) RotD100/RotD50 ratios

### Correlation Models Implemented
- Baker & Jayaram (2008) - Spectral acceleration correlations
- Al Atik (2015) - Updated correlation model
- Cimellaro et al. (2013) - Alternative approach
- Akkar et al. (2014) - ASA model

## Data Requirements

The analysis requires the following data files (included in `data/raw/`):
- `NGA_W2_corr_meta_data.mat`: NGA-West2 flatfile metadata
- `Z1_Z25_updated.mat`: Corrected Z1.0 and Z2.5 values (updated Jan 2021)
- `durObs.mat`: Observed 5-75% and 5-95% duration data

## Updates

### January 2021

The data were updated to correct errors in the Z1 and region values associated with the ground motion database, and the code now produces figures that are slightly updated relative to the originally published results. Thank you to Eduardo Miranda and Alan Poulos for identifying the errors and helping resolve them.

### June 2025

The project was reorganized for usability, with functions put into subfolders by category, and the input and output data and figures put in distinct folders. The computations are unchanged.

## Requirements

- MATLAB and standard toolboxes
- Statistics and Machine Learning Toolbox (for correlation functions)

## Citation

If you use this code, please cite:

```
Baker, J. W., and Bradley, B. A. (2017). "Intensity measure correlations
observed in the NGA-West2 database, and dependence of correlations on
rupture and site parameters." Earthquake Spectra, 33(1), 145–156.
```

## License

This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENSE file for details.

## Support

This code has been provided to document the calculations in the referenced paper and facilitate related research. It is not fully documented or actively maintained.

For questions about the methodology, please refer to the original publication. For technical issues with the code, please check that you have:
1. Run the `setup` script
2. All required data files in `data/raw/`
3. Sufficient computational resources for `getResiduals.m`
