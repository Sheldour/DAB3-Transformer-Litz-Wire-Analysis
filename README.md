# DAB3-Transformer-Litz-Wire-Analysis
MATLAB implementations of analytical models for Litz wire optimization and AC resistance calculation in high-frequency transformers.
# High-Frequency Transformer Litz Wire Analytical Models

This repository contains MATLAB scripts used during the winding design process of a high-frequency transformer for a Dual Active Bridge (DAB3) converter. 

The codes implement classical analytical models for optimizing Litz wire parameters and evaluating winding AC resistance (including skin and proximity effects).

## Features

1. **Optimum Strand Count Calculation (`calculate_nopt_sullivan.m`)**
   * Implements the closed-form physical model proposed by C.R. Sullivan.
   * Calculates the optimal number of Litz wire strands (`n_opt`) and the theoretically optimal strand diameter under specific window geometries and switching frequencies.
   * **Reference:** C. R. Sullivan, "Optimal choice for number of strands in a litz-wire transformer winding," *IEEE Transactions on Power Electronics*, 1999.

2. **AC Resistance Evaluation (`calculate_Rac_Tourkhani_simple.m`)**
   * Implements the analytical AC resistance model by F. Tourkhani and P. Viarouge.
   * Compares computational results using both exact Kelvin functions and Taylor series approximations.
   * **Reference:** F. Tourkhani and P. Viarouge, "Accurate analytical model of winding losses in round Litz wire windings," *IEEE Transactions on Magnetics*, 2001.

## Usage
Simply run the scripts in MATLAB. User inputs (such as switching frequency $f_s$, turns $N$, and core window allocation parameters $b_b, h$) are defined at the top of each script. 

## License
Provided under the [MIT License](LICENSE).
