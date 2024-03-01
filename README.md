# 2024_C_riparius_TMX_DEB-TKTD

DEBtool calibration files for:

(a) the physiological DEB model of Chironomus riparius fitted to literature data and new temperature-dependent life-history data 
(b) TKTD parameters fitted to new temperature-dependent lethal and sublethal effect data of thiamethoxam in C. riparius

as presented by Koch et al. 2024

<insert full publication data upon publication>

DEBtool uses a solving algorithm that iteratively adjusts parameter values to minimize the difference between the experimental data and the model predictions. The applied solving algorithm in this study was the Nelder-Mead simplex method, which is DEBtool’s default solving algorithm. The difference between the data and the model predictions was computed using a loss function, evaluating all data simultaneously. The loss function selected for this study was the multiplicative symmetric bounded loss function, which is DEBtool’s default loss function.
To execute the code, a licensed version of MATLAB and the free DEBtool package are required.
The model code consists of four separate files for the calibration:
(1) the ‘mydata’ file containing the data (including references) used for calibration;
(2) the ‘pars_init’ file containing the parameter values;
(3) the ‘predict’ file containing the model equations needed to make life-history predictions based on the parameters, and link them to the observed data;
(4) the ‘run’ file, which calls the necessary functions and which allows for custom settings concerning the estimation procedure.
To use the code, a path to the DEBtool package must be set in the MATLAB environment.
To run the model, the ‘Run’ button must be pressed with the ‘run’ file in the foreground.
To simply run the calibrated model without a new parameter estimation, ‘method’ (in the ‘run’ file) must be set to ‘no’. 
To perform a parameter estimation, the option ‘method’ (in the ‘run’ file) must be set to ‘nm’ and the property ‘free’ of the parameters which should be estimated (in the ‘pars-init’ file) must be set to 1.
When the simulation has finished, the data predictions and the estimated parameters will be shown.

More information about the use of DEBtool can be found via:
https://debtool.debtheory.org/docs/index.html

IMPORTANT:
Please note that these files were never
