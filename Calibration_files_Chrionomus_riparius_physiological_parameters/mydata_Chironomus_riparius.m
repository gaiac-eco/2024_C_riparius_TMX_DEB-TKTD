function [data, auxData, metaData, txtData, weights] = mydata_Chironomus_riparius

%% set metaData
metaData.phylum     = 'Arthropoda';
metaData.class      = 'Insecta';
metaData.order      = 'Diptera';
metaData.family     = 'Chironomidae';
metaData.species    = 'Chironomus_riparius';
metaData.species_en = 'Harlequin fly';

metaData.ecoCode.climate = {'Cfb', 'Dfb', 'Dfc'};
metaData.ecoCode.ecozone = {'TH'};
metaData.ecoCode.habitat = {'0eFl', '0eFp', '0eFm', 'eiTg'};
metaData.ecoCode.embryo  = {'Fs'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'bjD'};
metaData.ecoCode.gender  = {'D'};
metaData.ecoCode.reprod  = {'Os'};

metaData.T_typical  = C2K(20); % K, body temp
metaData.data_0     = {'ab'; 'aj_T'; 'ae'; 'am'; 'Lb'; 'Lj'; 'Wd0'; 'Wdj'; 'Wwe'; 'Wde'; 'Ni_f'};
metaData.data_1     = {'t-L_fT'; 'L-Wd'; 't-S'; 'Wd-JO'};

metaData.COMPLETE = 3; % using criteria of LikaKear2011

metaData.author   = {'Starrlight Augustine'};
metaData.date_subm = [2017 09 27];
metaData.email    = {'sta@akvaplan.niva.no'};
metaData.address  = {'Akvaplan-niva, Tromso Norway'};

metaData.author_mod_1   = {'Andre Gergs'};
metaData.date_mod_1     = [2019 08 22];
metaData.email_mod_1    = {'andre.gergs@bayer.com'};
metaData.address_mod_1  = {'Bayer AG, Monheim, Germany'};

metaData.curator     = {'Bas Kooijman'};
metaData.email_cur   = {'bas.kooijman@vu.nl'};
metaData.date_acc    = [2019 09 21];

%% set data
% zero-variate data
data.ab = 2;    units.ab = 'd';    label.ab = 'age at birth';           bibkey.ab = 'PeryMons2002';
temp.ab = C2K(21);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
comment.ab = 'i assume here that the 2-d old individuals used at the start of the experiment were newly hatched 1st instar larvae, SahrRafa2010 also report 2-d at 26 deg C';

% data not included due to direct contradiction with gaiac 2022 data:

% data.tj = 11.8;  units.tj = 'd';    label.tj = 'time since birth at pupation'; bibkey.tj = 'SahrRafa2010';   
%   temp.tj = C2K(21);  units.temp.tj = 'K'; label.temp.tj = 'temperature';
%   comment.tj = 'range 12-33';
% data.tTj_21 = 15;  units.tTj_21 = 'd';    label.tTj_21 = 'time since birth at pupation'; bibkey.tTj_21 = 'PeryGarr2006';   
%   temp.tTj_21 = C2K(21);  units.temp.tTj_21 = 'K'; label.temp.tTj_21 = 'temperature';
% data.tTj_15 = 25;  units.tTj_15 = 'd';    label.tTj_15 = 'time since birth at pupation'; bibkey.tTj_15 = 'PeryGarr2006';   
%   temp.tTj_15 = C2K(15);  units.temp.tTj_15 = 'K'; label.temp.tTj_15 = 'temperature';

data.te = 1;     units.te = 'd';     label.te = 'time since pupation at emergence'; bibkey.te = 'SahrRafa2010';
temp.te = C2K(26);  units.temp.te = 'K'; label.temp.te = 'temperature';
data.am = 1;    units.am = 'd';     label.am = 'life span as imago';       bibkey.am = 'SahrRafa2010';
temp.am = C2K(26);  units.temp.am = 'K'; label.temp.am = 'temperature';
comment.am = 'range 1-3';

data.Lb = 0.17; units.Lb = 'cm'; label.Lb = 'length at birth';      bibkey.Lb = 'PeryMons2002';
data.Lj = 1.38; units.Lj = 'cm'; label.Lj = 'female length of 4th instar larvae before pupation';      bibkey.Lj = 'PeryMons2002';

data.Wd0 = 0.99 *1e-6; units.Wd0 = 'g'; label.Wd0 = 'initial egg ash free dry weight';      bibkey.Wd0 = 'PentHolo1995';
comment.Wd0 = 'the authors work with ash-free dry weight';
data.Wwj = 10*1e-3; units.Wwj = 'g'; label.Wwj = 'fem. max observed wet weight of 4th instar larvae';      bibkey.Wwj = 'SildCran2000';
comment.Wwj = 'max value of fem. control wet wt. from fig. 1';
data.Wde = 1.1*1e-3; units.Wde = 'g'; label.Wde = 'fem. dry weight of imago';     bibkey.Wde = 'RodrGrav2015';
comment.Wde = 'fig. 3d';

% data not included due to direct contradiction with gaiac 2022 data:

% data.Ni  = 509; units.Ni  = '#';  label.Ni  = 'total number of eggs';       bibkey.Ni  = 'SahrRafa2010';
% temp.Ni = C2K(20); units.temp.Ni = 'K'; label.temp.Ni = 'temperature';
% comment.Ni = '1.4 mg/tetramin per ind - Nolte 1993 observed 1800 as max nb eggs';

data.N2  = 151.6; units.N2  = '#';  label.N2  = 'total number of eggs';       bibkey.N2  = 'PeryMons2002';
temp.N2 = C2K(20); units.temp.N2 = 'K'; label.temp.N2 = 'temperature';
comment.N2 = '0.2 mg tetramin/ ind';
data.N3  = 195.588; units.N3  = '#';  label.N3  = 'total number of eggs';       bibkey.N3  = 'PeryMons2002';
temp.N3 = C2K(20); units.temp.N3 = 'K'; label.temp.N3 = 'temperature';
comment.N3 = '0.3 mg tetramin/ ind';
data.N4  = 273.852; units.N4  = '#';  label.N4  = 'total number of eggs';       bibkey.N4  = 'PeryMons2002';
temp.N4 = C2K(20); units.temp.N4 = 'K'; label.temp.N4 = 'temperature';
comment.N4 = '0.4 mg tetramin/ ind';

% data not included due to direct contradiction with gaiac 2022 data:

% % T-aj data
% data.Taj = [ ... % temperature (degC), age at pupation (d)
%     18   35.5  % 'range 30-41'
%     22   21.5  % 'range 13-30'
%     24   21.5  % 'range 12-31'
%     26   22.5  % 'range 12-33'
%     28   35.0  % 'range 14-56'
%     30   26.0  % 'range 14-38'
%     ];
% units.Taj = {'C', 'd'};  label.Taj = {'temperature', 'age at pupation'};
% bibkey.Taj = 'SahrRafa2010';

%--------------------------------------------------------------------------
%---------- Start gaiac 2022 data -----------------------------------------
%--------------------------------------------------------------------------

% zero-variate data

data.N12  = 345.4; units.N12  = '#';  label.N12  = 'total number of eggs 12 degC';  bibkey.N12  = 'Gaia2022';
temp.N12 = C2K(11.8); units.temp.N12 = 'K'; label.temp.N12 = 'temperature';
data.N15  = 367.4; units.N15  = '#';  label.N15  = 'total number of eggs 15 degC';  bibkey.N15  = 'Gaia2022';
temp.N15 = C2K(15.2); units.temp.N15 = 'K'; label.temp.N15 = 'temperature';
data.N20  = 347.6; units.N20  = '#';  label.N20  = 'total number of eggs 20 degC';  bibkey.N20  = 'Gaia2022';
temp.N20 = C2K(20.2); units.temp.N20 = 'K'; label.temp.N20 = 'temperature';

% data not included since so far there is no submodel for temperature stress on repro:

% data.N23  = 301.4; units.N23  = '#';  label.N23  = 'total number of eggs 23 degC';  bibkey.N23  = 'Gaia2022';
% temp.N23 = C2K(23.2); units.temp.N23 = 'K'; label.temp.N23 = 'temperature';

data.ae12  = 42.12; units.ae12  = '#';  label.ae12  = 'age at emergence 12 degC';  bibkey.ae12  = 'Gaia2022';
temp.ae12 = C2K(11.8); units.temp.ae12 = 'd'; label.temp.ae12 = 'temperature';
data.ae15  = 23.87; units.ae15  = '#';  label.ae15  = 'age at emergence 15 degC';  bibkey.ae15  = 'Gaia2022';
temp.ae15 = C2K(15.2); units.temp.ae15 = 'd'; label.temp.ae15 = 'temperature';
data.ae20  = 17.47; units.ae20  = '#';  label.ae20  = 'age at emergence 20 degC';  bibkey.ae20  = 'Gaia2022';
temp.ae20 = C2K(20.2); units.temp.ae20 = 'd'; label.temp.ae20 = 'temperature';
data.ae23  = 15.67; units.ae23  = '#';  label.ae23  = 'age at emergence 23 degC';  bibkey.ae23  = 'Gaia2022';
temp.ae23 = C2K(23.2); units.temp.ae23 = 'd'; label.temp.ae23 = 'temperature';

data.ae = 19.29;    units.ae = 'd';    label.ae = 'age at emergence in extra 20 degC controls';  bibkey.ae = 'Gaia2022';
temp.ae = C2K(20);  units.temp.ae = 'K'; label.temp.ae = 'temperature';

% uni-variate data

% time - weight at 12 degC
data.tW12 = [ ... % time d - weight mg
    16+0.25	 0.039388889
    20+0.25  0.62915
    22+0.25	 0.41945
    26+0.25	 1.44935
    28+0.25	 1.67565
    32+0.25	 2.14155
    40+0.25	 2.379923077
    ];
units.tW12   = {'d', 'mg'};  label.tW12 = {'time since birth', 'dry weight'};
temp.tW12    = C2K(12);  units.temp.tW12 = 'K'; label.temp.tW12 = 'temperature';
bibkey.tW12 = 'Gaia2022';

% time - weight at 15 degC
data.tW15 = [ ... % time d - weight mg
    10+0.25  0.269263158
    12+0.25  0.486125
    15+0.25  0.86715
    17+0.25  1.370772727
    19+0.25  1.510142857
    21+0.25  1.756578947
    23+0.25  1.787846154
    25+0.25  1.851142857
    ];
units.tW15   = {'d', 'mg'};  label.tW15 = {'time since birth', 'dry weight'};
temp.tW15    = C2K(15);  units.temp.tW15 = 'K'; label.temp.tW15 = 'temperature';
bibkey.tW15 = 'Gaia2022';

% time - weight at 20 degC
data.tW20 = [ ... % time d - weight mg
    8+0.25	0.393789474
    10+0.25	0.70105
    11+0.25	0.5592
    12+0.25	0.830380952
    13+0.25	0.903318182
    15+0.25	1.319315789
    18+0.25	1.747384615
    ];
units.tW20   = {'d', 'mg'};  label.tW20 = {'time since birth', 'dry weight'};
temp.tW20    = C2K(20);  units.temp.tW20 = 'K'; label.temp.tW20 = 'temperature';
bibkey.tW20 = 'Gaia2022';

% time - weight at 23 degC
data.tW23 = [ ... % time d - weight mg
    4+0.25   0.0378
    6+0.25   0.122863636
    8+0.25	 0.267473684
    10+0.25  0.902578947
    11+0.25  0.961761905
    13+0.25  1.130833333
    14+0.25  1.1980625
    15+0.25  1.360888889
    ];
units.tW23   = {'d', 'mg'};  label.tW23 = {'time since birth', 'dry weight'};
temp.tW23    = C2K(23);  units.temp.tW23 = 'K'; label.temp.tW23 = 'temperature';
bibkey.tW23 = 'Gaia2022';

%--------------------------------------------------------------------------
%---------- End gaiac 2022 data -------------------------------------------
%--------------------------------------------------------------------------

% uni-variate data

% length - weight data
data.LW = [ ... % cubic length cm^3 - weight mg
    6.189e-001	2.811e-001
    8.045e-001	3.940e-001
    8.045e-001	3.692e-001
    8.341e-001	4.242e-001
    9.480e-001	4.711e-001
    1.166e+000	5.483e-001
    1.193e+000	6.390e-001
    1.403e+000	8.069e-001
    1.435e+000	8.179e-001
    1.421e+000	7.135e-001
    1.520e+000	7.851e-001
    1.527e+000	8.125e-001
    1.644e+000	8.402e-001
    1.819e+000	9.228e-001];
data.LW(:,1) = data.LW(:,1).^(1/3); % cm^3 to cm
units.LW   = {'cm', 'mg'};  label.LW = {'length', 'dry weight'};
temp.LW    = C2K(21);  units.temp.LW = 'K'; label.temp.LW = 'temperature';
bibkey.LW = 'PeryMons2002';
comment.LW = 'data obtained with fourth-instar larvae for two different diets (0.2, 0.3, and 1.4 mg/larvae/d)';

% time - length
data.tL = [ ...
    0   1.81
    1   2.76
    2   3.43
    3   4.922
    3.995	6.273
    5.000	8.540
    6.005	10.747
    6.987	12.803
    7.999	13.677
    8.995	13.823];
data.tL(:,2) = data.tL(:,2)/ 10; % convert cm to  mm
units.tL   = {'d', 'cm'};  label.tL = {'time since birth', 'length', 'female'};
temp.tL    = C2K(21);  units.temp.tL = 'K'; label.temp.tL = 'temperature';
bibkey.tL = 'PeryMons2002';
comment.tL = 'fig 3, fig 4, ad libitum 1.4 mg/ind tetramin, female';

data.tL_m = [ ...
    0   1.81
    1   2.76
    2   3.43
    3   4.922
    3.995	6.273
    5.000	7.965
    6.012	9.566
    7.001	11.106
    7.998	11.404
    9.002	11.399];
data.tL_m(:,2) = data.tL_m(:,2)/ 10; % convert cm to  mm
units.tL_m   = {'d', 'cm'};  label.tL_m = {'time since birth', 'length', 'male'};
temp.tL_m    = C2K(21);  units.temp.tL_m = 'K'; label.temp.tL_m = 'temperature';
bibkey.tL_m = 'PeryMons2002';
comment.tL_m = 'fig 3, fig 4, ad libitum 1.4 mg/ind tetramin, male';

% Eco-physiological information by PeryMons2002
% Chironomidae had been cultured prior to the experiments
% according to standard methods. The first day of the experiment,
% 2-d-old organisms were put in beakers. The volume of these
% beakers was 0.6 L, and the surface area was 14 cm2. They
% contained 0.11 L artificial sediment and 0.44 L water from an
% uncontaminated spring near our laboratory (pH 8.1 and conductivity
% 400 ms/cm).

% time-length data of 4 instars at different food levels

% 0.4 mg tetramin/larae/d
data.tL4 = [...
    0.013	1.795
    2.009	3.354
    3.013	4.608
    4.008	6.477
    5.034	8.064
    6.018	10.472
    7.993	12.236];
data.tL4(:,2) = data.tL4(:,2)/ 10; % convert cm to mm
% data.tL4(:,1) = data.tL4(:,1) + 2; % convert time to day since birth
units.tL4   = {'d', 'cm'};  label.tL4 = {'time since birth', 'length', '0.4 mg/d'};
temp.tL4    = C2K(21);  units.temp.tL4 = 'K'; label.temp.tL4 = 'temperature';
bibkey.tL4 = 'PeryMons2002';
comment.tL4 = 'fig 5';

% 0.3 mg tetramin/ind/d
data.tL3 = [ ...
    0.035	1.769
    2.009	3.251
    3.013	4.685
    4.007	6.015
    5.033	7.526
    6.027	9.472
    8.024	11.544];
data.tL3(:,2) = data.tL3(:,2)/ 10; % convert cm to mm
% data.tL3(:,1) = data.tL3(:,1) + 2; % convert time to day since birth
units.tL3   = {'d', 'cm'};  label.tL3 = {'time since birth', 'length', '0.3 mg/d'};
temp.tL3    = C2K(21);  units.temp.tL3 = 'K'; label.temp.tL3 = 'temperature';
bibkey.tL3 = 'PeryMons2002';
comment.tL3 = 'fig 5';

% 0.2 mg/tetramin/individual
data.tL2 = [ ... % time since start of experiment (d), length(mm)
    0.024	1.744
    2.020	3.354
    3.035	4.505
    3.996	5.605
    5.021	6.936
    6.026	8.472
    8.011	10.313
    10.018	11.923];
data.tL2(:,2) = data.tL2(:,2)/ 10; % convert cm to mm
% data.tL2(:,1) = data.tL2(:,1) + 2; % convert time to day since birth
units.tL2   = {'d', 'cm'};  label.tL2 = {'time since birth', 'length', '0.2 mg/d'};
temp.tL2    = C2K(21);  units.temp.tL2 = 'K'; label.temp.tL2 = 'temperature';
bibkey.tL2 = 'PeryMons2002';
comment.tL2 = 'fig 5';

% 0.1 mg tetramin/ ind/d
data.tL1 = [ ...
    1.999	3.303
    3.013	4.146
    4.006	5.503
    5.031	6.654
    6.024	7.318
    8.020	8.903
    10.015	9.744];
data.tL1(:,2) = data.tL1(:,2)/ 10; % convert cm to mm
% data.tL1(:,1) = data.tL1(:,1) + 2; % convert time to day since birth
units.tL1   = {'d', 'cm'};  label.tL1 = {'time since birth', 'length', '0.1 mg/d'};
temp.tL1    = C2K(21);  units.temp.tL1 = 'K'; label.temp.tL1 = 'temperature';
bibkey.tL1 = 'PeryMons2002';
comment.tL1 = 'fig 5';

% 15 deg C
data.tL15 = [ ...
    0.034	1.424
    2.026	2.221
    4.000	3.090
    4.996	4.103
    6.043	4.828
    7.021	5.721
    8.034	6.807
    9.013	8.521
    9.991	9.534
    11.021	10.645
    12.017	11.345
    12.996	12.045
    14.043	12.648];
data.tL15(:,2) = data.tL15(:,2)/ 10; % convert cm to mm
% data.tL15(:,1) = data.tL15(:,1) + 2; % convert time to day since birth
units.tL15   = {'d', 'cm'};  label.tL15 = {'time since birth', 'length', '15 deg C'};
temp.tL15    = C2K(15);  units.temp.tL15 = 'K'; label.temp.tL15 = 'temperature';
bibkey.tL15 = 'PeryGarr2006';
comment.tL15 = 'fig 1, ad libitum';

% 19.6 deg, ad libitum
data.tL196 = [ ...
    0.052	1.738
    2.026	3.331
    4.017	5.117
    5.013	6.324
    6.043	8.328
    7.039	9.824
    8.017	11.369
    9.013	12.045
    10.026	12.503];
data.tL196(:,2) = data.tL196(:,2)/ 10; % convert cm to mm
% data.tL196(:,1) = data.tL196(:,1) + 2; % convert time to day since birth
units.tL196   = {'d', 'cm'};  label.tL196 = {'time since birth', 'length', '19.6 deg C'};
temp.tL196    = C2K(19.6);  units.temp.tL196 = 'K'; label.temp.tL196 = 'temperature';
bibkey.tL196 = 'PeryGarr2006';
comment.tL196 = 'fig 1, ad libitum';

% 21 deg, ad libitum
data.tL21 = [ ...
    0.000	1.714
    2.009	3.355
    4.000	6.397
    5.013	8.303
    6.026	10.066
    7.021	12.093
    8.017	12.648
    9.030	12.576];
data.tL21(:,2) = data.tL21(:,2)/ 10; % convert cm to mm
% data.tL21(:,1) = data.tL21(:,1) + 2; % convert time to day since birth
units.tL21   = {'d', 'cm'};  label.tL21 = {'time since birth', 'length', '21 deg C'};
temp.tL21    = C2K(21);  units.temp.tL21 = 'K'; label.temp.tL21 = 'temperature';
bibkey.tL21 = 'PeryGarr2006';
comment.tL21 = 'fig 1, ad libitum';

% 24.4 ad libitum
data.tL244 = [ ...
    0.034	1.738
    2.043	3.741
    4.017	7.241
    5.013	10.259
    6.026	11.755
    7.021	12.286
    8.017	12.648
    9.030	12.576];
data.tL244(:,2) = data.tL244(:,2)/ 10; % convert cm to mm
% data.tL244(:,1) = data.tL244(:,1) + 2; % convert time to day since birth
units.tL244   = {'d', 'cm'};  label.tL244 = {'time since birth', 'length', '24.4 deg C'};
temp.tL244    = C2K(24.4);  units.temp.tL244 = 'K'; label.temp.tL244 = 'temperature';
bibkey.tL244 = 'PeryGarr2006';
comment.tL244 = 'fig 1, ad libitum';

% 26.7
data.tL267 = [ ...
    0.000	1.810
    2.043	3.693
    4.017	7.917
    5.013	10.452
    6.043	12.648
    7.073	12.648];
data.tL267(:,2) = data.tL267(:,2)/ 10; % convert cm to mm
% data.tL267(:,1) = data.tL267(:,1) + 2; % convert time to day since birth
units.tL267   = {'d', 'cm'};  label.tL267 = {'time since birth', 'length', '26.7 deg C'};
temp.tL267    = C2K(26.7);  units.temp.tL267 = 'K'; label.temp.tL267 = 'temperature';
bibkey.tL267 = 'PeryGarr2006';
comment.tL267 = 'fig 1, ad libitum';

% different temperatures at limiting food:

% 15 deg C, 0.2 mg tetramin per day
data.tL15f = [ ...
    0.017	1.402
    2.012	2.245
    3.992	3.066
    5.006	4.088
    6.019	4.755
    7.017	5.599
    8.016	7.022
    10.028	9.000
    12.040	11.023
    14.035	12.022];
data.tL15f(:,2) = data.tL15f(:,2)/ 10; % convert cm to mm
% data.tL15f(:,1) = data.tL15f(:,1) + 2; % convert time to day since birth
units.tL15f   = {'d', 'cm'};  label.tL15f = {'time since birth', 'length', '15 deg C'};
temp.tL15f    = C2K(15);  units.temp.tL15f = 'K'; label.temp.tL15f = 'temperature';
bibkey.tL15f = 'PeryGarr2006';
comment.tL15f = 'fig 1, 0.2 mg/d tetramin';

% 19.6 deg, 0.2 mg/tetramen/d
data.tL196f = [ ...
    0.017	1.447
    1.998	3.336
    4.025	5.002
    4.977	6.314
    6.007	7.604
    8.019	9.804
    10.062	11.315];
data.tL196f(:,2) = data.tL196f(:,2)/ 10; % convert cm to mm
% data.tL196f(:,1) = data.tL196f(:,1) + 2; % convert time to day since birth
units.tL196f   = {'d', 'cm'};  label.tL196f = {'time since birth', 'length', '19.6 deg C'};
temp.tL196f    = C2K(19.6);  units.temp.tL196f = 'K'; label.temp.tL196f = 'temperature';
bibkey.tL196f = 'PeryGarr2006';
comment.tL196f = 'fig 3, limiting food 0.2 mg/d';

% 21 deg, 0.2 mg/d
data.tL21f = [ ...
    0.033	1.447
    1.982	3.291
    3.995	6.360
    4.994	7.961
    6.040	8.984
    7.007	10.184
    7.989	11.095
    9.018	11.383];
data.tL21f(:,2) = data.tL21f(:,2)/ 10; % convert cm to mm
% data.tL21f(:,1) = data.tL21f(:,1) + 2; % convert time to day since birth
units.tL21f   = {'d', 'cm'};  label.tL21f = {'time since birth', 'length', '21 deg C'};
temp.tL21f    = C2K(21);  units.temp.tL21f = 'K'; label.temp.tL21f = 'temperature';
bibkey.tL21f = 'PeryGarr2006';
comment.tL21f = 'fig 3, limiting food level 0.2 mg/d';

% 26.7 deg C, limiting
data.tL267f = [ ...
    0.033	1.224
    1.998	3.714
    3.996	6.761
    5.011	8.429
    6.025	9.785
    7.023	10.429
    8.020	10.917
    9.002	11.806];
data.tL267f(:,2) = data.tL267f(:,2)/ 10; % convert cm to mm
% data.tL267f(:,1) = data.tL267f(:,1) + 2; % convert time to day since birth
units.tL267f   = {'d', 'cm'};  label.tL267f = {'time since birth', 'length', '26.7 deg C'};
temp.tL267f    = C2K(26.7);  units.temp.tL267f = 'K'; label.temp.tL267f = 'temperature';
bibkey.tL267f = 'PeryGarr2006';
comment.tL267f = 'fig 3, 0.2 mg/d';

% survival from starvation experiment
data.tS  =[... % t in days and S in number of surviving individuals, -
    0	1
    3	0.44
    4	0.26
    5	0.04
    6	0.01
    7   0];
units.tS = {'d', '-'}; label.tS = {'starvation period', 'fraction surviving', 'fed'};
temp.tS  = C2K(20); units.temp.tS = 'K'; label.temp.tS = 'temperature';
bibkey.tS = {'Baye2019'};
comment.tS = 'mortality of starved chironomids';

% survival from starvation experiment
data.tS2  =[... % t in days and S in number of surviving individuals, -
    0	1
    3	0.96
    4	0.96
    5	0.92
    6   0.86
    7   0.84];
units.tS2 = {'d', '-'}; label.tS2 = {'starvation period', 'fraction surviving', 'starved'};
temp.tS2  = C2K(20); units.temp.tS2 = 'K'; label.temp.tS2 = 'temperature';
bibkey.tS2 = {'Baye2019'};
comment.tS2 = 'background mortality for fed chironomids';

data.WdJO_20  =[... % dry weight in mg - oxygen consumption in mug per mg and hour
    0.12562	6.30013
    0.10866	6.04895
    0.22639	6.17987
    0.26837	6.11979
    0.19685	5.97663
    0.23882	5.89262
    0.24716	5.78487
    0.29756	5.74868
    0.25968	5.64121
    0.24275	5.43789
    0.30561	5.1623
    0.33067	4.92284
    0.29261	4.50426
    0.51948	4.47897
    0.57413	4.52651
    0.69168	4.35828
    0.81369	4.6208
    0.95233	4.60801
    0.98183	4.73945
    0.77541	3.85522
    0.82563	3.51988
    0.9854	3.69841
    1.02305	3.42297
    1.01866	3.09992
    1.09841	2.96782
    1.11961	3.29077
    1.13651	3.44622
    1.12406	3.69758
    1.12832	3.80524
    1.19577	4.16381
    1.25838	3.48139
    1.51056	3.6115
    1.60269	3.10839
    1.72885	3.28712
    1.83368	2.92752
    1.93864	2.79527];
units.WdJO_20 = {'mg', 'mug/h.mg'}; label.WdJO_20 = {'dry weight', 'oxygen consumption', '20 C'};
temp.WdJO_20  = C2K(20); units.temp.WdJO_20 = 'K'; label.temp.WdJO_20 = 'temperature';
bibkey.WdJO_20 = {'Edwa1956'};

data.WdJO_10  =[... % dry weight in mg - oxygen consumption in mug per mg and hour
    0.08386	2.20548
    0.10901	2.34853
    0.14256	2.39591
    0.17191	2.38361
    0.24319	2.43054
    0.25577	2.34677
    0.26834	2.22717
    0.27673	2.13151
    0.23061	2.06039
    0.20126	1.86961
    0.13836	2.02565
    0.17191	2.16859
    0.33543	1.92773
    0.2935	1.70127
    0.37736	1.6286
    0.45702	1.67543
    0.47799	1.85436
    0.53669	1.77004
    0.65409	1.8642
    0.62474	1.67342
    0.587	1.57831
    0.58281	1.47085
    0.72956	1.40937
    0.87631	1.5868
    1.09434	1.59614
    1.06498	1.47704
    1.06499	1.34564
    1.05241	1.23828
    0.95597	0.97664
    1.19078	1.28441
    1.24528	1.23598
    1.27883	1.41476
    1.38784	1.294
    1.47589	1.05404
    1.53459	1.31614
    1.6478	1.00421
    1.7652	1.06253];
units.WdJO_10 = {'mg', 'mug/h.mg'}; label.WdJO_10 = {'dry weight', 'oxygen consumption', '10 C'};
temp.WdJO_10  = C2K(10); units.temp.WdJO_10 = 'K'; label.temp.WdJO_10 = 'temperature';
bibkey.WdJO_10 = {'Edwa1956'};

%% set weights for all real data
weights = setweights(data, []);
% weight settings were determined experimentally after much estimations
% with different starting values

weights.ae = 10 * weights.ae;
weights.ae12 = 10 * weights.ae12;
weights.ae15 = 10 * weights.ae15;
weights.ae20 = 10 * weights.ae20;
weights.ae23 = 10 * weights.ae23;

weights.ab = 3 * weights.ab;
weights.WdJO_10 = 3 * weights.WdJO_10;
weights.WdJO_20 = 3 * weights.WdJO_20;
weights.Lj = 3 * weights.Lj;
weights.LW = 3 * weights.LW;

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);
% this was to help prevent p_M from going to too high values
% i am unsure of the effect of removing this
weights.psd.p_M = 10 * weights.psd.p_M;
weights.psd.v = 10 * weights.psd.v;

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

%% Group plots
set1 = {'tL', 'tL_m'}; comment1 = {'PeryMOns2002 data for females, males'};
set2 = {'tL4', 'tL3', 'tL2', 'tL1'}; comment2 = {'PeryMons2002 data 0.4, 0.3, 0.2 and 0.1 mg/d'};
set3 = {'tL267', 'tL244', 'tL21', 'tL196','tL15'}; comment3 = {'PeryMons2002 ad libitum 26.7, 24.4, 21, 19.6  and 15 deg C'};
set4 = {'tL267f', 'tL21f', 'tL196f','tL15f'}; comment4 = {'PeryMons2002 limited food 26.7, 21, 19.6  and 15 deg C'};
set5 = {'tS', 'tS2'}; comment5 = {'Mortality of fed and starved larvae'};
set6 = {'WdJO_20', 'WdJO_10'}; comment6 = {'Edwa1956 at 10C and 20C'};
set7 = {'tW12', 'tW15', 'tW20', 'tW23'}; comment7 = {'Dry weight at 12, 15, 20, 23 degC'};
metaData.grp.sets = {set1,set2,set3,set4,set5,set6,set7};
metaData.grp.comment = {comment1,comment2,comment3,comment4,comment5,comment6,comment7};

%% Facts
F1 = 'Its life cycle comprises aquatic stages (egg, four larval instars, and a pupal stage) and an aerial adult stage.';
metaData.bibkey.F1 = 'PeryMons2002';
F2 = 'Widely distributed in the northern hemisphere at temperate latitudes. Lentic and lotic environments, usually in organically enriched waters';
metaData.bibkey.F2 = 'PeryMons2002';
F3 = 'larvae, collectors and gatherers, feed on sediment-deposited detritus';
metaData.bibkey.F3 = 'PeryMons2002';
F4 = 'if head capsule width is not taken into account, C. riparius can be considered isomorphic during the larval development.';
metaData.bibkey.F4 = 'PeryMons2002';
F5 = 'Adult females produce 1 egg mass';
metaData.bibkey.F5 = 'SahrRafa2010';
metaData.facts = struct('F1',F1,'F2',F2,'F3',F3, 'F4', F4, 'F5', F5);

%% Discussion points
D1 = 'we assume that all the data from PeryMons2002 and PeryGarr2006 are in time since birth, and that it the 2-d old individuals at the start of the experiment just hatched';
D2 = 'males are assumed to differ from females by {p_Am} only';
D3 = 'I assume time since birth in the time axis of the data';
D4 = 'the maintenance is high, which is a stark contract to the assumption of 0 maintenance costs by the authors of the data used here: PeryMons2002 PeryGarr2006. The mydata file contains heat production data (as well as the reference) which is not yet implemented. One might consider implementing this to see if such high maintenance is consistent with that additional information.';
D5 = 'mod_1: time-Survival and Weight-respiration data added';
metaData.discussion = struct('D1', D1, 'D2', D2, 'D3', D3, 'D4', D4, 'D5', D5);

%% Acknowledgment
metaData.acknowledgment = 'The creation of this entry was supported by the European Food Safety Authority (grant number OC/EFSA/SCER/2015/01)';

%% Links
metaData.links.id_CoL = '8BGDM'; % Cat of Life
metaData.links.id_ITIS = '129313'; % ITIS
metaData.links.id_EoL = '745125'; % Ency of Life
metaData.links.id_Wiki = 'Chironomus_riparius'; % Wikipedia
metaData.links.id_ADW = 'Chironomus_riparius'; % ADW
metaData.links.id_Taxo = '315139'; % Taxonomicon
metaData.links.id_WoRMS = ''; % WoRMS
metaData.links.id_diptera = 'Chironomus+riparius'; % Diptera

%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
    'howpublished = {\url{https://en.wikipedia.org/wiki/Chironomus_riparius}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
    'author = {Kooijman, S.A.L.M.}, ' ...
    'year = {2010}, ' ...
    'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
    'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
    'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
    'howpublished = {\url{../../../bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'PeryMons2002'; type = 'Article'; bib = [ ...
    'author = {P\''{e}ry, Alexandre R. R. and Mons, Rapha\"{e}l and Flammarion, Patrick and Lagadic, Laurent and Garric, Jeanne}, ' ...
    'year = {2002}, ' ...
    'title = {A modeling approach to link food availability, growth, emergence, and reproduction for the midge \emph{Chironomus riparius}}, ' ...
    'journal = {Environmental Toxicology and Chemistry}, ' ...
    'volume = {21}, ' ...
    'number = {11}, ' ...
    'pages = {2507--2513}, '...
    'doi = {10.1002/etc.5620211133}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'PeryGarr2006'; type = 'Article'; bib = [ ...
    'author = {P\''{e}ry, Alexandre and R. R. Garric, Jeanne}, ' ...
    'year = {2006}, ' ...
    'title = {Modelling Effects of Temperature and Feeding Level on the Life Cycle of the Midge \emph{Chironomus riparius}: {A}n Energy-Based Modelling Approach}, ' ...
    'journal = {Hydrobiologia}, ' ...
    'volume = {553}, ' ...
    'number = {1}, ' ...
    'pages = {59}, '...
    'doi = {10.1007/s10750-005-1284-0}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'SahrRafa2010'; type = 'Article'; bib = [ ...
    'author = {A. Sahragard and M. Rafatifard}, ' ...
    'year = {2006}, ' ...
    'title = {BIOLOGY AND EFFECT OF TEMPERATURE ON LARVAL DEVELOPMENT TIME OF \emph{Chironomus riparius} {M}EIGEN ({D}IPTERA: {C}HIRONOMIDAE) UNDER LABORATORY CONDITIONS}, ' ...
    'journal = {Munis Entomology \& Zoology}, ' ...
    'volume = {5}, ' ...
    'pages = {1025--1033}, '...
    'url = {http://www.munisentzool.org/yayin/vol5/suppl/1025-1033.pdf}';];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'SildCran2000'; type = 'Article'; bib = [ ...
    'author = {Sildanchandra, W. and Crane, M.}, ' ...
    'year = {2000}, ' ...
    'title = {Influence of sexual dimorphism in \emph{Chironomus riparius} {M}eigen on toxic effects of cadmium}, ' ...
    'journal = {Environmental Toxicology and Chemistry}, ' ...
    'volume = {19}, ' ...
    'number = {9}, ' ...
    'pages = {2309--2313}, '...
    'doi = {10.1002/etc.5620190921}';];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'RodrGrav2015'; type = 'Article'; bib = [ ...
    'author = {A. C. M.. Rodrigues and C. Gravato and C. Quintaneiro and C. Barata and A. M. V. M. Soares and J. L. T. Pestana}, ' ...
    'year = {2015}, ' ...
    'title = {Sub-lethal toxicity of environmentally relevant concentrations of esfenvalerate to \emph{Chironomus riparius}}, ' ...
    'journal = {Environmental Pollution}, ' ...
    'volume = {207}, ' ...
    'number = {9}, ' ...
    'pages = {273--279}, '...
    'doi = {10.1016/j.envpol.2015.09.035}';];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'PentHolo1995'; type = 'Article'; bib = [ ...
    'author = {Penttinen, O.-P. and Holopainen, I. J.}, ' ...
    'year = {1995}, ' ...
    'title = {Physiological energetics of a midge, \emph{Chironomus riparius} {M}eigen ({I}nsecta, {D}iptera): normoxic heat output over the whole life cycle and response of larva to hypoxia and anoxia}, ' ...
    'journal = {Oecologia}, ' ...
    'volume = {103}, ' ...
    'number = {4}, ' ...
    'pages = {419--424}, '...
    'doi = {10.1007/BF00328679}';];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Edwa1956'; type = 'Article'; bib = [ ...
    'author = {Edwards, R.W.}, ' ...
    'year = {1958}, ' ...
    'title = {The relation of oxygen consumption to body Size and to temperature in the larvae of \emph{Chironomus riparius} {M}eigen}, ' ...
    'journal = {Journal of Experimental Biology}, ' ...
    'volume = {35}, ' ...
    'number = {-}, ' ...
    'pages = {383--395}, '...
    'doi = {-}';];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Baye2019'; type = 'Misc'; bib = [ ...
    'year = {2019}, ' ...
    'howpublished = {unpblished data}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Gaia2022'; type = 'Misc'; bib = [ ...
    'year = {2022}, ' ...
    'howpublished = {unpblished data}']; % Can be updated once the paper is published
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

