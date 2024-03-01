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
metaData.data_0     = {};
metaData.data_1     = {};

%% set data

%--------------------------------------------------------------------------
%---------- Start gaiac 2022 data -----------------------------------------
%--------------------------------------------------------------------------

% zero-variate data

% emergence 12 degC
data.ae12 = 42.12;   units.ae12 = 'd';   label.ae12 = 'age at emergence';   bibkey.ae12 = 'Gaia2022';
temp.ae12 = C2K(11.8);   units.temp.ae12 = 'K';   label.temp.ae12 = 'temperature';
% emergence 15 degC
data.ae15 = 23.87;   units.ae15 = 'd';   label.ae15 = 'age at emergence';   bibkey.ae15 = 'Gaia2022';
temp.ae15 = C2K(15.2);   units.temp.ae15 = 'K';   label.temp.ae15 = 'temperature';
data.ae15_c5 = 23.66;   units.ae15_c5 = 'd';   label.ae15_c5 = 'age at emergence';   bibkey.ae15_c5 = 'Gaia2022';
temp.ae15_c5 = C2K(15.2);   units.temp.ae15_c5 = 'K';   label.temp.ae15_c5 = 'temperature';
% emergence 20 degC
data.ae20 = 17.47;   units.ae20 = 'd';   label.ae20 = 'age at emergence';   bibkey.ae20 = 'Gaia2022';
temp.ae20 = C2K(20.2);   units.temp.ae20 = 'K';   label.temp.ae20 = 'temperature';
data.ae20_c4 = 18.62;   units.ae20_c4 = 'd';   label.ae20_c4 = 'age at emergence';   bibkey.ae20_c4 = 'Gaia2022';
temp.ae20_c4 = C2K(20.2);   units.temp.ae20_c4 = 'K';   label.temp.ae20_c4 = 'temperature';
data.ae20_c5 = 18.87;   units.ae20_c5 = 'd';   label.ae20_c5 = 'age at emergence';   bibkey.ae20_c5 = 'Gaia2022';
temp.ae20_c5 = C2K(20.2);   units.temp.ae20_c5 = 'K';   label.temp.ae20_c5 = 'temperature';
% emergence 23 degC
data.ae23 = 15.67;   units.ae23 = 'd';   label.ae23 = 'age at emergence';   bibkey.ae23 = 'Gaia2022';
temp.ae23 = C2K(23.2);   units.temp.ae23 = 'K';   label.temp.ae23 = 'temperature';
data.ae23_c5 = 14.47;   units.ae23_c5 = 'd';   label.ae23_c5 = 'age at emergence';   bibkey.ae23_c5 = 'Gaia2022';
temp.ae23_c5 = C2K(23.2);   units.temp.ae23_c5 = 'K';   label.temp.ae23_c5 = 'temperature';

% brood size 12 degC
data.Ni12 = 345.4;   units.Ni12 = '#';   label.Ni12 = 'brood size';   bibkey.Ni12 = 'Gaia2022';
temp.Ni12 = C2K(11.8);   units.temp.Ni12 = 'K';   label.temp.Ni12 = 'temperature';
data.Ni12_c5 = 308;   units.Ni12_c5 = '#';   label.Ni12_c5 = 'brood size';   bibkey.Ni12_c5 = 'Gaia2022';
temp.Ni12_c5 = C2K(11.8);   units.temp.Ni12_c5 = 'K';   label.temp.Ni12_c5 = 'temperature';
% brood size 15 degC
data.Ni15 = 367.4;   units.Ni15 = '#';   label.Ni15 = 'brood size';   bibkey.Ni15 = 'Gaia2022';
temp.Ni15 = C2K(15.2);   units.temp.Ni15 = 'K';   label.temp.Ni15 = 'temperature';
data.Ni15_c5 = 320.0;   units.Ni15_c5 = '#';   label.Ni15_c5 = 'brood size';   bibkey.Ni15_c5 = 'Gaia2022';
temp.Ni15_c5 = C2K(15.2);   units.temp.Ni15_c5 = 'K';   label.temp.Ni15_c5 = 'temperature';
% brood size 20 degC
data.Ni20 = 347.6;   units.Ni20 = '#';   label.Ni20 = 'brood size';   bibkey.Ni20 = 'Gaia2022';
temp.Ni20 = C2K(20.2);   units.temp.Ni20 = 'K';   label.temp.Ni20 = 'temperature';
data.Ni20_c4 = 311.4;   units.Ni20_c4 = '#';   label.Ni20_c4 = 'brood size';   bibkey.Ni20_c4 = 'Gaia2022';
temp.Ni20_c4 = C2K(20.2);   units.temp.Ni20_c4 = 'K';   label.temp.Ni20_c4 = 'temperature';
data.Ni20_c5 = 323.8;   units.Ni20_c5 = '#';   label.Ni20_c5 = 'brood size';   bibkey.Ni20_c5 = 'Gaia2022';
temp.Ni20_c5 = C2K(20.2);   units.temp.Ni20_c5 = 'K';   label.temp.Ni20_c5 = 'temperature';
% brood size 23 degC
data.Ni23 = 301.4;   units.Ni23 = '#';   label.Ni23 = 'brood size';   bibkey.Ni23 = 'Gaia2022';
temp.Ni23 = C2K(23.2);   units.temp.Ni23 = 'K';   label.temp.Ni23 = 'temperature';
data.Ni23_c5 = 273.2;   units.Ni23_c5 = '#';   label.Ni23_c5 = 'brood size';   bibkey.Ni23_c5 = 'Gaia2022';
temp.Ni23_c5 = C2K(23.2);   units.temp.Ni23_c5 = 'K';   label.temp.Ni23_c5 = 'temperature';

% uni-variate data

% t-S chronic data at 12 degC
data.tSc_12 = [ ... % time since hatch (d), survival absolute (#) in c_tmx = 0, 18, 25, 32, 39, 46 mug/L
0   20  20  20  20  20  20  
16	18	22	21	17	9	5
20	20	19	20	10	7	7
22	20	22	17	13	8	0
26	20	21	18	5	1	5
28	20	22	18	2	1	3
32	20	18	11	9	0	0
40	19	13	0	0	2	0
44	20	2	4	0	0	0
   ];
data.tSc_12(:,2:7) = data.tSc_12(:,2:7) / 20; % convert absolute to relative
data.tSc_12_c1 = [data.tSc_12(:,1), data.tSc_12(:,7)]; % 46 mug/L
 units.tSc_12_c1 = {'d', '-'};  label.tSc_12_c1 = {'time', 'survival'};
 temp.tSc_12_c1 = C2K(11.8);  units.temp.tSc_12_c1 = 'K';  label.temp.tSc_12_c1 = 'temperature';
 bibkey.tSc_12_c1 = 'Gaia2022';  comment.tSc_12_c1 = '46 mug/L';
data.tSc_12_c2 = [data.tSc_12(:,1), data.tSc_12(:,6)]; % 39 mug/L
 units.tSc_12_c2 = {'d', '-'};  label.tSc_12_c2 = {'time', 'survival'};
 temp.tSc_12_c2 = C2K(11.8);  units.temp.tSc_12_c2 = 'K';  label.temp.tSc_12_c2 = 'temperature';
 bibkey.tSc_12_c2 = 'Gaia2022';  comment.tSc_12_c2 = '39 mug/L';
data.tSc_12_c3 = [data.tSc_12(:,1), data.tSc_12(:,5)]; % 32 mug/L
 units.tSc_12_c3 = {'d', '-'};  label.tSc_12_c3 = {'time', 'survival'};
 temp.tSc_12_c3 = C2K(11.8);  units.temp.tSc_12_c3 = 'K';  label.temp.tSc_12_c3 = 'temperature';
 bibkey.tSc_12_c3 = 'Gaia2022';  comment.tSc_12_c3 = '32 mug/L';
data.tSc_12_c4 = [data.tSc_12(:,1), data.tSc_12(:,4)]; % 25 mug/L
 units.tSc_12_c4 = {'d', '-'};  label.tSc_12_c4 = {'time', 'survival'};
 temp.tSc_12_c4 = C2K(11.8);  units.temp.tSc_12_c4 = 'K';  label.temp.tSc_12_c4 = 'temperature';
 bibkey.tSc_12_c4 = 'Gaia2022';  comment.tSc_12_c4 = '25 mug/L';
data.tSc_12_c5 = [data.tSc_12(:,1), data.tSc_12(:,3)]; % 18 mug/L
 units.tSc_12_c5 = {'d', '-'};  label.tSc_12_c5 = {'time', 'survival'};
 temp.tSc_12_c5 = C2K(11.8);  units.temp.tSc_12_c5 = 'K';  label.temp.tSc_12_c5 = 'temperature';
 bibkey.tSc_12_c5 = 'Gaia2022';  comment.tSc_12_c5 = '18 mug/L';
data.tSc_12 = [data.tSc_12(:,1), data.tSc_12(:,2)]; % control
 units.tSc_12 = {'d', '-'};  label.tSc_12 = {'time', 'survival'};
 temp.tSc_12 = C2K(11.8);  units.temp.tSc_12 = 'K';  label.temp.tSc_12 = 'temperature';
 bibkey.tSc_12 = 'Gaia2022';  comment.tSc_12 = 'control treatment';
 
 % t-S chronic data at 15 degC
data.tSc_15 = [ ... % time since hatch (d), survival absolute (#) in c_tmx = 0, 11, 18, 25, 32, 39 mug/L
0   20  20  20  20  20  20  
10	19	19	19	20	18	15
12	24	19	18	20	16	18
15	20	19	18	20	6	3
17	22	18	18	8	4	2
19	21	19	15	10	6	1
21	19	17	18	9	4	0
23	20	17	14	5	3	0
25	21	12	12	0	0	0
   ];
data.tSc_15(:,2:7) = data.tSc_15(:,2:7) / 20; % convert absolute to relative
data.tSc_15_c1 = [data.tSc_15(:,1), data.tSc_15(:,7)]; % 39 mug/L
 units.tSc_15_c1 = {'d', '-'};  label.tSc_15_c1 = {'time', 'survival'};
 temp.tSc_15_c1 = C2K(15.2);  units.temp.tSc_15_c1 = 'K';  label.temp.tSc_15_c1 = 'temperature';
 bibkey.tSc_15_c1 = 'Gaia2022';  comment.tSc_15_c1 = '39 mug/L';
data.tSc_15_c2 = [data.tSc_15(:,1), data.tSc_15(:,6)]; % 32 mug/L
 units.tSc_15_c2 = {'d', '-'};  label.tSc_15_c2 = {'time', 'survival'};
 temp.tSc_15_c2 = C2K(15.2);  units.temp.tSc_15_c2 = 'K';  label.temp.tSc_15_c2 = 'temperature';
 bibkey.tSc_15_c2 = 'Gaia2022';  comment.tSc_15_c2 = '32 mug/L';
data.tSc_15_c3 = [data.tSc_15(:,1), data.tSc_15(:,5)]; % 25 mug/L
 units.tSc_15_c3 = {'d', '-'};  label.tSc_15_c3 = {'time', 'survival'};
 temp.tSc_15_c3 = C2K(15.2);  units.temp.tSc_15_c3 = 'K';  label.temp.tSc_15_c3 = 'temperature';
 bibkey.tSc_15_c3 = 'Gaia2022';  comment.tSc_15_c3 = '25 mug/L';
data.tSc_15_c4 = [data.tSc_15(:,1), data.tSc_15(:,4)]; % 18 mug/L
 units.tSc_15_c4 = {'d', '-'};  label.tSc_15_c4 = {'time', 'survival'};
 temp.tSc_15_c4 = C2K(15.2);  units.temp.tSc_15_c4 = 'K';  label.temp.tSc_15_c4 = 'temperature';
 bibkey.tSc_15_c4 = 'Gaia2022';  comment.tSc_15_c4 = '18 mug/L';
data.tSc_15_c5 = [data.tSc_15(:,1), data.tSc_15(:,3)]; % 11 mug/L
 units.tSc_15_c5 = {'d', '-'};  label.tSc_15_c5 = {'time', 'survival'};
 temp.tSc_15_c5 = C2K(15.2);  units.temp.tSc_15_c5 = 'K';  label.temp.tSc_15_c5 = 'temperature';
 bibkey.tSc_15_c5 = 'Gaia2022';  comment.tSc_15_c5 = '11 mug/L';
data.tSc_15 = [data.tSc_15(:,1), data.tSc_15(:,2)]; % control
 units.tSc_15 = {'d', '-'};  label.tSc_15 = {'time', 'survival'};
 temp.tSc_15 = C2K(15.2);  units.temp.tSc_15 = 'K';  label.temp.tSc_15 = 'temperature';
 bibkey.tSc_15 = 'Gaia2022';  comment.tSc_15 = 'control treatment';
 
% t-S chronic data at 20 degC
data.tSc_20 = [ ... % time since hatch (d), survival absolute (#) in c_tmx = 0, 4, 11, 18, 25, 32 mug/L
0   20  20  20  20  20  20   
8	19	22	21	21	18	7
10	20	19	20	23	19	7
11	25	18	21	19	13	9
12	21	19	19	24	11	3
13	22	20	19	13	7	6
15	19	22	22	18	10	0
18	17	15	14	7	4	0
   ];
data.tSc_20(:,2:7) = data.tSc_20(:,2:7) / 20; % convert absolute to relative
data.tSc_20_c1 = [data.tSc_20(:,1), data.tSc_20(:,7)]; % 32 mug/L
 units.tSc_20_c1 = {'d', '-'};  label.tSc_20_c1 = {'time', 'survival'};
 temp.tSc_20_c1 = C2K(20.2);  units.temp.tSc_20_c1 = 'K';  label.temp.tSc_20_c1 = 'temperature';
 bibkey.tSc_20_c1 = 'Gaia2022';  comment.tSc_20_c1 = '32 mug/L';
data.tSc_20_c2 = [data.tSc_20(:,1), data.tSc_20(:,6)]; % 25 mug/L
 units.tSc_20_c2 = {'d', '-'};  label.tSc_20_c2 = {'time', 'survival'};
 temp.tSc_20_c2 = C2K(20.2);  units.temp.tSc_20_c2 = 'K';  label.temp.tSc_20_c2 = 'temperature';
 bibkey.tSc_20_c2 = 'Gaia2022';  comment.tSc_20_c2 = '25 mug/L';
data.tSc_20_c3 = [data.tSc_20(:,1), data.tSc_20(:,5)]; % 18 mug/L
 units.tSc_20_c3 = {'d', '-'};  label.tSc_20_c3 = {'time', 'survival'};
 temp.tSc_20_c3 = C2K(20.2);  units.temp.tSc_20_c3 = 'K';  label.temp.tSc_20_c3 = 'temperature';
 bibkey.tSc_20_c3 = 'Gaia2022';  comment.tSc_20_c3 = '18 mug/L';
data.tSc_20_c4 = [data.tSc_20(:,1), data.tSc_20(:,4)]; % 11 mug/L
 units.tSc_20_c4 = {'d', '-'};  label.tSc_20_c4 = {'time', 'survival'};
 temp.tSc_20_c4 = C2K(20.2);  units.temp.tSc_20_c4 = 'K';  label.temp.tSc_20_c4 = 'temperature';
 bibkey.tSc_20_c4 = 'Gaia2022';  comment.tSc_20_c4 = '11 mug/L';
data.tSc_20_c5 = [data.tSc_20(:,1), data.tSc_20(:,3)]; % 4 mug/L
 units.tSc_20_c5 = {'d', '-'};  label.tSc_20_c5 = {'time', 'survival'};
 temp.tSc_20_c5 = C2K(20.2);  units.temp.tSc_20_c5 = 'K';  label.temp.tSc_20_c5 = 'temperature';
 bibkey.tSc_20_c5 = 'Gaia2022';  comment.tSc_20_c5 = '4 mug/L';
data.tSc_20 = [data.tSc_20(:,1), data.tSc_20(:,2)]; % control
 units.tSc_20 = {'d', '-'};  label.tSc_20 = {'time', 'survival'};
 temp.tSc_20 = C2K(20.2);  units.temp.tSc_20 = 'K';  label.temp.tSc_20 = 'temperature';
 bibkey.tSc_20 = 'Gaia2022';  comment.tSc_20 = 'control treatment';
 
% t-S chronic data at 23 degC
data.tSc_23 = [ ... % time since hatch (d), survival absolute (#) in c_tmx = 0, 4, 11, 18, 25, 32 mug/L
0   20  20  20  20  20  20  
4	20	19	18	19	18	1
6	22	20	19	20	8	4
8	19	20	20	18	6	0
10	19	20	19	2	0	0
11	21	20	11	5	0	0
13	20	18	9	3	0	0
14	20	16	4	3	0	0
15	20	13	2	1	0	0
   ];
data.tSc_23(:,2:7) = data.tSc_23(:,2:7) / 20; % convert absolute to relative
data.tSc_23_c1 = [data.tSc_23(:,1), data.tSc_23(:,7)]; % 32 mug/L
 units.tSc_23_c1 = {'d', '-'};  label.tSc_23_c1 = {'time', 'survival'};
 temp.tSc_23_c1 = C2K(23.2);  units.temp.tSc_23_c1 = 'K';  label.temp.tSc_23_c1 = 'temperature';
 bibkey.tSc_23_c1 = 'Gaia2022';  comment.tSc_23_c1 = '32 mug/L';
data.tSc_23_c2 = [data.tSc_23(:,1), data.tSc_23(:,6)]; % 25 mug/L
 units.tSc_23_c2 = {'d', '-'};  label.tSc_23_c2 = {'time', 'survival'};
 temp.tSc_23_c2 = C2K(23.2);  units.temp.tSc_23_c2 = 'K';  label.temp.tSc_23_c2 = 'temperature';
 bibkey.tSc_23_c2 = 'Gaia2022';  comment.tSc_23_c2 = '25 mug/L';
data.tSc_23_c3 = [data.tSc_23(:,1), data.tSc_23(:,5)]; % 18 mug/L
 units.tSc_23_c3 = {'d', '-'};  label.tSc_23_c3 = {'time', 'survival'};
 temp.tSc_23_c3 = C2K(23.2);  units.temp.tSc_23_c3 = 'K';  label.temp.tSc_23_c3 = 'temperature';
 bibkey.tSc_23_c3 = 'Gaia2022';  comment.tSc_23_c3 = '18 mug/L';
data.tSc_23_c4 = [data.tSc_23(:,1), data.tSc_23(:,4)]; % 11 mug/L
 units.tSc_23_c4 = {'d', '-'};  label.tSc_23_c4 = {'time', 'survival'};
 temp.tSc_23_c4 = C2K(23.2);  units.temp.tSc_23_c4 = 'K';  label.temp.tSc_23_c4 = 'temperature';
 bibkey.tSc_23_c4 = 'Gaia2022';  comment.tSc_23_c4 = '11 mug/L';
data.tSc_23_c5 = [data.tSc_23(:,1), data.tSc_23(:,3)]; % 4 mug/L
 units.tSc_23_c5 = {'d', '-'};  label.tSc_23_c5 = {'time', 'survival'};
 temp.tSc_23_c5 = C2K(23.2);  units.temp.tSc_23_c5 = 'K';  label.temp.tSc_23_c5 = 'temperature';
 bibkey.tSc_23_c5 = 'Gaia2022';  comment.tSc_23_c5 = '4 mug/L';
data.tSc_23 = [data.tSc_23(:,1), data.tSc_23(:,2)]; % control
 units.tSc_23 = {'d', '-'};  label.tSc_23 = {'time', 'survival'};
 temp.tSc_23 = C2K(23.2);  units.temp.tSc_23 = 'K';  label.temp.tSc_23 = 'temperature';
 bibkey.tSc_23 = 'Gaia2022';  comment.tSc_23 = 'control treatment';
 
 % t-W chronic data at 12 degC
data.tWc_12 = [ ... % time since hatch (d), dry weight (mg) in c_tmx = 0, 18, 25, 32, 39, 46 mug/L
16	0.039	0.572	0.126	0.121	0.087	0.574
20	0.629	0.895	0.639	0.147	0.033	0.562
22	0.419	0.898	0.921	0.159	0.015	NaN
26	1.449	1.595	1.351	0.461	NaN     0.150
28	1.676	1.738	1.405	0.640	NaN 	0.718
32	2.142	1.386	1.301	0.796	NaN 	NaN
40	2.380	1.669	NaN	    NaN     NaN 	NaN
44	NaN	    NaN     1.174	NaN     NaN 	NaN
   ];
data.tWc_12_c1 = [data.tWc_12([1 2 4 5],1), data.tWc_12([1 2 4 5],7)]; % 46 mug/L
 units.tWc_12_c1 = {'d', 'mg'};  label.tWc_12_c1 = {'time', 'dry weight'};
 temp.tWc_12_c1 = C2K(11.8);  units.temp.tWc_12_c1 = 'K';  label.temp.tWc_12_c1 = 'temperature';
 bibkey.tWc_12_c1 = 'Gaia2022';  comment.tWc_12_c1 = '46 mug/L';
data.tWc_12_c2 = [data.tWc_12(1:3,1), data.tWc_12(1:3,6)]; % 39 mug/L
 units.tWc_12_c2 = {'d', 'mg'};  label.tWc_12_c2 = {'time', 'dry weight'};
 temp.tWc_12_c2 = C2K(11.8);  units.temp.tWc_12_c2 = 'K';  label.temp.tWc_12_c2 = 'temperature';
 bibkey.tWc_12_c2 = 'Gaia2022';  comment.tWc_12_c2 = '39 mug/L';
data.tWc_12_c3 = [data.tWc_12(1:6,1), data.tWc_12(1:6,5)]; % 32 mug/L
 units.tWc_12_c3 = {'d', 'mg'};  label.tWc_12_c3 = {'time', 'dry weight'};
 temp.tWc_12_c3 = C2K(11.8);  units.temp.tWc_12_c3 = 'K';  label.temp.tWc_12_c3 = 'temperature';
 bibkey.tWc_12_c3 = 'Gaia2022';  comment.tWc_12_c3 = '32 mug/L';
data.tWc_12_c4 = [data.tWc_12([1:6 8],1), data.tWc_12([1:6 8],4)]; % 25 mug/L
 units.tWc_12_c4 = {'d', 'mg'};  label.tWc_12_c4 = {'time', 'dry weight'};
 temp.tWc_12_c4 = C2K(11.8);  units.temp.tWc_12_c4 = 'K';  label.temp.tWc_12_c4 = 'temperature';
 bibkey.tWc_12_c4 = 'Gaia2022';  comment.tWc_12_c4 = '25 mug/L';
data.tWc_12_c5 = [data.tWc_12(1:7,1), data.tWc_12(1:7,3)]; % 18 mug/L
 units.tWc_12_c5 = {'d', 'mg'};  label.tWc_12_c5 = {'time', 'dry weight'};
 temp.tWc_12_c5 = C2K(11.8);  units.temp.tWc_12_c5 = 'K';  label.temp.tWc_12_c5 = 'temperature';
 bibkey.tWc_12_c5 = 'Gaia2022';  comment.tWc_12_c5 = '18 mug/L';
data.tWc_12 = [data.tWc_12(1:7,1), data.tWc_12(1:7,2)]; % control
 units.tWc_12 = {'d', 'mg'};  label.tWc_12 = {'time', 'dry weight'};
 temp.tWc_12 = C2K(11.8);  units.temp.tWc_12 = 'K';  label.temp.tWc_12 = 'temperature';
 bibkey.tWc_12 = 'Gaia2022';  comment.tWc_12 = 'control treatment';

% t-W chronic data at 15 degC
data.tWc_15 = [ ... % time since hatch (d), dry weight (mg) in c_tmx = 0, 11, 18, 25, 32, 39 mug/L
10	0.269	0.116	0.175	0.153	0.198	0.094
12	0.486	0.226	0.235	0.304	0.272	0.150
15	0.867	0.808	0.965	0.690	0.205	0.251
17	1.371	1.475	1.073	0.845	0.194	0.434
19	1.510	1.333	1.255	0.952	0.659	NaN
21	1.757	1.283	0.937	0.712	0.509	NaN
23	1.788	1.734	1.366	0.966	0.612	NaN
25	1.851	2.130	1.321	NaN	    NaN 	NaN
   ];
data.tWc_15_c1 = [data.tWc_15(1:4,1), data.tWc_15(1:4,7)]; % 39 mug/L
 units.tWc_15_c1 = {'d', 'mg'};  label.tWc_15_c1 = {'time', 'dry weight'};
 temp.tWc_15_c1 = C2K(15.2);  units.temp.tWc_15_c1 = 'K';  label.temp.tWc_15_c1 = 'temperature';
 bibkey.tWc_15_c1 = 'Gaia2022';  comment.tWc_15_c1 = '39 mug/L';
data.tWc_15_c2 = [data.tWc_15(1:7,1), data.tWc_15(1:7,6)]; % 32 mug/L
 units.tWc_15_c2 = {'d', 'mg'};  label.tWc_15_c2 = {'time', 'dry weight'};
 temp.tWc_15_c2 = C2K(15.2);  units.temp.tWc_15_c2 = 'K';  label.temp.tWc_15_c2 = 'temperature';
 bibkey.tWc_15_c2 = 'Gaia2022';  comment.tWc_15_c2 = '32 mug/L';
data.tWc_15_c3 = [data.tWc_15(1:7,1), data.tWc_15(1:7,5)]; % 25 mug/L
 units.tWc_15_c3 = {'d', 'mg'};  label.tWc_15_c3 = {'time', 'dry weight'};
 temp.tWc_15_c3 = C2K(15.2);  units.temp.tWc_15_c3 = 'K';  label.temp.tWc_15_c3 = 'temperature';
 bibkey.tWc_15_c3 = 'Gaia2022';  comment.tWc_15_c3 = '25 mug/L';
data.tWc_15_c4 = [data.tWc_15(:,1), data.tWc_15(:,4)]; % 18 mug/L
 units.tWc_15_c4 = {'d', 'mg'};  label.tWc_15_c4 = {'time', 'dry weight'};
 temp.tWc_15_c4 = C2K(15.2);  units.temp.tWc_15_c4 = 'K';  label.temp.tWc_15_c4 = 'temperature';
 bibkey.tWc_15_c4 = 'Gaia2022';  comment.tWc_15_c4 = '18 mug/L';
data.tWc_15_c5 = [data.tWc_15(:,1), data.tWc_15(:,3)]; % 11 mug/L
 units.tWc_15_c5 = {'d', 'mg'};  label.tWc_15_c5 = {'time', 'dry weight'};
 temp.tWc_15_c5 = C2K(15.2);  units.temp.tWc_15_c5 = 'K';  label.temp.tWc_15_c5 = 'temperature';
 bibkey.tWc_15_c5 = 'Gaia2022';  comment.tWc_15_c5 = '11 mug/L';
data.tWc_15 = [data.tWc_15(:,1), data.tWc_15(:,2)]; % control
 units.tWc_15 = {'d', 'mg'};  label.tWc_15 = {'time', 'dry weight'};
 temp.tWc_15 = C2K(15.2);  units.temp.tWc_15 = 'K';  label.temp.tWc_15 = 'temperature';
 bibkey.tWc_15 = 'Gaia2022';  comment.tWc_15 = 'control treatment';
 
% t-W chronic data at 20 degC
data.tWc_20 = [ ... % time since hatch (d), dry weight (mg) in c_tmx = 0, 4, 11, 18, 25, 32 mug/L
8	0.394	0.471	0.407	0.176	0.348	0.343
10	0.701	0.884	0.277	0.636	0.794	0.525
11	0.559	1.073	0.986	0.759	0.645	0.299
12	0.830	1.076	0.881	0.867	0.645	0.267
13	0.903	1.123	1.071	0.733	0.788	0.335
15	1.319	1.409	0.798	0.980	0.826	NaN
18	1.747	1.194	1.094	0.995	1.101	NaN
   ];
data.tWc_20_c1 = [data.tWc_20(1:5,1), data.tWc_20(1:5,7)]; % 32 mug/L
 units.tWc_20_c1 = {'d', 'mg'};  label.tWc_20_c1 = {'time', 'dry weight'};
 temp.tWc_20_c1 = C2K(20.2);  units.temp.tWc_20_c1 = 'K';  label.temp.tWc_20_c1 = 'temperature';
 bibkey.tWc_20_c1 = 'Gaia2022';  comment.tWc_20_c1 = '32 mug/L';
data.tWc_20_c2 = [data.tWc_20(:,1), data.tWc_20(:,6)]; % 25 mug/L
 units.tWc_20_c2 = {'d', 'mg'};  label.tWc_20_c2 = {'time', 'dry weight'};
 temp.tWc_20_c2 = C2K(20.2);  units.temp.tWc_20_c2 = 'K';  label.temp.tWc_20_c2 = 'temperature';
 bibkey.tWc_20_c2 = 'Gaia2022';  comment.tWc_20_c2 = '25 mug/L';
data.tWc_20_c3 = [data.tWc_20(:,1), data.tWc_20(:,5)]; % 18 mug/L
 units.tWc_20_c3 = {'d', 'mg'};  label.tWc_20_c3 = {'time', 'dry weight'};
 temp.tWc_20_c3 = C2K(20.2);  units.temp.tWc_20_c3 = 'K';  label.temp.tWc_20_c3 = 'temperature';
 bibkey.tWc_20_c3 = 'Gaia2022';  comment.tWc_20_c3 = '18 mug/L';
data.tWc_20_c4 = [data.tWc_20(:,1), data.tWc_20(:,4)]; % 11 mug/L
 units.tWc_20_c4 = {'d', 'mg'};  label.tWc_20_c4 = {'time', 'dry weight'};
 temp.tWc_20_c4 = C2K(20.2);  units.temp.tWc_20_c4 = 'K';  label.temp.tWc_20_c4 = 'temperature';
 bibkey.tWc_20_c4 = 'Gaia2022';  comment.tWc_20_c4 = '11 mug/L';
data.tWc_20_c5 = [data.tWc_20(:,1), data.tWc_20(:,3)]; % 4 mug/L
 units.tWc_20_c5 = {'d', 'mg'};  label.tWc_20_c5 = {'time', 'dry weight'};
 temp.tWc_20_c5 = C2K(20.2);  units.temp.tWc_20_c5 = 'K';  label.temp.tWc_20_c5 = 'temperature';
 bibkey.tWc_20_c5 = 'Gaia2022';  comment.tWc_20_c5 = '4 mug/L';
data.tWc_20 = [data.tWc_20(:,1), data.tWc_20(:,2)]; % control
 units.tWc_20 = {'d', 'mg'};  label.tWc_20 = {'time', 'dry weight'};
 temp.tWc_20 = C2K(20.2);  units.temp.tWc_20 = 'K';  label.temp.tWc_20 = 'temperature';
 bibkey.tWc_20 = 'Gaia2022';  comment.tWc_20 = 'control treatment';
 
% t-W chronic data at 23 degC
data.tWc_23 = [ ... % time since hatch (d), dry weight (mg) in c_tmx = 0, 4, 11, 18, 25, 32 mug/L
4	0.038	0.014	0.039	0.021	0.024	NaN
6	0.123	0.093	0.124	0.044	0.159	0.130
8	0.267	0.486	0.420	0.184	0.145	NaN
10	0.903	1.060	0.813	0.256	NaN 	NaN
11	0.962	0.966	0.754	0.515	NaN 	NaN
13	1.131	1.056	0.518	0.223	NaN 	NaN
14	1.198	0.987	0.693	0.438	NaN 	NaN
15	1.361	1.221	0.604	NaN 	NaN 	NaN
   ];
data.tWc_23_c1 = [data.tWc_23(2,1), data.tWc_23(2,7)]; % 32 mug/L
 units.tWc_23_c1 = {'d', 'mg'};  label.tWc_23_c1 = {'time', 'dry weight'};
 temp.tWc_23_c1 = C2K(23.2);  units.temp.tWc_23_c1 = 'K';  label.temp.tWc_23_c1 = 'temperature';
 bibkey.tWc_23_c1 = 'Gaia2022';  comment.tWc_23_c1 = '32 mug/L';
data.tWc_23_c2 = [data.tWc_23(1:3,1), data.tWc_23(1:3,6)]; % 25 mug/L
 units.tWc_23_c2 = {'d', 'mg'};  label.tWc_23_c2 = {'time', 'dry weight'};
 temp.tWc_23_c2 = C2K(23.2);  units.temp.tWc_23_c2 = 'K';  label.temp.tWc_23_c2 = 'temperature';
 bibkey.tWc_23_c2 = 'Gaia2022';  comment.tWc_23_c2 = '25 mug/L';
data.tWc_23_c3 = [data.tWc_23(1:7,1), data.tWc_23(1:7,5)]; % 18 mug/L
 units.tWc_23_c3 = {'d', 'mg'};  label.tWc_23_c3 = {'time', 'dry weight'};
 temp.tWc_23_c3 = C2K(23.2);  units.temp.tWc_23_c3 = 'K';  label.temp.tWc_23_c3 = 'temperature';
 bibkey.tWc_23_c3 = 'Gaia2022';  comment.tWc_23_c3 = '18 mug/L';
data.tWc_23_c4 = [data.tWc_23(:,1), data.tWc_23(:,4)]; % 11 mug/L
 units.tWc_23_c4 = {'d', 'mg'};  label.tWc_23_c4 = {'time', 'dry weight'};
 temp.tWc_23_c4 = C2K(23.2);  units.temp.tWc_23_c4 = 'K';  label.temp.tWc_23_c4 = 'temperature';
 bibkey.tWc_23_c4 = 'Gaia2022';  comment.tWc_23_c4 = '11 mug/L';
data.tWc_23_c5 = [data.tWc_23(:,1), data.tWc_23(:,3)]; % 4 mug/L
 units.tWc_23_c5 = {'d', 'mg'};  label.tWc_23_c5 = {'time', 'dry weight'};
 temp.tWc_23_c5 = C2K(23.2);  units.temp.tWc_23_c5 = 'K';  label.temp.tWc_23_c5 = 'temperature';
 bibkey.tWc_23_c5 = 'Gaia2022';  comment.tWc_23_c5 = '4 mug/L';
data.tWc_23 = [data.tWc_23(:,1), data.tWc_23(:,2)]; % control
 units.tWc_23 = {'d', 'mg'};  label.tWc_23 = {'time', 'dry weight'};
 temp.tWc_23 = C2K(23.2);  units.temp.tWc_23 = 'K';  label.temp.tWc_23 = 'temperature';
 bibkey.tWc_23 = 'Gaia2022';  comment.tWc_23 = 'control treatment';

%--------------------------------------------------------------------------
%---------- End gaiac 2022 data -------------------------------------------
%--------------------------------------------------------------------------

%% set weights for all real data
weights = setweights(data, []);
% weight settings were determined experimentally after much estimations
% with different starting values

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
set8 = {'tSc_12', 'tSc_12_c1', 'tSc_12_c2', 'tSc_12_c3', 'tSc_12_c4', 'tSc_12_c5'}; comment8 = {'survival chron test at 12 degC'};
set9 = {'tSc_15', 'tSc_15_c1', 'tSc_15_c2', 'tSc_15_c3', 'tSc_15_c4', 'tSc_15_c5'}; comment9 = {'survival chron test at 15 degC'};
set10 = {'tSc_20', 'tSc_20_c1', 'tSc_20_c2', 'tSc_20_c3', 'tSc_20_c4', 'tSc_20_c5'}; comment10 = {'survival chron test at 20 degC'};
set11 = {'tSc_23', 'tSc_23_c1', 'tSc_23_c2', 'tSc_23_c3', 'tSc_23_c4', 'tSc_23_c5'}; comment11 = {'survival chron test at 23 degC'};
set12 = {'tWc_12', 'tWc_12_c1', 'tWc_12_c2', 'tWc_12_c3', 'tWc_12_c4', 'tWc_12_c5'}; comment12 = {'dry weight at 12 degC'};
set13 = {'tWc_15', 'tWc_15_c1', 'tWc_15_c2', 'tWc_15_c3', 'tWc_15_c4', 'tWc_15_c5'}; comment13 = {'dry weight at 15 degC'};
set14 = {'tWc_20', 'tWc_20_c1', 'tWc_20_c2', 'tWc_20_c3', 'tWc_20_c4', 'tWc_20_c5'}; comment14 = {'dry weight at 20 degC'};
set15 = {'tWc_23', 'tWc_23_c1', 'tWc_23_c2', 'tWc_23_c3', 'tWc_23_c4', 'tWc_23_c5'}; comment15 = {'dry weight at 23 degC'};
metaData.grp.sets = {set8,set9,set10,set11,set12,set13,set14,set15};
metaData.grp.comment = {comment8,comment9,comment10,comment11,comment12,comment13,comment14,comment15};

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
bibkey = 'Gaia2022'; type = 'Misc'; bib = [ ...
    'year = {2022}, ' ...
    'howpublished = {unpblished data}']; % Can be updated once the paper is published
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

