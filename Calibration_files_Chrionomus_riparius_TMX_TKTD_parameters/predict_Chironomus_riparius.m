function [prdData, info] = predict_Chironomus_riparius(par, data, auxData)

% unpack par, data, auxData
cPar = parscomp_st(par); vars_pull(par);
vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

global glo

filterChecks = E_Hp < E_Hb || ...
    z_b < 0 || h_b < 0 || z_s < 0;

if filterChecks
    prdData = []; info = 0; return
end

% compose vector of temperature parameters
pars_T = [T_A; T_H; T_AH]; 

% zero-variate data

% life cycle under default hax assumptions
v_Rj = kap/ (1 - kap) * E_Rj/ E_G; % scaled repro buffer density at pupation
pars_tj_hax = [g, k, v_Hb, v_Hp, v_Rj, v_He, kap, kap_V];
[t_j, t_e, t_p, t_b, l_j, l_e, l_p, l_b, l_i, rho_j, rho_B, u_Ee, info] = get_tj_hax(pars_tj_hax, f);
if ~info
    prdData = []; return
end

% birth
L_b = L_m * l_b; % cm, structural length at birth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% important modification to DEB hax model %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% at different f, we assume that pupation is not trigerred by a fixed repro
% buffer density (i.e. E_R/V) but a threshold of repro buffer by structural 
% length L (i.e. E_R/L)

L_j = L_m * l_j; % cm, structural length at pupation at f
E_Rj_L = E_Rj * L_j^2; % J/cm, repro buffer by length at pupation
% Here assumed to be independent of f
% E_Rj_L will be used in custom ODE system instead of E_Rj

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
%------------------- Start gaiac 2022 data predictions --------------------
%--------------------------------------------------------------------------

options = odeset;
par_LEHR = [p_Am, v, kap, kap_R, p_M, k_J, E_G, E_Hb, E_Hp, E_Rj_L, E_He, kap_V]; % pack parameters
E_b = p_Am/ v * L_b^3; % reserve at birth set at maximum reserve density (parent generation fed ad libitum)
LEHRDS_0 = [L_b; E_b; E_Hb; 0; 1; L_b; L_b; 0; 0; L_b]; % pack initial conditions
par_DS = [k_d, h_b, b_b, z_s, z_b, c_T, s_shrink, kap_G]; % pack tox parameters

%% emergence 12 degC
TC = tempcorr(temp.ae12, T_ref, pars_T); TC_kd = tempcorr(temp.ae12, T_ref, T_Akd);

k1 = -log(0.5)/(77.323 * exp(-0.072 * K2C(temp.ae12))); % # 1/d, k1
par_LEHRDS = [TC, 1, par_LEHR, par_DS, TC_kd]; % combine TC, f = 1, and parameters

glo.c_dyn = [1.25 0 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS, te, LEHRDSe] = odehax(LEHRDS_0, par_LEHRDS, 60); % ODE simulation with life events
prdData.ae12 = te(2);
prdData.Ni12 = LEHRDS(end,9);

glo.c_dyn = [1.25 19.59 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS, te, LEHRDSe] = odehax(LEHRDS_0, par_LEHRDS, 60); % ODE simulation with life events
prdData.Ni12_c5 = LEHRDS(end,9);

%% emergence 15 degC
TC = tempcorr(temp.ae15, T_ref, pars_T); TC_kd = tempcorr(temp.ae15, T_ref, T_Akd);

k1 = -log(0.5)/(77.323 * exp(-0.072 * K2C(temp.ae15))); % # 1/d, k1
par_LEHRDS = [TC, 1, par_LEHR, par_DS, TC_kd]; % combine TC, f = 1, and parameters

glo.c_dyn = [1.25 0 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS, te, LEHRDSe] = odehax(LEHRDS_0, par_LEHRDS, 90); % ODE simulation with life events
prdData.ae15 = te(2);
prdData.Ni15 = LEHRDS(end,9);

glo.c_dyn = [1.25 10.65 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS, te, LEHRDSe] = odehax(LEHRDS_0, par_LEHRDS, 90); % ODE simulation with life events
prdData.ae15_c5 = te(2);
prdData.Ni15_c5 = LEHRDS(end,9);

%% emergence 20 degC
TC = tempcorr(temp.ae20, T_ref, pars_T); TC_kd = tempcorr(temp.ae20, T_ref, T_Akd);

k1 = -log(0.5)/(77.323 * exp(-0.072 * K2C(temp.ae20))); % # 1/d, k1
par_LEHRDS = [TC, 1, par_LEHR, par_DS, TC_kd]; % combine TC, f = 1, and parameters

glo.c_dyn = [1.25 0 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS, te, LEHRDSe] = odehax(LEHRDS_0, par_LEHRDS, 60); % ODE simulation with life events
prdData.ae20 = te(2);
prdData.Ni20 = LEHRDS(end,9);

glo.c_dyn = [1.25 12.37 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS, te, LEHRDSe] = odehax(LEHRDS_0, par_LEHRDS, 60); % ODE simulation with life events
prdData.ae20_c4 = te(2);
prdData.Ni20_c4 = LEHRDS(end,9);

glo.c_dyn = [1.25 4.35 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS, te, LEHRDSe] = odehax(LEHRDS_0, par_LEHRDS, 60); % ODE simulation with life events
prdData.ae20_c5 = te(2);
prdData.Ni20_c5 = LEHRDS(end,9);

%% emergence 23 degC
TC = tempcorr(temp.ae23, T_ref, pars_T); TC_kd = tempcorr(temp.ae23, T_ref, T_Akd);

k1 = -log(0.5)/(77.323 * exp(-0.072 * K2C(temp.ae23))); % # 1/d, k1
par_LEHRDS = [TC, 1, par_LEHR, par_DS, TC_kd]; % combine TC, f = 1, and parameters

glo.c_dyn = [1.25 0 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS, te, LEHRDSe] = odehax(LEHRDS_0, par_LEHRDS, 60); % ODE simulation with life events
prdData.ae23 = te(2);
prdData.Ni23 = LEHRDS(end,9);

glo.c_dyn = [1.25 6.49 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS, te, LEHRDSe] = odehax(LEHRDS_0, par_LEHRDS, 60); % ODE simulation with life events
prdData.ae23_c5 = te(2);
prdData.Ni23_c5 = LEHRDS(end,9);

%% time - survival and dry weight; chronic tests at 12 degC
TC = tempcorr(temp.tSc_12, T_ref, pars_T); TC_kd = tempcorr(temp.tSc_12, T_ref, T_Akd);

par_LEHRDS = [TC, 1, par_LEHR, par_DS, TC_kd]; % combine TC, f = 1, and parameters
k1 = -log(0.5)/(77.323 * exp(-0.072 * K2C(temp.tSc_12))); % # 1/d, k1

glo.c_dyn = [1.25 0 k1]; % [t0(d) c(mug/L) k1(1/d)]
% [t, LEHRDS, te, LEHRDSe] = odehax(LEHRDS_0, par_LEHRDS, 60); % ODE simulation with life events
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc12 = (interp1(t(:), LEHRDS(:,5), tSc_12(:,1)+0.25)); % -, survival
prdData.tSc_12 = Sc12;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_12(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_12(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_12(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc12 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_12 = Wc12;

glo.c_dyn = [1.25 50.63 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc12c1 = (interp1(t(:), LEHRDS(:,5), tSc_12_c1(:,1)+0.25)); % -, survival
prdData.tSc_12_c1 = Sc12c1;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_12_c1(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_12_c1(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_12_c1(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc12c1 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_12_c1 = Wc12c1;

glo.c_dyn = [1.25 43.89 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc12c2 = (interp1(t(:), LEHRDS(:,5), tSc_12_c2(:,1)+0.25)); % -, survival
prdData.tSc_12_c2 = Sc12c2;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_12_c2(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_12_c2(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_12_c2(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc12c2 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_12_c2 = Wc12c2;

glo.c_dyn = [1.25 36.44 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc12c3 = (interp1(t(:), LEHRDS(:,5), tSc_12_c3(:,1)+0.25)); % -, survival
prdData.tSc_12_c3 = Sc12c3;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_12_c3(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_12_c3(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_12_c3(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc12c3 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_12_c3 = Wc12c3;

glo.c_dyn = [1.25 28.63 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc12c4 = (interp1(t(:), LEHRDS(:,5), tSc_12_c4(:,1)+0.25)); % -, survival
prdData.tSc_12_c4 = Sc12c4;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_12_c4(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_12_c4(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_12_c4(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc12c4 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_12_c4 = Wc12c4;

glo.c_dyn = [1.25 19.59 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc12c5 = (interp1(t(:), LEHRDS(:,5), tSc_12_c5(:,1)+0.25)); % -, survival
prdData.tSc_12_c5 = Sc12c5;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_12_c5(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_12_c5(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_12_c5(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc12c5 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_12_c5 = Wc12c5;

%% time - survival and dry weight; chronic tests at 15 degC
TC = tempcorr(temp.tSc_15, T_ref, pars_T); TC_kd = tempcorr(temp.tSc_15, T_ref, T_Akd);

par_LEHRDS = [TC, 1, par_LEHR, par_DS, TC_kd]; % combine TC, f = 1, and parameters
k1 = -log(0.5)/(77.323 * exp(-0.072 * K2C(temp.tSc_15))); % # 1/d, k1

glo.c_dyn = [1.25 0 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc15 = (interp1(t(:), LEHRDS(:,5), tSc_15(:,1)+0.25)); % -, survival
prdData.tSc_15 = Sc15;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_15(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_15(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_15(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc15 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_15 = Wc15;

glo.c_dyn = [1.25 38.14 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc15c1 = (interp1(t(:), LEHRDS(:,5), tSc_15_c1(:,1)+0.25)); % -, survival
prdData.tSc_15_c1 = Sc15c1;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_15_c1(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_15_c1(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_15_c1(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc15c1 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_15_c1 = Wc15c1;

glo.c_dyn = [1.25 29.68 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc15c2 = (interp1(t(:), LEHRDS(:,5), tSc_15_c2(:,1)+0.25)); % -, survival
prdData.tSc_15_c2 = Sc15c2;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_15_c2(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_15_c2(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_15_c2(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc15c2 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_15_c2 = Wc15c2;

glo.c_dyn = [1.25 24.69 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc15c3 = (interp1(t(:), LEHRDS(:,5), tSc_15_c3(:,1)+0.25)); % -, survival
prdData.tSc_15_c3 = Sc15c3;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_15_c3(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_15_c3(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_15_c3(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc15c3 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_15_c3 = Wc15c3;

glo.c_dyn = [1.25 18.28 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc15c4 = (interp1(t(:), LEHRDS(:,5), tSc_15_c4(:,1)+0.25)); % -, survival
prdData.tSc_15_c4 = Sc15c4;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_15_c4(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_15_c4(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_15_c4(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc15c4 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_15_c4 = Wc15c4;

glo.c_dyn = [1.25 10.65 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc15c5 = (interp1(t(:), LEHRDS(:,5), tSc_15_c5(:,1)+0.25)); % -, survival
prdData.tSc_15_c5 = Sc15c5;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_15_c5(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_15_c5(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_15_c5(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc15c5 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_15_c5 = Wc15c5;

%% time - survival and dry weight; chronic tests at 20 degC
TC = tempcorr(temp.tSc_20, T_ref, pars_T); TC_kd = tempcorr(temp.tSc_20, T_ref, T_Akd);

par_LEHRDS = [TC, 1, par_LEHR, par_DS, TC_kd]; % combine TC, f = 1, and parameters
k1 = -log(0.5)/(77.323 * exp(-0.072 * K2C(temp.tSc_20))); % # 1/d, k1

glo.c_dyn = [1.25 0 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc20 = (interp1(t(:), LEHRDS(:,5), tSc_20(:,1)+0.25)); % -, survival
prdData.tSc_20 = Sc20;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_20(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_20(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_20(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc20 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_20 = Wc20;

glo.c_dyn = [1.25 36.62 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc20c1 = (interp1(t(:), LEHRDS(:,5), tSc_20_c1(:,1)+0.25)); % -, survival
prdData.tSc_20_c1 = Sc20c1;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_20_c1(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_20_c1(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_20_c1(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc20c1 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_20_c1 = Wc20c1;

glo.c_dyn = [1.25 28.10 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc20c2 = (interp1(t(:), LEHRDS(:,5), tSc_20_c2(:,1)+0.25)); % -, survival
prdData.tSc_20_c2 = Sc20c2;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_20_c2(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_20_c2(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_20_c2(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc20c2 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_20_c2 = Wc20c2;

glo.c_dyn = [1.25 20.47 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc20c3 = (interp1(t(:), LEHRDS(:,5), tSc_20_c3(:,1)+0.25)); % -, survival
prdData.tSc_20_c3 = Sc20c3;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_20_c3(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_20_c3(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_20_c3(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc20c3 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_20_c3 = Wc20c3;

glo.c_dyn = [1.25 12.37 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc20c4 = (interp1(t(:), LEHRDS(:,5), tSc_20_c4(:,1)+0.25)); % -, survival
prdData.tSc_20_c4 = Sc20c4;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_20_c4(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_20_c4(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_20_c4(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc20c4 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_20_c4 = Wc20c4;

glo.c_dyn = [1.25 4.35 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc20c5 = (interp1(t(:), LEHRDS(:,5), tSc_20_c5(:,1)+0.25)); % -, survival
prdData.tSc_20_c5 = Sc20c5;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_20_c5(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_20_c5(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_20_c5(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc20c5 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_20_c5 = Wc20c5;

%% time - survival and dry weight; chronic tests at 23 degC
TC = tempcorr(temp.tSc_23, T_ref, pars_T); TC_kd = tempcorr(temp.tSc_23, T_ref, T_Akd);

par_LEHRDS = [TC, 1, par_LEHR, par_DS, TC_kd]; % combine TC, f = 1, and parameters
k1 = -log(0.5)/(77.323 * exp(-0.072 * K2C(temp.tSc_23))); % # 1/d, k1

glo.c_dyn = [1.25 0 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc23 = (interp1(t(:), LEHRDS(:,5), tSc_23(:,1)+0.25)); % -, survival
prdData.tSc_23 = Sc23;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_23(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_23(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_23(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc23 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_23 = Wc23;

glo.c_dyn = [1.25 50.28 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc23c1 = (interp1(t(:), LEHRDS(:,5), tSc_23_c1(:,1)+0.25)); % -, survival
prdData.tSc_23_c1 = Sc23c1;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_23_c1(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_23_c1(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_23_c1(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc23c1 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_23_c1 = Wc23c1;

glo.c_dyn = [1.25 39.81 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc23c2 = (interp1(t(:), LEHRDS(:,5), tSc_23_c2(:,1)+0.25)); % -, survival
prdData.tSc_23_c2 = Sc23c2;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_23_c2(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_23_c2(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_23_c2(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc23c2 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_23_c2 = Wc23c2;

glo.c_dyn = [1.25 28.99 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc23c3 = (interp1(t(:), LEHRDS(:,5), tSc_23_c3(:,1)+0.25)); % -, survival
prdData.tSc_23_c3 = Sc23c3;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_23_c3(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_23_c3(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_23_c3(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc23c3 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_23_c3 = Wc23c3;

glo.c_dyn = [1.25 17.99 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc23c4 = (interp1(t(:), LEHRDS(:,5), tSc_23_c4(:,1)+0.25)); % -, survival
prdData.tSc_23_c4 = Sc23c4;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_23_c4(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_23_c4(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_23_c4(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc23c4 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_23_c4 = Wc23c4;

glo.c_dyn = [1.25 6.49 k1]; % [t0(d) c(mug/L) k1(1/d)]
[t, LEHRDS] = ode45(@get_LEHRDS, [0 60], LEHRDS_0, options, par_LEHRDS); % ODE simulation w/o life events
Sc23c5 = (interp1(t(:), LEHRDS(:,5), tSc_23_c5(:,1)+0.25)); % -, survival
prdData.tSc_23_c5 = Sc23c5;
WL3 = (interp1(t(:), LEHRDS(:,1), tWc_23_c5(:,1)+0.25)).^3 .* d_V; % g, dry weight of structure
WE = (interp1(t(:), LEHRDS(:,2), tWc_23_c5(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reserve
WR = (interp1(t(:), LEHRDS(:,4), tWc_23_c5(:,1)+0.25)).* w_E./ mu_E; % g, dry weight of reproduction buffer
Wc23c5 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
prdData.tWc_23_c5 = Wc23c5;

%--------------------------------------------------------------------------
%-------------------- End gaiac 2022 data predictions ---------------------
%--------------------------------------------------------------------------

end

%%
% Overarching function to merge consecutive ode45-runs of the different hax life stages
function [t, LEHRDS, te, LEHRDSe] = odehax(LEHRDS_0, par_LEHRDS, tfinal)
tstart = 0;
refine = 4;

p_Am  = par_LEHRDS(3); % J / d.cm^2, Surface-area specific maximum assimilation rate
v     = par_LEHRDS(4); % cm/d, Energy conductance
kap   = par_LEHRDS(5); % -, Allocation fraction to soma
kap_R = par_LEHRDS(6); % -, Reproduction efficiency
p_M   = par_LEHRDS(7); % J/d.cm^3, Vol-spec somatic maint
k_J   = par_LEHRDS(8); % 1/d, Maturity maint rate coefficient
E_G   = par_LEHRDS(9); % J/cm^3, Spec cost for structure
E_Hb  = par_LEHRDS(10); % J, Maturity at birth
kap_V = par_LEHRDS(14); % -, Conversion efficiency from larval reserve to larval structure, back to imago reserve
k_M   = p_M / E_G;         % 1/d, somatic maintenance rate coefficient
k     = k_J / k_M;         % -, maintenance ratio
E_m   = p_Am / v;          % J/cm^3, [E_m], reserve capacity
g     = E_G / kap / E_m ;  % -, energy investment ratio
U_Hb  = E_Hb / p_Am;       % d.cm^2, Scaled maturity at birth male (d.cm^2)
V_Hb  = U_Hb / (1 - kap); v_Hb = V_Hb * g^2 * k_M^3 / v^2;
L_m   = v/ k_M/ g;         % cm, maximum length

% simulate until pupa stage
options = odeset('Events', @pupationEvent, 'Refine', refine);
[t, LEHRDS, te, LEHRDSe] = ode45(@get_LEHRDS, [tstart tfinal], LEHRDS_0, options, par_LEHRDS); % ODE simulation
% simulate pupa until emergence
if (t(end) < tfinal)
    LEHRDSet = LEHRDSe; % initial states for pupa stage
%     e_b  = LEHRDSet(2) * v / (LEHRDSet(1)^3 * p_Am); % -, scaled reserve density at birth of new eggs
    e_b = 1; % Assuming no maternal effects related to pMoA
    u_E0 = get_ue0([g, k, v_Hb], e_b); % -, scaled initial reserve of new eggs
    E_0  = u_E0 * g * E_m * L_m^3;   % J, initial energy content of new eggs
    LEHRDSet(2) = LEHRDSe(2) + LEHRDSe(1)^3 * E_G * kap_V; % Reserve expanded by transformed structure
    LEHRDSet(1) = 0.000001; LEHRDSet(6) = LEHRDSet(1); LEHRDSet(7) = LEHRDSet(1); LEHRDSet(10) = LEHRDSet(1); % Structure to zero
    LEHRDSet(3) = 0; % Maturity to zero
    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation. 'refine' is 4 by default.
    nt = length(t);
    options = odeset('Events', @emergenceEvent, 'Refine', refine, 'InitialStep', t(nt)-t(nt-refine), 'MaxStep',t(nt)-t(1));
    [t2, LEHRDS2, te2, LEHRDSe2] = ode45(@get_LEHRDS, [te tfinal], LEHRDSet, options, par_LEHRDS); % ODE simulation
    % Accumulate output
    nt = length(t2);
    t = [t; t2(2:nt)];
    LEHRDS = [LEHRDS; LEHRDS2(2:nt,:)];
    te = [te; te2]; % time at event
    LEHRDSe = [LEHRDSe; LEHRDSe2]; % states at event
end
% simulate imago after emergence
if (t(end) < tfinal)
    LEHRDSe2t = LEHRDSe2; % initial states for pupa stage
    LEHRDSe2t(9) = LEHRDSe2(4) * kap_R / E_0; % Make eggs from repro buffer
    LEHRDSe2t(4) = 0; % Repro buffer to zero
    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation. 'refine' is 4 by default.
    nt = length(t2);
    options = odeset('Refine', refine, 'InitialStep', t2(nt)-t2(nt-refine),'MaxStep',t2(nt)-t2(1));
    [t3, LEHRDS3] = ode45(@get_LEHRDS, [te2 tfinal], LEHRDSe2t, options, par_LEHRDS); % ODE simulation
    % Accumulate output
    nt = length(t3);
    t = [t; t3(2:nt)];
    LEHRDS = [LEHRDS; LEHRDS3(2:nt,:)];
end
end

%%
% Event function to determine the moment of pupation
function [value, isterminal, direction] = pupationEvent(t, LEHRDS, par_LEHRDS)
L       = LEHRDS(1); % cm, Structural length
E_R     = LEHRDS(4); % # J, Cumulative energy in reproduction buffer
E_Rj_L = par_LEHRDS(12); % J/cm, Energy in repro buffer by length at pupation
% E_Rj    = par_LEHRDS(12); % J/cm^3, Reproduction buffer density at pupation

value = E_R - (E_Rj_L * L); % The value that we want to be zero
% value = E_R - (E_Rj * L^3); % The value that we want to be zero
isterminal = 1;  % Halt integration
direction = 0;   % The zero can be approached from either direction
end

%%
% Event function to determine the moment of emergence
function [value, isterminal, direction] = emergenceEvent(t, LEHRDS, par_LEHRDS)
E_H     = LEHRDS(3); % J, Energy invested in maturity
E_R     = LEHRDS(4); % J, Cumulative energy in reproduction buffer
E_He    = par_LEHRDS(13); % J, Maturity at emergence

value = (E_H - E_He) + (E_R == 0) ; % The value that we want to be zero
isterminal = 1;  % Halt integration
direction = 0;   % The zero can be approached from either direction
end

%%
% ODE system of the hax model; different life-stages are covered by if-clauses
% The system is designed to work with the event functions for pupa and emergence
% If executed without these, it should behave like an abj model w/o egg laying
function [dLEHRDS] = get_LEHRDS(t, LEHRDS, par_LEHRDS)

% Read parameters
global glo
TC    = par_LEHRDS(1); % -, Temperature correction
f     = par_LEHRDS(2); % -, Functional response
p_Am  = par_LEHRDS(3); % J / d.cm^2, Surface-area specific maximum assimilation rate
v     = par_LEHRDS(4); % cm/d, Energy conductance
kap   = par_LEHRDS(5); % -, Allocation fraction to soma
% kap_R = par_LEHRDS(6); % -, Reproduction efficiency
p_M   = par_LEHRDS(7); % J/d.cm^3, Vol-spec somatic maint
k_J   = par_LEHRDS(8); % 1/d, Maturity maint rate coefficient
E_G   = par_LEHRDS(9); % J/cm^3, Spec cost for structure
E_Hb  = par_LEHRDS(10); % J, Maturity at birth
E_Hp  = par_LEHRDS(11); % J, Maturity at puberty
E_Rj_L = par_LEHRDS(12); % J/cm, Energy in repro buffer by length at pupation
E_He  = par_LEHRDS(13); % J, Maturity at emergence
% kap_V  = par_LEHRDS(14); % -, Conversion efficiency from larval reserve to larval structure, back to imago reserve
k_d   = par_LEHRDS(15); % 1/d, Dominant rate constant
h_b   = par_LEHRDS(16); % mg/L, Background mortality rate (SD)
b_b   = par_LEHRDS(17); % L/(mg d), Killing rate (SD)
z_s   = par_LEHRDS(18); % mg/L, Population-wide tolerance threshold (SD)
z_b   = par_LEHRDS(19); % mg/L, mug/L, Population-wide tolerance threshold (sublethal)
c_T   = par_LEHRDS(20); % mg/L, Tolerance concentration (sublethal)
s_shrink = par_LEHRDS(21); % -, Shrinking stress coefficient
kap_G = par_LEHRDS(22); % -, Growth efficiency
TC_kd = par_LEHRDS(23); % -, Temperature correction for dominant rate constant

% initialize state variables
L       = LEHRDS(1); % cm, Structural length
E       = LEHRDS(2); % J, Energy in reserve
E_H     = LEHRDS(3); % # J, Energy invested in maturity
E_R     = LEHRDS(4); % # J, Cumulative energy in reproduction buffer
S       = LEHRDS(5); % -, Surviving fraction
Lsa     = LEHRDS(6); % cm, Min(Current Length, Length at start of acceleration)
Lea     = LEHRDS(7); % cm, Min(Current Length, Length at end of acceleration)
Dmg     = LEHRDS(8); % mg/L, Damage state
eggs    = LEHRDS(9); % #, Produced eggs
maxL    = LEHRDS(10); % cm, Maximum length reached so far

% Metabolic acceleration (Type M)
M = Lea / Lsa; % Acceleration factor
p_Am_M = p_Am * M; % Metabolic acceleration affects p_Am
v_M = v * M; % Metabolic acceleration affects v

% Scaled reserve density (to check for starvation)
e = E * v / (L^3 * p_Am);
% Scaled structural length (to check for starvation)
l = L * p_M / (kap * p_Am);

% Toxicodynamics (DEB-TKTD)
s = (1 / c_T) * max(0, Dmg - z_b); % Stress factor DEB-TKTD

% Mode of action: effect on feeding
f = f * max(0, 1 - s);
% E_G = E_G * max(0, 1 + s);

% Assimilation flux (from birth until pupation)
if (E_H < E_Hb) && (E_R == 0) ||... % Before birth
        (E_R > (E_Rj_L * L)) && (E_H < E_He) ||... % During pupa stage
        (eggs > 0) % After emergence
    p_A = 0; % J/d increase in maturation
else
    p_A = f * p_Am_M * L^2; % DEB Book 3rd ed. Eq. 2.3
end

% Somatic maintenance flux (always)
p_S = p_M * L^3;

% Growth rate (from embryo until emergence unless starving)
if (eggs == 0) && (e > (l/M))
    r = (E / L^3 * v_M / L - (p_S / L^3 / kap)) / (E / L^3 + E_G / kap);
elseif (L <= 0.8 * maxL) % threshold for shrinking
    r = 0;
else
    r = (E / L^3 * v_M / L - (p_S / L^3 / kap)) / (E / L^3 + E_G * kap_G^2 / kap); % shrinking
end

% Maturity maintenance flux (always)
p_J = k_J * E_H; % DEB Book 3rd ed. Eq. 2.19

% Mobilization flux (always; changes to only maintenance when emerged)
if (eggs == 0)
    p_C = E * (v_M / L - r); % DEB Book 3rd ed. Eq. 2.12
else
    p_C = p_S + p_J; % DEB Book 3rd ed. Eq. 2.12
end

% Reproduction flux (from embryo until emergence)
if (eggs == 0) && (e > (l/M))
    p_R = (1 - kap) * p_C - p_J; % DEB Book 3rd ed. Eq. 2.18
else
    p_R = max(0, p_C - p_S - p_J); % No kappa rule if growth is zero
end

% Toxicodynamics (GUTS-RED-SD)
L_m = kap * p_Am / p_M; % cm, maximum structural length
h = b_b * max(0, Dmg - z_s / (L/L_m)); % calculate the hazard rate

% Changes in state variables
dL = (L * r / 3) * TC; % DEB Book 3rd ed. Chap. 2.6.1
dE = (p_A - p_C) * TC; % DEB Book 3rd ed. Chap. 2.3

if (E_H < E_Hp) || (E_R > (E_Rj_L * L)) && (E_H < E_He)
    dE_H = p_R * TC;
    dE_R = 0;
else
    dE_H = 0;
    dE_R = p_R * TC;
end

if (E_H < E_Hb) && (E_R < (E_Rj_L * L))
    dLsa = dL;
else
    dLsa = 0;
end

if (E_H < E_Hp) && (E_R < (E_Rj_L * L))
    dLea = dL;
else
    dLea = 0;
end

% Toxicokinetics
c_w0 = glo.c_dyn(2); % mug/L, initial environmental concentration
k1 = glo.c_dyn(3); % 1/d, k1
if (t < glo.c_dyn(1)) % time between birth and start of exposure
    c_w = 0;
else % after first exposure
    c_w = c_w0 * exp(-k1 * (t - glo.c_dyn(1)));
end
dDmg = k_d * (c_w - Dmg) * TC_kd;

deggs = 0 * TC;

% shrinking stress
dmaxL = max(0, dL); % cm, max stuctural length
k_M = p_M/ E_G;  % 1/d, maturity maintenance rate coefficient
h_sh = s_shrink * k_M * max(0, maxL - L)/maxL; % hazard due to shrinking

dS = -(h + h_b + h_sh) * S ; % with tox
% dS = -(h_b + h_sh) * S ; % without tox

dLEHRDS = [dL; dE; dE_H; dE_R; dS; dLsa; dLea; dDmg; deggs; dmaxL]; % collect derivatives
end

