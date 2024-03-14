function [prdData, info] = predict_Chironomus_riparius(par, data, auxData)

% unpack par, data, auxData
cPar = parscomp_st(par); vars_pull(par);
vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

filterChecks = z_m > z || E_Hp < E_Hb || ...
    ~reach_birth(g, k, v_Hb, f1) || f1 > f2 ||...
    ~reach_birth(g, k, v_Hb, f2) || f2 > f3 ||...
    ~reach_birth(g, k, v_Hb, f3) || f3 > f4 ||...
    ~reach_birth(g, k, v_Hb, f4) || f4 > 1 ||...
    T_H <= T_ref || T_L >= T_ref;

if filterChecks
    prdData = []; info = 0; return
end

% compute temperature correction factors
% pars_T = [T_A; T_L; T_H; T_AL; T_AH];
pars_T = [T_A; T_H; T_AH];

TC_ab      = tempcorr(C2K(21), T_ref, pars_T);
TC_am      = tempcorr(C2K(26), T_ref, pars_T);
TC_te      = tempcorr(C2K(26), T_ref, pars_T);
TC         = tempcorr(C2K(21), T_ref, pars_T); % all experiments in PeryMons2002 were at 21 deg C
TC15       = tempcorr(C2K(15), T_ref, pars_T);
TC196      = tempcorr(C2K(19.6), T_ref, pars_T);
TC21       = tempcorr(C2K(21), T_ref, pars_T);
TC244      = tempcorr(C2K(24.4), T_ref, pars_T);
TC267      = tempcorr(C2K(26.7), T_ref, pars_T);
TC_20      = tempcorr(C2K(20), T_ref, pars_T);
TC_10      = tempcorr(C2K(10), T_ref, pars_T);

% zero-variate data

% life cycle
v_Rj = kap/ (1 - kap) * E_Rj/ E_G;
pars_tj_hax = [g, k, v_Hb, v_Hp, v_Rj, v_He, kap, kap_V];
[t_j, t_e, t_p, t_b, l_j, l_e, l_p, l_b, l_i, rho_j, rho_B, u_Ee, info] = get_tj_hax(pars_tj_hax, f);
if ~info
    prdData = []; return
end

% initial
u_E0 = get_ue0([g, k, v_Hb], f);    % -, scaled initial reserve
E_0 = u_E0 * g * E_m * L_m^3;       % J, initial energy content
Wd_0 = E_0/ mu_E * w_E;             % g, initial dry weight

% birth
aT_b = t_b/ k_M/ TC_ab;             % d, age at birth at f and T
L_b = L_m * l_b;                    % cm, structural length at birth
Lw_b = L_b/ del_M; % cm, length at birth

% pupation
L_j = L_m * l_j;                    % cm, structural length at pupation
Lw_j = L_j/ del_M;                  % cm, physical length at pupation
Ww_j = L_j^3 * (1 + f * w);         % g, wet weight at pupation, excluding reprod buffer at pupation
Ww_Rj = E_Rj * L_j^3 * w_E/ mu_E/ d_E; % g, wet weight repod buffer at pupation
Ww_j = Ww_j + Ww_Rj;                % g, wet weight including reprod buffer

% emergence
L_e = L_m * l_e;                    % cm, structural length at emergence
Wd_e = d_V * L_e^3 + d_E * Ww_Rj + u_Ee * g * E_m * L_m^3 * w_E/ mu_E; % g, dry weight including reprod buffer
tT_e = (t_e - t_j)/ k_M/ TC_te;            % d, time since pupation at emergence

% life span
pars_tm = [g; l_T; h_a/ k_M^2; s_G]; % compose parameter vector at T_ref
t_m = get_tm_s(pars_tm, f, l_b);     % -, scaled mean life span at T_ref
aT_m = (t_m-t_e)/ k_M/ TC_am;        % d, mean life span as imago

% males
p_Am_m = z_m * p_M/ kap;             % J/d.cm^2, {p_Am} spec assimilation flux
E_m_m = p_Am_m/ v;                   % J/cm^3, reserve capacity [E_m]
g_m = E_G/ (kap* E_m_m);             % -, energy investment ratio
m_Em_m = y_E_V * E_m_m/ E_G;         % mol/mol, reserve capacity
w_m = m_Em_m * w_E/ w_V;             % -, contribution of reserve to weight
L_m_m = kap * p_Am_m/ p_M;           % cm, max stuct length for males
U_Hb_m = E_Hb/ p_Am_m; % d.cm^2, Scaled maturity at birth male (d.cm^2)
U_Hp_m = E_Hp/ p_Am_m; % d.cm^2, Scaled maturity at metam male (d.cm^2)
U_He_m = E_He/ p_Am_m; % d.cm^2, Scaled maturity at puberty male (d.cm^2)
V_Hb_m = U_Hb_m / (1 - kap);  v_Hb_m = V_Hb_m * g_m^2 * k_M^3 / v^2;
V_Hp_m = U_Hp_m / (1 - kap);  v_Hp_m = V_Hp_m * g_m^2 * k_M^3 / v^2;
V_He_m = U_He_m / (1 - kap);  v_He_m = V_He_m * g_m^2 * k_M^3 / v^2;
% life cycle males
pars_tj_hax_m = [g_m, k, v_Hb_m, v_Hp_m, v_Rj, v_He_m, kap, kap_V];
[t_j_m, t_e_m, t_p_m, t_b_m, l_j_m, l_e_m, l_p_m, l_b_m, l_i_m, rho_j_m, rho_B_m, u_Ee_m, info_m] = get_tj_hax(pars_tj_hax_m, f);
L_b_m = L_m_m * l_b_m;

% pack to output
prdData.ab = aT_b;
prdData.te = tT_e;
prdData.am = aT_m;
prdData.Lb = Lw_b;
prdData.Lj = Lw_j;
prdData.Wd0 = Wd_0;
prdData.Wwj = Ww_j;
prdData.Wde = Wd_e;

%--------------------------------------------------------------------------
%-------------- Start zero var gaiac 2022 data predictions ----------------
%--------------------------------------------------------------------------

aT_e = (t_e - t_b)/ k_M/ TC_20; % d, age at emergence
prdData.ae = aT_e;

% life cycle with f = 1
pars_tj_hax = [g, k, v_Hb, v_Hp, v_Rj, v_He, kap, kap_V];
[t_j, t_e, t_p, t_b, l_j, l_e, l_p, l_b, l_i, rho_j, rho_B, u_Ee, info] = get_tj_hax(pars_tj_hax, f);

ae_12 = (t_e - t_b)./ k_M./ (tempcorr(C2K(11.8), T_ref, pars_T)); prdData.ae12 = ae_12; % d, time since birth at emergence
ae_15 = (t_e - t_b)./ k_M./ (tempcorr(C2K(15.2), T_ref, pars_T)); prdData.ae15 = ae_15; % d, time since birth at emergence
ae_20 = (t_e - t_b)./ k_M./ (tempcorr(C2K(20.2), T_ref, pars_T)); prdData.ae20 = ae_20; % d, time since birth at emergence
ae_23 = (t_e - t_b)./ k_M./ (tempcorr(C2K(23.2), T_ref, pars_T)); prdData.ae23 = ae_23; % d, time since birth at emergence
L_j = L_m * l_j; % cm, structural length at pupation
u_E0 = get_ue0([g, k, v_Hb], f); % -, scaled initial reserve
E_0 = u_E0 * g * E_m * L_m^3; % J, initial energy content
N_i = kap_R * E_Rj * L_j^3/ E_0; % #, ultimate reproduction
prdData.N12 = N_i;
prdData.N15 = N_i;
prdData.N20 = N_i;

%--------------------------------------------------------------------------
%--------------- End zero var gaiac 2022 data predictions -----------------
%--------------------------------------------------------------------------

% uni-variate data

% PeryMons: 21 deg C, different food levels

% Length-dry weight L-Wd
% it includes the dry weight of the eggs
pars_R = [kap; kap_R; g; k_J * TC; k_M * TC; 0; v * TC; U_Hb/ TC; U_Hp/ TC; U_Hp/TC + 1e-8];

L_p = l_p * L_m; L_i = l_i * L_m; % cm, structural lengths at puberty and ultimate at f
r_B = rho_B * k_M * TC; % 1/d, von Bert growth rate at f and T
a_b = t_b/ k_M/ TC; % d, age at birth at T and f
a_p = t_p/ k_M/ TC; % d, age at puberty at f and T
L    = LW(:,1) * del_M; % cm, structural length
t = 1/ r_B * log((L_i - L_p)./ (L_i - L)); % d, time since puberty
[N, Lpred, U_E0] = cum_reprod_j(t + (a_p - a_b), f, pars_R); % cum reproductive output at T
E_0  = U_E0 * p_Am * TC; % J, energy in the egg
Wd_0 = w_E/ mu_E * E_0; % g, dry mass of the egg
EWd = L.^3 * d_V * (1 + f * w); % g, dry weight  excluding the reproduction buffer
EWd = (EWd +  N.* Wd_0) * 1e3; % mg, dry weight including the reproduction buffer

% These predictions use the default hax assumption for the pupation cue and its food dependency:
[t_j, t_e, t_p, t_b, l_j, l_e, l_p, l_b, l_i, rho_j, rho_B, u_Ee, info] = get_tj_hax(pars_tj_hax, f2);
L_j = l_j * L_m; % cm,
EN2   = kap_R * E_Rj * L_j^3/ E_0;             % #/d, ultimate reproduction rate at T
[t_j, t_e, t_p, t_b, l_j, l_e, l_p, l_b, l_i, rho_j, rho_B, u_Ee, info] = get_tj_hax(pars_tj_hax, f3);
L_j = l_j * L_m; % cm,
EN3   = kap_R * E_Rj * L_j^3/ E_0;             % #/d, ultimate reproduction rate at T
[t_j, t_e, t_p, t_b, l_j, l_e, l_p, l_b, l_i, rho_j, rho_B, u_Ee, info] = get_tj_hax(pars_tj_hax, f4);
L_j = l_j * L_m; % cm,
EN4   = kap_R * E_Rj * L_j^3/ E_0;             % #/d, ultimate reproduction rate at T

prdData.N2 = EN2;
prdData.N3 = EN3;
prdData.N4 = EN4;

% time-length and time-weight data for larvae at f 
% redone by Koch in 2022 to also include the reproduction buffer in weight data
% using a more complex ODE system that has been thoroughly tested

options = odeset;
par_CS = [s_shrink, kap_G, h_b_shrink]; % pack survival parameters
E_b = p_Am/ v * L_b^3; % reserve at birth set at maximum reserve density (parent generation fed ad libitum)
E_b_m = p_Am_m/ v * L_b_m^3; % reserve at birth set at maximum reserve density (parent generation fed ad libitum)

% males
par_LEHR = [p_Am_m, v, kap, kap_R, p_M, k_J, E_G, E_Hb, E_Hp, E_Rj, E_He, kap_V]; % pack parameters, male
LEHRS_0 = [L_b_m; E_b_m; E_Hb; 0; 100; L_b_m; L_b_m; 0; L_b_m]; % pack initial conditions
par_LEHRS = [TC, f, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL_m = (interp1(t(:), LEHRS(:,1), tL(:,1)))/del_M;

% females
par_LEHR = [p_Am, v, kap, kap_R, p_M, k_J, E_G, E_Hb, E_Hp, E_Rj, E_He, kap_V]; % pack parameters, female
LEHRS_0 = [L_b; E_b; E_Hb; 0; 100; L_b; L_b; 0; L_b]; % pack initial conditions
par_LEHRS = [TC, f, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL = (interp1(t(:), LEHRS(:,1), tL(:,1)))/del_M;

par_LEHRS = [TC, f1, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL1 = (interp1(t(:), LEHRS(:,1), tL1(:,1)))/del_M;
par_LEHRS = [TC, f2, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL2 = (interp1(t(:), LEHRS(:,1), tL2(:,1)))/del_M;
par_LEHRS = [TC, f3, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL3 = (interp1(t(:), LEHRS(:,1), tL3(:,1)))/del_M;
par_LEHRS = [TC, f4, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL4 = (interp1(t(:), LEHRS(:,1), tL4(:,1)))/del_M;

par_LEHRS = [TC15, f, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL15 = (interp1(t(:), LEHRS(:,1), tL15(:,1)))/del_M;
par_LEHRS = [TC196, f, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL196 = (interp1(t(:), LEHRS(:,1), tL196(:,1)))/del_M;
par_LEHRS = [TC21, f, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL21 = (interp1(t(:), LEHRS(:,1), tL21(:,1)))/del_M;
par_LEHRS = [TC244, f, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL244 = (interp1(t(:), LEHRS(:,1), tL244(:,1)))/del_M;
par_LEHRS = [TC267, f, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL267 = (interp1(t(:), LEHRS(:,1), tL267(:,1)))/del_M;

par_LEHRS = [TC15, f2, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL15f = (interp1(t(:), LEHRS(:,1), tL15f(:,1)))/del_M;
par_LEHRS = [TC196, f2, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL196f = (interp1(t(:), LEHRS(:,1), tL196f(:,1)))/del_M;
par_LEHRS = [TC21, f2, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL21f = (interp1(t(:), LEHRS(:,1), tL21f(:,1)))/del_M;
par_LEHRS = [TC267, f2, par_LEHR, par_CS]; [t, LEHRS, ~, ~] = odehax(LEHRS_0, par_LEHRS, 15, E_0); % ODE simulation
EL267f = (interp1(t(:), LEHRS(:,1), tL267f(:,1)))/del_M;

% t-S data for larvae
LEHRS_0 = [L_p; E_m * L_p^3; E_Hp; 0; 100; L_b; L_p; 0; L_p]; % pack initial conditions
% starvation data; starvation experiment initiated after 2-3 days after birth:
par_LEHRS = [TC_20, 0, par_LEHR, par_CS]; % combine TC, f, and parameters
[t, LEHRS] = ode45(@get_LEHRS, tS(:,1), LEHRS_0 , options, par_LEHRS);
ES = LEHRS(:,5)/100;
% background mortality for starvation data
par_LEHRS = [TC_20, 1, par_LEHR, par_CS]; % combine TC, f, and parameters
[t, LEHRS] = ode45(@get_LEHRS, tS(:,1), LEHRS_0 ,options, par_LEHRS);
ES2 = LEHRS(:,5)/100;

% weight - respiration
options = odeset;
E_b = p_Am/ v * L_b^3; % reserve at birth set at maximum reserve density (parent generation fed ad libitum)
par_LEHR = [p_Am, v, kap, kap_R, p_M, k_J, E_G, E_Hb, E_Hp, E_Rj, E_He, kap_V]; % pack parameters
par_CS = [s_shrink, kap_G, h_b_shrink]; % pack survival parameters
LEHRS_0 = [L_b; E_b; E_Hb; 0; 100; L_b; L_b; 0; L_b]; % pack initial conditions

par_LEHRS = [TC_20, f5, par_LEHR, par_CS]; % combine TC, f, and parameters
[t, LEHRS] = ode45(@get_LEHRS, [0 150], LEHRS_0, options, par_LEHRS); % ODE simulation w.o. events
WL3 = LEHRS(:,1).^3 .* d_V; % g, dry weight of structure
WE = LEHRS(:,2).* w_E./ mu_E; % g, dry weight of reserve
WR = LEHRS(:,4).* w_E./ mu_E; % g, dry weight of reproduction buffer
W = (WL3 + WE + WR) .* 1e3;  % mg, total dry weight
L = (interp1(W(:), LEHRS(:,1), WdJO_20(:,1))); % cm, length at Edwa1956 data points
E_H = (interp1(W(:), LEHRS(:,3), WdJO_20(:,1))); % cm, maturity at Edwa1956 data points
Wd =  WdJO_20(:,1); % g, total dry weight at Edwa1956 data points
M = min(L, L_p)/L_b;
pA = f5 * p_Am * M .* L.^2; % J/d, assimilation flux
E = f5 * (p_Am / v) .* L.^3;
pS = p_M * L.^3; % J/d, somatic maintenance flux
r = (E ./ L.^3 .* v .* M ./ L - (pS ./ L.^3 ./ kap)) ./ (E ./ L.^3 + E_G ./ kap);
pC = E .* (v .* M ./ L - r); % J/d, mobilization flux
pG = kap * pC - pS;   % J/d, growth flux
pJ = k_J * E_H;  % J/d, maturity maintenance flux
pR = (1 - kap) .* pC - pJ; % J/d, maturation/reproduction
pD = pS + pJ + (1 - kap_R) * pR; % J/d, dissipation flux
JM = -(n_M\ n_O) * eta_O * [pA, pD, pG]'; % mol/d, mineral fluxes (J_C, J_H, J_O, J_N in rows)
EJT_O20 = -1e3 *  JM(3,:)' * TC_20/ 24/ 32e-6 ./Wd; % µg O2/h.mg, O2 consumption

par_LEHRS = [TC_10, f5, par_LEHR, par_CS]; % combine TC, f, and parameters
[t, LEHRS] = ode45(@get_LEHRS, [0 150], LEHRS_0, options, par_LEHRS); % ODE simulation w.o. events
WL3 = LEHRS(:,1).^3 .* d_V; % g, dry weight of structure
WE = LEHRS(:,2).* w_E./ mu_E; % g, dry weight of reserve
WR = LEHRS(:,4).* w_E./ mu_E; % g, dry weight of reproduction buffer
W = (WL3 + WE + WR) .* 1e3;  % mg, total dry weight
L = (interp1(W(:), LEHRS(:,1), WdJO_10(:,1))); % cm, length at Edwa1956 data points
E_H = (interp1(W(:), LEHRS(:,3), WdJO_10(:,1))); % cm, maturity at Edwa1956 data points
Wd =  WdJO_10(:,1); % g, total dry weight at Edwa1956 data points
M = min(L, L_p)/L_b;
pA = f5 * p_Am * M .* L.^2; % J/d, assimilation flux
E = f5 * (p_Am / v) * L.^3;
pS = p_M * L.^3; % J/d, somatic maintenance flux
r = (E ./ L.^3 .* v .* M ./ L - (pS ./ L.^3 ./ kap)) ./ (E ./ L.^3 + E_G ./ kap);
pC = E .* (v .* M ./ L - r); % J/d, mobilization flux
pG = kap * pC - pS;   % J/d, growth flux
pJ = k_J * E_H;  % J/d, maturity maintenance flux
pR = (1 - kap) .* pC - pJ; % J/d, maturation/reproduction
pD = pS + pJ + (1 - kap_R) * pR; % J/d, dissipation flux
JM = -(n_M\ n_O) * eta_O * [pA, pD, pG]'; % mol/d, mineral fluxes (J_C, J_H, J_O, J_N in rows)
EJT_O10 = -1e3 *  JM(3,:)' * TC_10/ 24/ 32e-6 ./Wd; % µg O2/h.mg, O2 consumption

% pack to output
prdData.tL   = EL;
prdData.tL_m = EL_m;
prdData.LW   = EWd;
prdData.tL1  = EL1;
prdData.tL2  = EL2;
prdData.tL3  = EL3;
prdData.tL4  = EL4;

prdData.tL15  = EL15;
prdData.tL196 = EL196;
prdData.tL21  = EL21;
prdData.tL244 = EL244;
prdData.tL267 = EL267;

prdData.tL15f  = EL15f;
prdData.tL196f = EL196f;
prdData.tL21f  = EL21f;
prdData.tL267f = EL267f;

prdData.tS  = ES;
prdData.tS2 = ES2;

prdData.WdJO_20 = EJT_O20;
prdData.WdJO_10 = EJT_O10;

%--------------------------------------------------------------------------
%------------------ Start gaiac 2022 data predictions ---------------------
%--------------------------------------------------------------------------

par_LEHR = [p_Am, v, kap, kap_R, p_M, k_J, E_G, E_Hb, E_Hp, E_Rj, E_He, kap_V]; % pack parameters
pars_UE0 = [V_Hb, g, k_J, k_M, v]; % compose parameter vector for E_0
par_CS = [s_shrink, kap_G, h_b_shrink]; % pack survival parameters
LEHRS_0 = [L_b; E_b; E_Hb; 0; 100; L_b; L_b; 0; L_b]; % pack initial conditions

  TC12 = tempcorr(C2K(11.8), T_ref, pars_T);
  par_LEHRS = [TC12, f, par_LEHR, par_CS]; % combine TC, f, and parameters
%   [t, LEHRS, te, LEHRSe] = odehax(LEHRS_0, par_LEHRS, 40, E_0); % ODE simulation with life events
  [t, LEHRS] = ode45(@get_LEHRS, [0 80], LEHRS_0, options, par_LEHRS); % ODE simulation w/o life events
  WL3 = (interp1(t(:), LEHRS(:,1), tW12(:,1))).^3 .* d_V; % g, dry weight of structure
  WE = (interp1(t(:), LEHRS(:,2), tW12(:,1))).* w_E./ mu_E; % g, dry weight of reserve
  WR = (interp1(t(:), LEHRS(:,4), tW12(:,1))).* w_E./ mu_E; % g, dry weight of reproduction buffer
  W12 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
  prdData.tW12 = W12;

  TC15 = tempcorr(C2K(15.2), T_ref, pars_T);
  par_LEHRS = [TC15, f, par_LEHR, par_CS]; % combine TC, f, and parameters
%   [t, LEHRS, te, LEHRSe] = odehax(LEHRS_0, par_LEHRS, 40, E_0); % ODE simulation with life events
  [t, LEHRS] = ode45(@get_LEHRS, [0 80], LEHRS_0, options, par_LEHRS); % ODE simulation w/o life events
  WL3 = (interp1(t(:), LEHRS(:,1), tW15(:,1))).^3 .* d_V; % g, dry weight of structure
  WE = (interp1(t(:), LEHRS(:,2), tW15(:,1))).* w_E./ mu_E; % g, dry weight of reserve
  WR = (interp1(t(:), LEHRS(:,4), tW15(:,1))).* w_E./ mu_E; % g, dry weight of reproduction buffer
  W15 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
  prdData.tW15 = W15;

  TC20 = tempcorr(C2K(20.2), T_ref, pars_T);
  par_LEHRS = [TC20, f, par_LEHR, par_CS]; % combine TC, f, and parameters
%   [t, LEHRS, te, LEHRSe] = odehax(LEHRS_0, par_LEHRS, 40, E_0); % ODE simulation with life events
  [t, LEHRS] = ode45(@get_LEHRS, [0 80], LEHRS_0, options, par_LEHRS); % ODE simulation w/o life events
  WL3 = (interp1(t(:), LEHRS(:,1), tW20(:,1))).^3 .* d_V; % g, dry weight of structure
  WE = (interp1(t(:), LEHRS(:,2), tW20(:,1))).* w_E./ mu_E; % g, dry weight of reserve
  WR = (interp1(t(:), LEHRS(:,4), tW20(:,1))).* w_E./ mu_E; % g, dry weight of reproduction buffer
  W20 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
  prdData.tW20 = W20;

  TC23 = tempcorr(C2K(23.2), T_ref, pars_T);
  par_LEHRS = [TC23, f, par_LEHR, par_CS]; % combine TC, f, and parameters
%   [t, LEHRS, te, LEHRSe] = odehax(LEHRS_0, par_LEHRS, 40, E_0); % ODE simulation with life events
  [t, LEHRS] = ode45(@get_LEHRS, [0 80], LEHRS_0, options, par_LEHRS); % ODE simulation w/o life events
  WL3 = (interp1(t(:), LEHRS(:,1), tW23(:,1))).^3 .* d_V; % g, dry weight of structure
  WE = (interp1(t(:), LEHRS(:,2), tW23(:,1))).* w_E./ mu_E; % g, dry weight of reserve
  WR = (interp1(t(:), LEHRS(:,4), tW23(:,1))).* w_E./ mu_E; % g, dry weight of reproduction buffer
  W23 = (WL3 + WE + WR) .* 1e3 ; % mg, dry weight
  prdData.tW23 = W23;

%--------------------------------------------------------------------------
%------------------- End gaiac 2022 data predictions ----------------------
%--------------------------------------------------------------------------

end

%%
% Overarching function to merge consecutive ode45-runs of the different hax life stages
function [t, LEHRS, te, LEHRSe] = odehax(LEHRS_0, par_LEHRS, tfinal, E_0)
tstart = 0;
refine = 4;

E_G = par_LEHRS(9); % J/cm^3, Spec cost for structure
kap_V = par_LEHRS(14); % -, Conversion efficiency from larval reserve to larval structure, back to imago reserve
kap_R = par_LEHRS(6); % -, Reproduction efficiency

% simulate until pupa stage
options = odeset('Events', @pupationEvent, 'Refine', refine);
[t, LEHRS, te, LEHRSe] = ode45(@get_LEHRS, [tstart tfinal], LEHRS_0, options, par_LEHRS); % ODE simulation
% simulate pupa until emergence
if (t(end) < tfinal)
    LEHRSet = LEHRSe; % initial states for pupa stage
    LEHRSet(2) = LEHRSe(2) + LEHRSe(1)^3 * E_G * kap_V; % Reserve expanded by transformed structure
    LEHRSet(1) = 0.000001; LEHRSet(6) = LEHRSet(1); LEHRSet(7) = LEHRSet(1); LEHRSet(9) = LEHRSet(1); % Structure to zero
    LEHRSet(3) = 0; % Maturity to zero
    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation. 'refine' is 4 by default.
    nt = length(t);
    options = odeset('Events', @emergenceEvent, 'Refine', refine, 'InitialStep', t(nt)-t(nt-refine), 'MaxStep',t(nt)-t(1));
    [t2, LEHRS2, te2, LEHRSe2] = ode45(@get_LEHRS, [te tfinal], LEHRSet, options, par_LEHRS); % ODE simulation
    % Accumulate output
    nt = length(t2);
    t = [t; t2(2:nt)];
    LEHRS = [LEHRS; LEHRS2(2:nt,:)];
    te = [te; te2]; % time at event
    LEHRSe = [LEHRSe; LEHRSe2]; % states at event
end
% simulate imago after emergence
if (t(end) < tfinal)
    LEHRSe2t = LEHRSe2; % initial states for pupa stage
    LEHRSe2t(8) = LEHRSe2(4) * kap_R / E_0; % Make eggs from repro buffer
    LEHRSe2t(4) = 0; % Repro buffer to zero
    % A good guess of a valid first timestep is the length of the last valid
    % timestep, so use it for faster computation. 'refine' is 4 by default.
    nt = length(t2);
    options = odeset('Refine', refine, 'InitialStep', t2(nt)-t2(nt-refine),'MaxStep',t2(nt)-t2(1));
    [t3, LEHRS3] = ode45(@get_LEHRS, [te2 tfinal], LEHRSe2t, options, par_LEHRS); % ODE simulation
    % Accumulate output
    nt = length(t3);
    t = [t; t3(2:nt)];
    LEHRS = [LEHRS; LEHRS3(2:nt,:)];
end
end

%%
% Event function to determine the moment of pupation
function [value, isterminal, direction] = pupationEvent(t, LEHRS, par_LEHRS)
L       = LEHRS(1); % cm, Structural length
E_R     = LEHRS(4); % # J, Cumulative energy in reproduction buffer
E_Rj    = par_LEHRS(12); % J/cm^3, Reproduction buffer density at pupation

value = E_R - (E_Rj * L^3); % The value that we want to be zero
isterminal = 1;  % Halt integration
direction = 0;   % The zero can be approached from either direction
end

%%
% Event function to determine the moment of emergence
function [value, isterminal, direction] = emergenceEvent(t, LEHRS, par_LEHRS)
E_H     = LEHRS(3); % J, Energy invested in maturity
E_R     = LEHRS(4); % J, Cumulative energy in reproduction buffer
E_He    = par_LEHRS(13); % J, Maturity at emergence

value = (E_H - E_He) + (E_R == 0) ; % The value that we want to be zero
isterminal = 1;  % Halt integration
direction = 0;   % The zero can be approached from either direction
end

%%
% ODE system of the hax model; different life-stages are covered by if-clauses
% The system is designed to work with the event functions for pupa and emergence
% If executed without these, it should behave like an abj model w/o egg laying
function [dLEHRS] = get_LEHRS(t, LEHRS, par_LEHRS)

% Read parameters
TC    = par_LEHRS(1); % -, Temperature correction
f     = par_LEHRS(2); % -, Functional response
p_Am  = par_LEHRS(3); % J / d.cm^2, Surface-area specific maximum assimilation rate
v     = par_LEHRS(4); % cm/d, Energy conductance
kap   = par_LEHRS(5); % -, Allocation fraction to soma
% kap_R = par_LEHRS(6); % -, Reproduction efficiency
p_M   = par_LEHRS(7); % J/d.cm^3, Vol-spec somatic maint
k_J   = par_LEHRS(8); % 1/d, Maturity maint rate coefficient
E_G   = par_LEHRS(9); % J/cm^3, Spec cost for structure
E_Hb  = par_LEHRS(10); % J, Maturity at birth
E_Hp  = par_LEHRS(11); % J, Maturity at puberty
E_Rj  = par_LEHRS(12); % J/cm^3, Reproduction buffer density at pupation
E_He  = par_LEHRS(13); % J, Maturity at emergence
% kap_V  = par_LEHRS(14); % -, Conversion efficiency from larval reserve to larval structure, back to imago reserve
s_shrink = par_LEHRS(15); % -, Shrinking stress coefficient
kap_G = par_LEHRS(16); % -, Growth efficiency
h_b_shrink = par_LEHRS(17); % -, Background hazard rate coefficient for starvation experiment;

% initialize state variables
L       = LEHRS(1); % cm, Structural length
E       = LEHRS(2); % J, Energy in reserve
E_H     = LEHRS(3); % # J, Energy invested in maturity
E_R     = LEHRS(4); % # J, Cumulative energy in reproduction buffer
S       = LEHRS(5); % -, Surviving fraction
Lsa     = LEHRS(6); % cm, Min(Current Length, Length at start of acceleration)
Lea     = LEHRS(7); % cm, Min(Current Length, Length at end of acceleration)
eggs    = LEHRS(8); % #, Produced eggs
maxL    = LEHRS(9); % cm, Maximum length reached so far

% Metabolic acceleration (Type M)
M = Lea / Lsa; % Acceleration factor
p_Am_M = p_Am * M; % Metabolic acceleration affects p_Am
v_M = v * M; % Metabolic acceleration affects v

% Scaled reserve density (to check for starvation)
e = E * v / (L^3 * p_Am);
% Scaled structural length (to check for starvation)
l = L * p_M / (kap * p_Am);

% Assimilation flux (from birth until pupation)
if (E_H < E_Hb) && (E_R == 0) ||... % Before birth
        (E_R > (E_Rj * L^3)) && (E_H < E_He) ||... % During pupa stage
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
else
    r = (E / L^3 * v_M / L - (p_S / L^3 / kap)) / (E / L^3 + E_G * kap_G^2 / kap);
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

% Changes in state variables
dL = L * r / 3; % DEB Book 3rd ed. Chap. 2.6.1
dE = p_A - p_C; % DEB Book 3rd ed. Chap. 2.3

if (E_H < E_Hp) || (E_R > (E_Rj * L^3)) && (E_H < E_He)
    dE_H = p_R;
    dE_R = 0;
else
    dE_H = 0;
    dE_R = p_R;
end

if (E_H < E_Hb) && (E_R < (E_Rj * L^3))
    dLsa = dL;
else
    dLsa = 0;
end

if (E_H < E_Hp) && (E_R < (E_Rj * L^3))
    dLea = dL;
else
    dLea = 0;
end

deggs = 0;

% Shrinking stress
dmaxL = max(0, dL); % cm, max stuctural length
k_M = p_M/ E_G;  % 1/d, maturity maintenance rate coefficient
h_sh = s_shrink * k_M * max(0, maxL - L)/maxL; % hazard due to shrinking

dS = - S * (h_sh + h_b_shrink); % 1/d, survival fraction

dLEHRS = [dL; dE; dE_H; dE_R; dS/ TC; dLsa; dLea; deggs; dmaxL] * TC; % collect derivatives
end
