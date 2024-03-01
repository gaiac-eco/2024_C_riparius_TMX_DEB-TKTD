function [par, metaPar, txtPar] = pars_init_Chironomus_riparius(metaData)

metaPar.model = 'hax'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 7377;       free.T_A   = 1;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.z = 0.09386;      free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.01537;      free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.3995;     free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.kap_V = 0.7;      free.kap_V = 0;   units.kap_V = '-';        label.kap_V = 'conversion efficiency E -> V -> E'; 
par.p_M = 3474;       free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 4400;       free.E_G   = 0;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 0.005305;  free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hp = 0.01744;   free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.E_He = 0.003788;  free.E_He  = 1;   units.E_He = 'J';         label.E_He = 'maturity at emergence'; 
par.E_Rj = 3.541e+04; free.E_Rj  = 1;   units.E_Rj = 'J/cm^3';    label.E_Rj = 'reproduction buffer density at pupation'; 
par.h_a = 0.0001;     free.h_a   = 0;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.T_AH = 5949;      free.T_AH  = 1;   units.T_AH = 'K';         label.T_AH = 'Arrhenius temperature for upper boundary';
par.T_H = 293.6;      free.T_H   = 1;   units.T_H = 'K';          label.T_H = 'upper boundary'; 
par.T_AL = 0;         free.T_AL  = 0;   units.T_AL = 'K';         label.T_AL = 'Arrhenius temperature for lower boundary'; 
par.T_L = 10;         free.T_L   = 0;   units.T_L = 'K';          label.T_L = 'lower boundary';
par.del_M = 0.04925;  free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coefficient for larva length'; 
par.f = 1;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var and ad libitum data'; 
par.f1 = 0.7105;      free.f1    = 1;   units.f1 = '-';           label.f1 = 'sc. funct. res.  0.1 mg tetramin per ind'; 
par.f2 = 0.7127;      free.f2    = 1;   units.f2 = '-';           label.f2 = 'sc. funct. res.  0.2 mg tetramin per ind'; 
par.f3 = 0.7517;      free.f3    = 1;   units.f3 = '-';           label.f3 = 'sc. funct. res.  0.3 mg tetramin per ind'; 
par.f4 = 0.8958;      free.f4    = 1;   units.f4 = '-';           label.f4 = 'sc. funct. res.  0.4 mg tetramin per ind'; 
par.f5 = 0.5799;      free.f5    = 1;   units.f5 = '-';           label.f5 = 'scaled functional response for respiration data'; 
par.z_m = 0.09377;    free.z_m   = 1;   units.z_m = '-';          label.z_m = 'zoom factor for males';
par.h_b_shrink = 0.02295;  free.h_b_shrink = 1;  units.h_b_shrink = '-';  label.h_b_shrink = 'background hazard rate coefficient for starv data';
par.s_shrink = 4.131;      free.s_shrink = 1;    units.s_shrink = '-';    label.s_shrink = 'shrinking stress coefficient'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
