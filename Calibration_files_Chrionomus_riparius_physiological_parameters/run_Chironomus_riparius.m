close all; 
global pets

pets = {'Chironomus_riparius'};
check_my_pet(pets);

estim_options('default'); 
estim_options('max_step_number', 500); 
estim_options('max_fun_evals', 10000);

estim_options('pars_init_method', 2);
% get initial parameters from:
% 0 = automatized computation
% 1 = results_my_pet.mat file
% 2 = pars_init_my_pet.m file
estim_options('results_output', -2);

% 1 = show output on screen and write file
estim_options('method', 'no');
% 'no' = do not estimate
% 'nm' = use Nelder-Mead method
estim_pars;

% % write results_my_pet.mat to pars_init_my_pet.m
% mat2pars_init('Chironomus_riparius')