%%-----------------------------------------------------------------------------------------------------------------------------%%
%%--         Automated Inverse Finite Element Method Matlab Script Using Bayesian Optimization                               --%%
%%--                        Date : 23 November 2021                                                                               --%%
%%--                       Authors : Sourabh Shende, Teng Long                                                               --%%
%% 

%% This MATLAB script accompanies the following manuscript:
%% Teng Long, Sourabh Shende, Chia-Ying Lin, Kumar Vemaganti, "Experiments and Hyperelastic Modeling of Porcine Meniscus Show Heterogeneity at High Strains" submitted to Biomechanics and Mechanobiology (2022).
%%-----------------------------------------------------------------------------------------------------------------------------%%

clc
clear
close all
format long

training_set = 1;
%%-----------------------------------------------------------------------%%
fclose('all');

currentfolderpath = pwd; % path to current folder
dynapath = 'C:\LSDYNA\program\ls-dyna_smp_d_R10.0_winx64_ifort160.exe'; % path to LS-Dyna .exe file
K_Address = 'Test31824.k'; % base FE Input filename


ini_training_point = 5; % # initial training points
iter_max = 50;          % stopping criteria == # of FE evaluations   
command = strcat('"',dynapath,'"',' i=',currentfolderpath,'\',K_Address,' ncpu=12 memory=20m'); % Command to launch LS-Dyna 
 
nd = 2;   % number of dimensions
x_scal = 3.82e6;        % Scaling factor to reduce domain to [0,1]
lb_acf = [2.1e5,0]./x_scal; % lower bound of the parameters
ub_acf = [3.82e6,8.8e5]./x_scal; % upper bound of the parameters 

%----------------------------------------------------------------------%
%---------      Loading Experimental Data          --------------------%
%----------------------------------------------------------------------%
load('Anterior_Mean_101_20P.mat')
load('Anterior_SD_101_20P.mat')

Exp_Mean(1,1) = 0; 
Exp_SD(1,1) = 100;

% This is matelab processed data. The row experiment data has been uploaded 
% to https://doi.org/10.7945/1k0e-ew43
%----------------------------------------------------------------------%

%----------------------------------------------------------------------%
%---------      Initial guesses for  BO            --------------------%
%----------------------------------------------------------------------%
rng default    % Random number generation
tic
t_ini_train_start = toc;
load('Material_Matrix_Total.mat') % File with initial 30 FE simulations 
r = randperm(30);
training_set_count = 1;
r = r((training_set_count-1)*ini_training_point+1:training_set_count*ini_training_point);    

x_train_i = Material_Matrix_Total(r,1:2)./x_scal;   
y_train_i = zeros(ini_training_point,1);      
for ii= 1:ini_training_point
    FE_response = Material_Matrix_Total(r(ii),5:105);        
    y_train_i(ii,1) = mean_square_error(Exp_Mean,Exp_SD,FE_response);
end        
t_ini_train_end = toc;
time_ini_train = t_ini_train_end - t_ini_train_start;   
%----------------------------------------------------------------------%


%----------------------------------------------------------------------%
%---------------      BO hyperparameters           --------------------%
%----------------------------------------------------------------------%
kernel_length_scale = 0.5*sqrt(nd);
kernel_scaling_parameter = 1*sqrt(1e3); 
acq_para = 3; % kappa for lower confidence bound (LCB) acquisition function
output_noise = 0.00001;
%----------------------------------------------------------------------%

%----------------------------------------------------------------------%
%------------    Initializing output variables      -------------------%
%----------------------------------------------------------------------%
time_acq_func = zeros(iter_max,1);
time_hyper_opt = zeros(iter_max,1);
time_obj_func = zeros(iter_max,1);
acq_func_eval = zeros(iter_max,1);
acq_func_val = zeros(iter_max,1);
mu_val = zeros(iter_max,1);
condition_num_Kf = zeros(iter_max,1);
determinant_Kf = zeros(iter_max,1);
parameters = zeros(iter_max,4);
%----------------------------------------------------------------------%


%----------------------------------------------------------------------%
%--------------    Generating output files        ---------------------%
%----------------------------------------------------------------------%
tic
x_train = x_train_i;
y_train = y_train_i;
file_name_1 = strcat('obj_funct_ts_',num2str(training_set_count),'_2D_param_est');
file_name_2 = strcat('design_var_ts_',num2str(training_set_count),'_2D_param_est');
file_name_3 = strcat('Output_info_ts_',num2str(training_set_count),'_2D_param_est');
file_name_4 = strcat('parameters_ts_',num2str(training_set_count),'_2D_param_est');
file_name_5 = strcat('obj_funct_ts_',num2str(training_set_count),'_2D_param_est.txt');
fileID_5 = fopen(file_name_5,'a');
fprintf(fileID_5,'%12.8f\n',y_train);
file_name_6 = strcat('design_var_ts_',num2str(training_set_count),'_2D_param_est.txt');
fileID_6 = fopen(file_name_6,'a');
fprintf(fileID_6,'%d %d\n',x_train.');
file_name_7 = strcat('Output_info_ts_',num2str(training_set_count),'_2D_param_est.txt');
fileID_7 = fopen(file_name_7,'a'); 
fprintf(fileID_7,'%12s %12s %12s %12s %12s %12s %12s %12s\n','condition_num_Kf','determinant_Kf','time_hyper_opt','acq_func_val'...
    ,'mu_val','acq_func_eval','time_acq_func','time_fem');
file_name_8 = strcat('parameters_ts_',num2str(training_set_count),'_2D_param_est.txt');
fileID_8 = fopen(file_name_8,'a');   
fprintf(fileID_8,'%12s %12s %12s %12s\n','kernel_length_scale','kernel_scaling_parameter','output_noise','acq_para');
%----------------------------------------------------------------------%


%%-----------------------------------------------------------------------%%
%%----------------- Main optimization iterations    ---------------------%%
%%-----------------------------------------------------------------------%%

for iter=1:iter_max
    disp(['Iteration = ' num2str(iter)])    
    t_hyper_opt_start = toc;        
    t_acq_start = toc;
    %----------------------------------------------------------------------%
    %--------------    Evaluating next sampling point ---------------------%
    %----------------------------------------------------------------------%
    [K] = squared_exponential_n_dim(x_train,x_train,kernel_length_scale,kernel_scaling_parameter);  % Constructing covariance function
    Kf = K + output_noise^2*eye(size(K,2));      
    dK = decomposition(Kf);
    Kf_inv_times_ytrain = dK\y_train; 
    [x_next_best,max_acq_func,mu_next_iter,count] = next_sample(x_train,kernel_length_scale,kernel_scaling_parameter,output_noise,nd,acq_para,lb_acf,ub_acf,K,dK,Kf_inv_times_ytrain);
    x_train(end+1,:) = x_next_best;
    %----------------------------------------------------------------------%
    parameters(iter,:) = [kernel_length_scale ...
    kernel_scaling_parameter output_noise acq_para];
    cond_Kf = cond(Kf);
    det_Kf = det(Kf);
    t_acq_end = toc;                                                    
    time_acq_func(iter,1) = t_acq_end - t_acq_start;
    acq_func_val(iter,1) = max_acq_func;
    mu_val(iter,1) = mu_next_iter;
    condition_num_Kf(iter,1) = cond_Kf;
    determinant_Kf(iter,1) = det_Kf;
    acq_func_eval(iter,1) = count;
    size_end = size(x_train,1);                              
    %----------------------------------------------------------------------%
    %-------    Evaluating expensive objective function  ------------------%
    %----------------------------------------------------------------------%
    t_obj_func_start = toc;    
    [obj_fval] = get_residual(x_next_best,Exp_Mean,Exp_SD,x_scal,command,K_Address);    
    y_train(end+1,1) = obj_fval;
    t_obj_func_end = toc;
    disp(['y_train = ',num2str(y_train(end))])
    time_obj_func(iter,1) = t_obj_func_end - t_obj_func_start;
    %----------------------------------------------------------------------%
    
    fprintf(fileID_5,'%12.8f\n',y_train(end));                          
    fprintf(fileID_6,'%d %d\n',x_train(end,:));        
    fprintf(fileID_7,...
        '%12.8f %5e %12.8f %12.8f %12.8f %d %12.8f %12.8f %12.8f\n',...
        condition_num_Kf(iter,1),determinant_Kf(iter,1),time_hyper_opt(iter,1),...
        acq_func_val(iter,1),mu_val(iter,1),acq_func_eval(iter,1),...
        time_acq_func(iter,1),time_obj_func(iter,1),...
        time_ini_train);
    fprintf(fileID_8,'%12.8f %12.8f %12.8f %12.8f\n',kernel_length_scale,...
        kernel_scaling_parameter,output_noise,acq_para);

end

fclose(fileID_5);
fclose(fileID_6);
fclose(fileID_7);
fclose(fileID_8);

save((file_name_1),'y_train');
save((file_name_2),'x_train');
output_info = [condition_num_Kf,determinant_Kf,time_hyper_opt,acq_func_val,mu_val,acq_func_eval,...
    num_of_generation,time_acq_func,time_obj_func,norm_D_loglikelihood_vec,D_optimality,time_ini_train*ones(size(time_obj_func,1))];
save((file_name_3),'output_info');                 
save((file_name_4),'parameters');
toc
t = toc;
disp(['Time elapsed = ',num2str(t)])
[y_max,I] = max(y_train);
x_max = x_train(I,:);
[val,idx] = min(y_train);
disp(['x_train = ',num2str(x_train(idx,:))])
disp(['obj funct value = ',num2str(val)])
%%-----------------------------------------------------------------------------------------------------------------------------%%
%----------------------                                End of code                               -------------------------------%
%%-----------------------------------------------------------------------------------------------------------------------------%%

