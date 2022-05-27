function [obj_fval] = get_residual(x_train,Exp_Mean,Exp_SD,x_scal,command,K_Address)
%----------------------------------------------------------------------%
%---------   Calculate the residual value             -----------------%
%----------------------------------------------------------------------%    
    
    cm1 = x_train(1,1)*x_scal; % scale back to real value for C10
    cm4 = x_train(1,2)*x_scal; % scale back to real value for C20
    cm6 = 0;
	
	[Line_ID] = GetDynaReplace(K_Address, cm1,cm4,cm6);  % update C10 and C20 values for FEA input file
    
    [status,cmdout] = system(command)
    FEA_RF_ADD = 'rwforc';  % 'rwforc' is and reaction force file from LS-Dyna in current directory
    [FEA_RF_T] = get_reactfrcFEA(FEA_RF_ADD); % extract the reaction force from FEA model
    
    [mse] = mean_square_error(Exp_Mean,Exp_SD,FEA_RF_T);% calculate the mean square error based on the experiment data and FEA result
    obj_fval = mse;
%----------------------------------------------------------------------%    
end


