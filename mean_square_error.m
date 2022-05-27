function [mse] = mean_square_error(EXP_RF_mean_interp1,EXP_RF_std_interp1,FE_response)
%-----------------------------------------------------------------------------------------------%
%-------   Calibrate the model on the mean and SD of experimeent data and FEA result   ---------%
%-----------------------------------------------------------------------------------------------%  
% EXP_RF_mean_interp1: mean of experimeent data
% FE_response: FEA result
% EXP_RF_std_interp1: SD of experimeent data
   mse = sum((EXP_RF_mean_interp1 - FE_response).^2./(2*EXP_RF_std_interp1.^2));
%-----------------------------------------------------------------------------------------------%  
end