function [mu,stdv]= surrogate_model(K_s,K_ss,dK,Kf_inv_times_ytrain) 
    mu = K_s.'*Kf_inv_times_ytrain;
    var = K_ss - K_s.'*(dK\K_s);
    stdv = sqrt(var);
end
