function [acq_f_max_Xtest_val,Fval,mu,count] = next_sample(x_train,kernel_length_scale,kernel_scaling_parameter,output_noise,nd,acq_para,lb_acf,ub_acf,K,dK,Kf_inv_times_ytrain)         
    function f = myfunc(x_test)
        [K_ss] = squared_exponential_n_dim(x_test,x_test,kernel_length_scale,kernel_scaling_parameter);
        [K_s] = squared_exponential_n_dim(x_train,x_test,kernel_length_scale,kernel_scaling_parameter);
        [mu,stdv] = surrogate_model(K_s,K_ss,dK,Kf_inv_times_ytrain);
        kappa = acq_para;  
        f = 1*lower_confidence_bound(mu,stdv,kappa);
    end    
    %----------------------------------------------------------------------%
    %-------    Genetic Algorithm for minimizing LCB  ---------------------%
    %----------------------------------------------------------------------%
    nonlcon = [];    
    A = [];
    b = [];
    Aeq = [];
    beq = [];    
    rng default
    [acq_f_max_Xtest_val,Fval,exitFlag,Output,population,Best]= ga(@(x_test) myfunc(x_test),nd,A,b,Aeq,beq,lb_acf,ub_acf,nonlcon);       
    count = Output.funccount;   
end
