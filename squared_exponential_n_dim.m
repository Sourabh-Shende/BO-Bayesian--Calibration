function K = squared_exponential_n_dim(xi,xj,kernel_length_scale,kernel_sigma)
    K = zeros(size(xi,1),size(xj,1));
    for i=1:size(xi,1)
        for j=1:size(xj,1) 
            r = sqrt(sum(xi(i,1:end).^2 - 2*xi(i,1:end).*xj(j,1:end) + xj(j,1:end).^2));
            K(i,j) =kernel_sigma.^2*exp(-0.5*((1/kernel_length_scale)*r)^2);                       
        end
    end
end