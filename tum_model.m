function [f, g, p, t] = tum_model(delta_vec, tau_vec, params)

    dt = params.dt;
    T  = params.T;
    lambdaA   = params.lambdaA; 
    lambdaK   = params.lambdaK;
    lambdaD   = params.lambdaD; 
    lambdaP   = params.lambdaP;
    g0        = params.g0; 
    f0        = params.f0;
    sigma     = params.sigma;
    % delta_vec = params.delta_vec;
    % kappa     = params.kappa  ; 
    % tau_vec   = params.tau_vec; 
    
    
    %since floating point error
    t = T/dt;
    p = zeros(t+1, 1);
    time = 0:dt:T;
    
    for ind = 1:t+1
        for i = 1:length(delta_vec)
            p(ind) = p(ind) + delta_vec(i)*(1/(sigma*(2*pi)^(1/2))) * exp((-1 * (abs(tau_vec(i) - time(ind))^2)/(2 * sigma^2)));
        end
    end
    
    
    g = zeros(t,1);
    f = zeros(t,1);
    
    g(1) = g0;
    f(1) = f0;
    
    for i = 1:t
        g(i+1) = g(i) + dt * (lambdaP * g(i) * (1 - g(i)) - lambdaA * g(i) - lambdaK * g(i) * f(i));
        f(i+1) = f(i) + dt * (-1 * lambdaD * f(i) + p(i));
    end                                                                               

    t = time;
end