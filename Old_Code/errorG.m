%% Problem 1: Error Computation.

function ek = errorG()
    
    global T g0 f0;
    
    old = 0;
    new = 0;
    e = [];
    
    for k = 0:9
        % Get the last value of the discretization.
        dt = 0.01 * T / 2^k;
        [~, g] = discretization(dt, g0, f0);
        old = g(end);
        
        % Get the last value of the discretization.
        dt = 0.01 * T / 2^(k+1);
        [~, g] = discretization(dt, g0, f0);
        new = g(end);
        
        e(k+1) = 100 * ((new - old) / old);
    end
    
    ek = e;
end
