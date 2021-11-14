%% Discretization Problem 1

function [fSolve, gSolve] = discretization(h, g0, f0, p)

    global lambdaA lambdaK lambdaD;
    
    n = 50 / h;
    
    g = zeros(n,1);
    f = zeros(n,1);
    
    g(1) = g0;
    f(1) = f0;
    
    for i = 1:n
        g(i+1) = g(i) + h * (lambdaD * g(i) * (1 - g(i)) - lambdaA * g(i) - lambdaK * g(i) * f(i));
        f(i+1) = f(i) + h * (-1 * lambdaD * f(i) + p(i * h));
    end
    
    fSolve = f;
    gSolve = g;
end
