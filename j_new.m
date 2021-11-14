%inputs:
%delta_vec=  ; a=0.0005; b=0.00005; c=0.001; T=50;
%g_vec and f_vec computed in main


function j_new = j_new(delta_tau_vec, params, tum)

[f, g, p, t] = tum(delta_tau_vec);

delta_vec = delta_tau_vec(1:4);

toxic_effects=sum(delta_vec.^2);                %first term of Jnew
cellPopulation1=params.a*g(end);                   %second term - index last time value (T)
cellPopulation2=params.b*trap_int(params.dt, g);          %third term
cellPopulation3=params.c*trap_int(params.dt,f);           %fourth term

j_new=toxic_effects+cellPopulation1+cellPopulation2+cellPopulation3;

end
