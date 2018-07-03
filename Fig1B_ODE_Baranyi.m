function [t, X] = ODE_Baranyi(t_end,dt,n_ode)
%Solves a system of coupled ODEs via MATLAB solver 
[t X] = ode45(@equations,0:dt:t_end,n_ode); 

end

function dx = equations(t,x)

    global mu N_max lambda n
    
    dx = zeros(1,1);
    
    dx(1) = mu*x(1)*(t^n/(lambda^n+t^n))*(1-(x(1)/N_max));
    
end