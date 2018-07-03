function [t, X] = Fig1A_ODE_solver(t_end,dt,N_Allee)
%Solves a system of coupled ODEs via MATLAB solver 
[t X] = ode45(@equations,0:dt:t_end,N_Allee);

end

function dx = equations(t,x)

    global r K Nc
    
    dx = zeros(1,1);
 
    dx(1) = r*x(1)*(1-(x(1)/K))*((x(1)/Nc)-1); 
    
end