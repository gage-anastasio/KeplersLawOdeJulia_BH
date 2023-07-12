

#A modular rk4 class that works for both non systems and systems of equations
#input: odefun, dt, u, t, a , b
#odefun is the ode function
#dt is the time step
#u is the solution vector --- u[1] is the initial condition
#t is the time vector
#a and b are the parameters

#output: u

#The class is called by the function rk4modular


function rk4modular(f, u0, t,orbiting_val, mG_var)
    #If u0 is a vector, then the system of equations is solved
    if isa(u0, Vector)
        dt = t[2] - t[1]
        u = zeros(length(u0), length(t))
        u[:,1] = u0
        for i = 2:length(t) 
            k1 = dt * f(u[:,i-1], t[i-1],orbiting_val, mG_var)
            k2 = dt * f(u[:,i-1] .+ (dt * k1)/2, t[i-1] + dt/2,orbiting_val, mG_var)
            k3 = dt * f(u[:,i-1] .+ (dt * k2)/2, t[i-1] + dt/2,orbiting_val, mG_var)
            k4 = dt * f(u[:,i-1] .+ (dt * k3), t[i-1] + dt,orbiting_val, mG_var)
            u[:,i] = u[:,i-1] + (k1 .+ 2*k2 .+ 2*k3 .+ k4) /6 
        end
        return u

    #If u0 is a float, then the non system of equations is solved
    elseif isa(u0, Float64)
        dt = t[2] - t[1]
        u = zeros(length(t))
        u[1] = u0
        for i = 1:length(t)-1
            k1 = dt * f(u[i], t[i])
            k2 = dt * f(u[i] + k1/2, t[i] + dt/2)
            k3 = dt * f(u[i] + k2/2, t[i] + dt/2)
            k4 = dt * f(u[i] + k3, t[i] + dt)
            u[i+1] = u[i] + (k1 + 2*k2 + 2*k3 + k4)/6
        end
        return u
    end

end



