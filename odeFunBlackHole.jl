
using LinearAlgebra

function odefunBlackHole(u::Vector{Float64}, t::Float64, orbiting_val, mG_var)
    #A is 4 x 4 matrix 
    #Define A to be a 4 x 4 matrix
    A = zeros(4,4)
    #x(t)
    A[1,3] = orbiting_val

    #y(t)
    A[2,4] = orbiting_val

    #radius cubed
    r_3 = norm([u[1], u[2]])^3
 
    #px(t)
    A[3,1] = -mG_var / (r_3) 

    #py(t)
    A[4,2] = -mG_var / (r_3) 

    #Define the derivative of the system of equations
    dydt = A * u 
    return dydt 
end

