using Plots
using LinearAlgebra
using Printf

pyplot()

#This file is used to computer rk4
include("modularRKFUNC.jl")

#File is used to plot the orbital shape
include("plottingOrbitalShapeFunc.jl")

#File is used to define the function that is used in rk4
include("odeFunBlackHole.jl")


function blackHoleOrbit(time, mass1, massRatiom1_m2, xInit, yInit, comX_Init, comY_Init, Gval,PlotDim)

    t0 = t[1]
    tfin = t[end]
    dt = t[2] - t[1]
    loopSize = length(t)

    #The variables defined here are set to const so that they cannot be changed
    #and so that they can be used in the function odefunpp without having to pass
    m1 = mass1
    massRatio = massRatiom1_m2
    x0 = xInit
    y0 = yInit
    comX = comX_Init
    comY = comY_Init
    G = Gval
    t0 = t0
    loopSize = length(t)

    #Use the mass ratio to compute the mass of the second black hole
    m2 = massRatio * m1

    #make the marker size a function of the mass
    marker1 = 3 * m1
    marker2 = 3 * m2

    #Using a if statement to check what planet should be orbiting
    if m1 < m2 
        orbiting_val = 1/m1
    else
        orbiting = 1/m2
    end

    #The value of mG_Var is a constant in the momentum equation so it is defined here
    mG_var = G * m1 * m2

    #Compute the radius of the black from the starting position
    r = norm([x0- comX, y0 - comY])

    #Compute the velocity of the black hole
    v = sqrt(G * (m1 + m2) / r)

    #initial x momentum
    px0 = 0

    #initial y momentum
    py0 = m1 * v 

    u = zeros(4, length(t))
    #u[1] is the x position of the black hole
    #u[2] is the y position of the black hole
    #u[3] is the x momentum of the black hole
    #u[4] is the y momentum of the black hole
    u[:,1] = [x0; y0; px0; py0]

    u0 = u[:,1]
    u = rk4modular(odefunBlackHole, u0, t,orbiting_val, mG_var)

    #Define a new figure to be plotted on
    p = plot()

    ############### Call to plottingOrbitalShapeFunc.jl to plot the the params of the orbit ###############################
    #find largest x value
    x_max = round(maximum(u[1,:]), digits = 0)
    
    #find largest y value
    y_max = round(maximum(u[2,:]), digits = 0) 
    
    x_annotate = x_max
    y_annotate = y_max /2 

    #Add try catch to find the orbit equation
    #try
    #    x_annotate, y_annotate, p =  find_orbit_equation(u,p) 
    #catch
    #    println("The orbit equation could not be found")
    #end

    #######################################################################################################################
    
    #Plot the result of the Runge Kutta method on p 
    plot!(p, u[1,:], u[2,:], label = "Black Hole Trajectory", xlabel = "x position", ylabel = "y position", title = "Black Hole Trajectory")

    #Plot the starting position
    scatter!(p,[x0], [y0], label = "Starting Position", markersize = marker1)


    #Plot the center of mass
    scatter!(p,[comX], [comY], label = "Center of Mass", markersize = marker2)

    #Plot the final position of the black hole
    scatter!(p,[u[1,end]], [u[2,end]], label = "Final Postion", markersize = marker1)

    #move of the page to the right of max x value

    #Moving the plotting legend

    y_max = round(maximum(u[2,:]),digits = 0)
    y_min =round(minimum(u[2,:]), digits = 0)

    x_max = round(maximum(u[1,:]), digits = 0)
    x_min = round(minimum(u[1,:]), digits = 0)

    #Move the legend to the top right of the plot with anchor 
    plot!(p,legend = :topright, legendfont = 10, legendtitle = "Legend",legendtitlefont = 10, legendlabel = ["Black Hole Trajectory", "Starting Position", "Center of Mass", "Final Position"], bbox_to_anchor = (x_max + 10, y_max + 10), size = PlotDim)

    #set the x and y ticks 

    xlims!(p,x_min-1, x_max+ 20)
    ylims!(p,y_min-1, y_max+ 5)

    #dynamic set the x and y ticks
    XtickIncrementer = (x_max - x_min)/ 10

    ytickIncrementer = (y_max- y_min)/10

    x_tick = x_min:XtickIncrementer:y_max
    y_tick = y_min:ytickIncrementer:y_max

    #set the x and y ticks
    xticks!(x_tick)
    yticks!(y_tick)

    function annotateInitConditions(p,x_annotate, y_annotate, r, m1, m2, dt, tfin, x0, y0, pX0, pY0 ,G)

        values = [r, m1, m2, dt, tfin,G]    
        StringLabels = ["r", "m1", "m2", "dt", "tfin", "G"]

        for i in 1:length(values)
            #Round the values to 2 decimal places
            values[i] = round(values[i], digits = 2)

            annotate!(p,x_annotate, y_annotate - 4 * i, text("$(StringLabels[i]) = $(values[i])", :left, 10))
        end
        #Annotate the starting x0, y0 value
        #Move x0 the right .2 units
        x0Plot = x0 + .2

        annotate!(p,[x0Plot], [y0], text("($x0,$y0)" , :left, 10))

        #Annotate the starting px0, py0 value
        round_px0 = round(pX0, digits = 0)
        round_py0 = round(pY0, digits = 2)
        annotate!(p,x_annotate, y_annotate - 4 * length(values)-2, text("Initial P ($round_px0,$round_py0)" , :left, 10))
    end

    #Change the size of the plot
    plot!(p, size = (plotDim[1],plotDim[2]))

    annotateInitConditions(p,x_annotate, y_annotate, r, m1, m2, dt, tfin, x0, y0, px0, py0,G)

    #Round the values of the mass, G, and r
    m1 = round(m1, digits = 0)
    m2 = round(m2, digits = 0)
    G = round(G, digits = 0)
    r = round(r, digits = 0)

    #Write dt in scientific notation using @sprintf
    dt = @sprintf("%.e0", dt)

    #write tfin in scientific notation
    tfin = @sprintf("%.0e", tfin)

    py0 = round(py0, digits = 0)

    #Create unique file name based on m1, m2, G, and r
    filename = "BlackHole_m1_$(m1)_m2_$(m2)_Py0_$(py0)_G$(G)_r$(r)_dt$(dt)_tfin$(tfin).png"

    #Save the plot as a png file
    savefig(p, "/Users/gagexavieranastasio/Desktop/UMassD_SU2023_REU/UMassD_S23_Julia/WK4/momentum_profScott_prob/BlackHolePlots/TestingCond1/$filename")
end