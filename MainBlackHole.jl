

include("BlackHoleFunc.jl")


using Plots
using LinearAlgebra
pyplot()

#Define the time vector
t0 = 0
tfin = 10000
dt = .0001
t = t0:dt:tfin

#Define mass of the black hole
mass1 = 1

#Define the mass ration
#   m2 = massRatio * m1
massRatiom1_m2 = 10

#Starting position of the orbiting black hole m1
y0 = 1#xInit
x0 = 1 #yInit

#Position of the center of mass
comX = 0 #comX_Init
comY = 0 #comY_Init

#Gravitational constant
G = 1 #Gval

#Define the plot dimensions
plotDim = [1000, 1000]

#Call the function blackHoleOrbit function from BlackHoleFunc.jl
#blackHoleOrbit(t, mass1, massRatiom1_m2, x0, y0, comX, comY, G, plotDim)


#G starts at 50 and goes down by 10 to 1
G = [1, 2 ,3]

for i in 1:length(G)
    run(`say \"I am on $i\"`)
    println("I am on $i")
    blackHoleOrbit(t, mass1, massRatiom1_m2, x0, y0, comX, comY, G[i], plotDim)
end


run(`say \"I am done\"`)