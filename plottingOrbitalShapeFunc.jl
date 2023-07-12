
#Define a function that will find an equation for the elipse shape of the orbit
function find_orbit_equation(u,p)
    #Use the only the first lap of the orbit to find the equation
    #Find the index of the first lap

    #Define x0
    x0 = u[1,1]

    #Define y0
    y0 = u[2,1]

    tol = 0.008
    #Find when planet return to the starting position x0
    index = findall(isapprox.(x0, u[1, :], atol=tol) .& isapprox.(y0, u[2, :], atol=tol))
    
    function find_first_large_gap(index)
        #Find the first large gap in the index
        for i in 1:length(index)-1
            if index[i+1] - index[i] > 1000
                return [index[i], index[i+1]]
            end
        end
    end

    largest_gap_index = find_first_large_gap(index)

    #function to find the slope of of the axis of the elipse
    function find_slope(xData, yData)

        #Find the smallestand largest x position 
        x_min = minimum(xData)
        x_max = maximum(xData)

        #Find the corresponding y values @ x_min and x_max
        y_min = yData[findall(x_min .== xData)]
        y_max = yData[findall(x_max .== xData)]

        #Plot Y min @ x min
        #scatter!([x_min], [y_min], label = "Y min @ x min", markersize = 5)

        #Use linear regression to find the slope of the axis of the elipse
        m = (y_max - y_min) / (x_max - x_min)

        # The equation for the line is y - y_min = m(x - x_min)
        #Find the y intercept
        b = y_min - m * x_min

        return m, b
    end

    #Slice the u matrix to only include x and y position of the first complete lap
    u_slice = u[1:2, 1:largest_gap_index[2]]

    #Make the last part of u slice equal to the first value of u slice
    u_slice[:,end] = u_slice[:,1]

    #Define the x and y data from the first lap of the orbit as a vector 
    xData = vec(u_slice[1,:]')
    yData = vec(u_slice[2,:]')

    #Find the slope of the axis of the elipse
    m, b = find_slope(xData, yData)

    #Find the of the line
    function find_y(x, m, b)
        return m .* x .+ b
    end

    
    theta = atan(m[1])

    #Find the y values of the line
    y = find_y(xData, m, b)

    #Plot the first point of y and xData
    #scatter!([xData[1]], [y[1]], label = "First Point", markersize = 5)


    #Return the x and y data and the slope and y intercept of the line from 0 to xmin to xmax

    #Plot the x and y data
    plot!(p, xData, y, label = "Orbit Data", linewidth = 2)

    #Find the midpoint between xData and y
    function findMidpoint(x, y)
        xMax = maximum(x)
        xMin = minimum(x)
        ymax = maximum(y)
        ymin = minimum(y)

        x_mid = (xMax + xMin) / 2
        y_mid = (ymax + ymin) / 2

        #println("\nThe midpoint is: ", [x_mid, y_mid])
        
        #find the distance from the midpoint to the xMax and yMax using norm
        midpoint = [x_mid, y_mid]

        #Find the poinit where x_mid and y_mid intersect with the elipse
        

        return x_mid, y_mid
    end

    #Find the line that is perp to the line of the axis of the elipse at the midpoint
    function find_perp_line(xData, x_mid, y_mid, m,yData)

        function find_perp_slope(m)
            if m != 0
                return -1 / m
            else
                return 0
            end
        end

        #Find the slope of the line that is perp to the line of the axis of the elipse at the midpoint
        m_perp = find_perp_slope(m[1])

        #Find the y intercept of the line that is perp to the line of the axis of the elipse at the midpoint
        b_perp = y_mid - m_perp * x_mid

        #Max x and y
        x_min = minimum(xData)
        x_max = maximum(xData)
        y_min = minimum(yData)
        y_max = maximum(yData)
        
        #Slice the values of xData and yData so it only plots within the elipse
        #xData = xData[findall(x_min .<= xData .<= x_max)]
        
        #Find the y values of the line that is perp to the line of the axis of the elipse at the midpoint
        y_perp = find_y(xData, m_perp, b_perp)

        #Find when y_perp intersects with the top of the elipse
        index = findall(isapprox.(y_perp, yData, atol=0.008))

    
        #Find when [xData, y_perp] approx interscts with the elipse [xData, yData] at the bottom
        #index_bottom = findall(isapprox.(y_perp, yData, atol=0.008))
        
        sliceVars = find_first_large_gap(index)

        #println("\nSliceVars: ", sliceVars)

        #find when y_perp first equals yMid
        index2 = findall(isapprox.(y_perp, y_mid, atol=0.008))


        #Slice the xData and y_perp to only include the values between the sliceVars
        #y_perp = y_perp[sliceVars_bottom[1]:sliceVars[2]]
        y_perp = y_perp[index2[1]:1000]


        #create a function to find the x vales from y_perp
        function find_x(y_perp, m_perp, b_perp)
            return (y_perp .- b_perp) ./ m_perp
        end

        #Find the x values from y_perp
        x_perp = find_x(y_perp, m_perp, b_perp)

        #The point at which the perp line intersects the top of the elipse
        topPoint = [xData[index[1]], yData[index[1]]]
        midpoint = [x_mid, y_mid]

        #Find the distance between the top point and the midpoint
        distanceB = norm(topPoint .- midpoint)
        distanceB = round(distanceB, digits = 2)

        #Add text to the plot to display the distance B
        
      
        #Plot point where y_perp intersects with the bottom of the elipse
        #plot!([xData[index_bottom[1]]], [yData[index_bottom[1]]], label = "Intersection", markersize = 5)

        return distanceB
    end


    x_mid, y_mid = findMidpoint(xData, yData)

    #Find the distance A from the midpoint to xData[1], y[1] using norm
    distanceA = [xData[1], y[1]] - [x_mid, y_mid]
    distanceA = norm(distanceA)
    distanceA = round(distanceA, digits = 2)

    #Plot the midpoint
    scatter!(p,[x_mid], [y_mid], label = "Midpoint", markersize = 5)

    #Find the y values of the line that is perp to the line of the axis of the elipse at the midpoint
    distanceB = find_perp_line(xData, x_mid, y_mid, m, yData)

    x_max = maximum(xData)

    annotate!(p,x_max + 10, y_mid, text("B = $distanceB", :left, 10))


    annotate!(p,x_max + 10, y_mid - 2, text("A = $distanceA", :left, 10))

    #Round theta to 2 digits
    theta = round(theta, digits = 2)

    annotate!(p,x_max + 10, y_mid - 4, text("θ = $theta rads", :left, 10))

    return x_max + 10, y_mid - 4, p
end