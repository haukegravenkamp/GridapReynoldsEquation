## manufactured solutions used in Sections 8.1 and 8.2 in the paper

# The chosen solutions are given by u_ex1 and u_ex2 below. The corresponding body loads have been computed using the Matlab script available at  https://github.com/haukegravenkamp/manufacturedSolutions and given by the functions f_ex1.jl and f_ex2.jl 

export f_ex1, u_ex1, selectExample

function selectExample(exampleNo)

    if exampleNo == 1
        return u_ex1, f_ex1, "./results/example81/", ub, u0
    elseif exampleNo == 2
        return u_ex2, f_ex2, "./results/example82/", ub, u0
    elseif exampleNo == 3
        return u_ex3, f_ex3, "./results/example83/", ub, u0
    end
end

function f_ex1(xy::VectorValue{2,Float64})::Float64

    x = xy[1]
    y = xy[2]
    ((0.5 * cos(x) + 1)^3 * (0.5 * pi^2 * cos(pi * y) * sin(x) * (atan(25.0 * sin(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / pi - 0.5) * (cos(2 * x) / 3 - 1 / 3) - (25.0 * pi * sin(pi * y)^2 * sin(x)^2 * (cos(2 * x) / 3 - 1 / 3)^2) / (625.0 * sin(x)^2 * (cos(pi * y) + 1)^2 * (cos(2 * x) / 3 - 1 / 3)^2 + 1) + (12.5 * pi * cos(pi * y) * sin(x)^2 * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)^2) / (625.0 * sin(x)^2 * (cos(pi * y) + 1)^2 * (cos(2 * x) / 3 - 1 / 3)^2 + 1) + (15625.0 * pi * sin(pi * y)^2 * sin(x)^4 * (cos(pi * y) + 1)^2 * (cos(2 * x) / 3 - 1 / 3)^4) / (625.0 * sin(x)^2 * (cos(y * pi) + 1)^2 * (cos(2 * x) / 3 - 1 / 3)^2 + 1)^2)) / 12 - 0.5 * sin(x) + ((0.5 * cos(x) + 1)^3 * (0.66666666666666666666666666666667 * cos(2 * x) * sin(x) * (atan(25.0 * sin(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / pi - 0.5) * (cos(pi * y) + 1) + 0.66666666666666666666666666666667 * sin(2 * x) * cos(x) * (atan(25.0 * sin(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / pi - 0.5) * (cos(pi * y) + 1) + 0.5 * sin(x) * (atan(25.0 * sin(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / pi - 0.5) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3) + (0.66666666666666666666666666666667 * sin(2 * x) * sin(x) * (25.0 * cos(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3) - 16.666666666666666666666666666667 * sin(2 * x) * sin(x) * (cos(pi * y) + 1)) * (cos(pi * y) + 1)) / (pi * (625.0 * sin(x)^2 * (cos(pi * y) + 1)^2 * (cos(2 * x) / 3 - 1 / 3)^2 + 1)) + (0.5 * sin(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3) * (25.0 * sin(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3) + 33.333333333333333333333333333333 * cos(2 * x) * sin(x) * (cos(pi * y) + 1) + 33.333333333333333333333333333333 * sin(2 * x) * cos(x) * (cos(pi * y) + 1))) / (pi * (625.0 * sin(x)^2 * (cos(pi * y) + 1)^2 * (cos(2 * x) / 3 - 1 / 3)^2 + 1)) - (1.0 * cos(x) * (25.0 * cos(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3) - 16.666666666666666666666666666667 * sin(2 * x) * sin(x) * (cos(pi * y) + 1)) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / (pi * (625.0 * sin(x)^2 * (cos(pi * y) + 1)^2 * (cos(2 * x) / 3 - 1 / 3)^2 + 1)) - (0.5 * sin(x) * (833.33333333333333333333333333333 * sin(2 * x) * sin(x)^2 * (cos(pi * y) + 1)^2 * (cos(2 * x) / 3 - 1 / 3) - 1250.0 * cos(x) * sin(x) * (cos(pi * y) + 1)^2 * (cos(2 * x) / 3 - 1 / 3)^2) * (25.0 * cos(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3) - 16.666666666666666666666666666667 * sin(2 * x) * sin(x) * (cos(pi * y) + 1)) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / (pi * (625.0 * sin(x)^2 * (cos(y * pi) + 1)^2 * (cos(2 * x) / 3 - 1 / 3)^2 + 1)^2))) / 12 + 0.125 * sin(x) * (0.5 * cos(x) + 1)^2 * (0.5 * cos(x) * (atan(25.0 * sin(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / pi - 0.5) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3) - 0.33333333333333333333333333333333 * sin(2 * x) * sin(x) * (atan(25.0 * sin(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / pi - 0.5) * (cos(pi * y) + 1) + (0.5 * sin(x) * (25.0 * cos(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3) - 16.666666666666666666666666666667 * sin(2 * x) * sin(x) * (cos(pi * y) + 1)) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / (pi * (625.0 * sin(x)^2 * (cos(pi * y) + 1)^2 * (cos(2 * x) / 3 - 1 / 3)^2 + 1))) + 0.25 * sin(x)^2 * (atan(25.0 * sin(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / pi + 0.5) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3) - 0.5 * cos(x) * (atan(25.0 * sin(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / pi + 0.5) * (0.5 * cos(x) + 1) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3) + 0.33333333333333333333333333333333 * sin(2 * x) * sin(x) * (atan(25.0 * sin(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / pi + 0.5) * (0.5 * cos(x) + 1) * (cos(pi * y) + 1) - (0.5 * sin(x) * (25.0 * cos(x) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3) - 16.666666666666666666666666666667 * sin(2 * x) * sin(x) * (cos(pi * y) + 1)) * (0.5 * cos(x) + 1) * (cos(pi * y) + 1) * (cos(2 * x) / 3 - 1 / 3)) / (pi * (625.0 * sin(x)^2 * (cos(pi * y) + 1)^2 * (cos(2 * x) / 3 - 1 / 3)^2 + 1))

end

function u_ex1(x::VectorValue{2,Float64})::Float64
    1 / 6 * (1 - cos(2x[1])) * sin(x[1]) * (1 + cos(π * x[2]))
end

function f_ex2(xy::VectorValue{2,Float64})::Float64

    x = xy[1]
    y = xy[2]

    ((0.5 * cos(x) + 1)^3 * ((atan(50.0 * (cos(pi * y) + 1) * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)) / pi + 1 / 2) * (0.03125 * cos(x / 2) - (625.0 * exp((50 * x) / pi)) / (pi^2 * (exp(100) - 1))) * (cos(pi * y) + 1) - (100.0 * (0.0625 * sin(x / 2) - (12.5 * exp((50 * x) / pi)) / (pi * (exp(100) - 1)))^2 * (cos(pi * y) + 1)^2) / (pi * (2500.0 * (cos(pi * y) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)^2 + 1)) + (250000.0 * (0.0625 * sin(x / 2) - (12.5 * exp((50 * x) / pi)) / (pi * (exp(100) - 1)))^2 * (cos(pi * y) + 1)^4 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)^2) / (pi * (2500.0 * (cos(y * pi) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * exp((50 * x) / pi) - 0.25) / (exp(100) - 1) - 0.125)^2 + 1)^2) + (50.0 * (0.03125 * cos(x / 2) - (625.0 * exp((50 * x) / pi)) / (pi^2 * (exp(100) - 1))) * (cos(pi * y) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)) / (pi * (2500.0 * (cos(pi * y) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)^2 + 1)))) / 12 - 0.5 * sin(x) + ((0.5 * cos(x) + 1)^3 * (pi^2 * cos(pi * y) * (atan(50.0 * (cos(pi * y) + 1) * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)) / pi + 1 / 2) * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125) - (100.0 * pi * sin(pi * y)^2 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)^2) / (2500.0 * (cos(pi * y) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)^2 + 1) + (250000.0 * pi * sin(pi * y)^2 * (cos(pi * y) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)^4) / (2500.0 * (cos(y * pi) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * exp((50 * x) / pi) - 0.25) / (exp(100) - 1) - 0.125)^2 + 1)^2 + (50.0 * pi * cos(pi * y) * (cos(pi * y) + 1) * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)^2) / (2500.0 * (cos(pi * y) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)^2 + 1))) / 12 - 0.125 * sin(x) * ((atan(50.0 * (cos(pi * y) + 1) * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)) / pi + 1 / 2) * (0.0625 * sin(x / 2) - (12.5 * exp((50 * x) / pi)) / (pi * (exp(100) - 1))) * (cos(pi * y) + 1) + (50.0 * (0.0625 * sin(x / 2) - (12.5 * exp((50 * x) / pi)) / (pi * (exp(100) - 1))) * (cos(pi * y) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)) / (pi * (2500.0 * (cos(pi * y) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)^2 + 1))) * (0.5 * cos(x) + 1)^2 + 0.5 * sin(x) * (atan(50.0 * (cos(pi * y) + 1) * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)) / pi - 1 / 2) * (cos(pi * y) + 1) * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125) + (atan(50.0 * (cos(pi * y) + 1) * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)) / pi - 1 / 2) * (0.0625 * sin(x / 2) - (12.5 * exp((50 * x) / pi)) / (pi * (exp(100) - 1))) * (0.5 * cos(x) + 1) * (cos(pi * y) + 1) + (50.0 * (0.0625 * sin(x / 2) - (12.5 * exp((50 * x) / pi)) / (pi * (exp(100) - 1))) * (0.5 * cos(x) + 1) * (cos(pi * y) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)) / (pi * (2500.0 * (cos(pi * y) + 1)^2 * (0.125 * cos(x / 2) + (0.25 * (exp((50 * x) / pi) - 1)) / (exp(100) - 1) - 0.125)^2 + 1))

end

function u_ex2(x::VectorValue{2,Float64})::Float64
    0.25 * ((1 - exp(100 * x[1] / (2 * π))) / (1 - exp(100)) - 1 + (cos(x[1] / 2) + 1) / 2) * (1 + cos(π * x[2]))
end



function f_ex3(xy::VectorValue{2,Float64})::Float64
    0.0
end

function u_ex3(xy::VectorValue{2,Float64})::Float64
    0.0
end


function u0(xy::VectorValue{2,Float64})::Float64
    1.0
end

function ub(xy::VectorValue{2,Float64})::Float64
    0.0
end