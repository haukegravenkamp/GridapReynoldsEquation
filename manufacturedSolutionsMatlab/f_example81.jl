# copied symbolic expression from Matlab file
fm(x,y) = ((0.5*cos(x) + 1)^3*(0.5*pi^2*cos(pi*y)*sin(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(2*x)/3 - 1/3) - (25.0*pi*sin(pi*y)^2*sin(x)^2*(cos(2*x)/3 - 1/3)^2)/(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1) + (12.5*pi*cos(pi*y)*sin(x)^2*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3)^2)/(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1) + (15625.0*pi*sin(pi*y)^2*sin(x)^4*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^4)/(625.0*sin(x)^2*(cos(y*pi) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1)^2))/12 - 0.5*sin(x) + ((0.5*cos(x) + 1)^3*(0.66666666666666666666666666666667*cos(2*x)*sin(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(pi*y) + 1) + 0.66666666666666666666666666666667*sin(2*x)*cos(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(pi*y) + 1) + 0.5*sin(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) + (0.66666666666666666666666666666667*sin(2*x)*sin(x)*(25.0*cos(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 16.666666666666666666666666666667*sin(2*x)*sin(x)*(cos(pi*y) + 1))*(cos(pi*y) + 1))/(pi*(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1)) + (0.5*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3)*(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) + 33.333333333333333333333333333333*cos(2*x)*sin(x)*(cos(pi*y) + 1) + 33.333333333333333333333333333333*sin(2*x)*cos(x)*(cos(pi*y) + 1)))/(pi*(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1)) - (1.0*cos(x)*(25.0*cos(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 16.666666666666666666666666666667*sin(2*x)*sin(x)*(cos(pi*y) + 1))*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/(pi*(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1)) - (0.5*sin(x)*(833.33333333333333333333333333333*sin(2*x)*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3) - 1250.0*cos(x)*sin(x)*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2)*(25.0*cos(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 16.666666666666666666666666666667*sin(2*x)*sin(x)*(cos(pi*y) + 1))*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/(pi*(625.0*sin(x)^2*(cos(y*pi) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1)^2)))/12 + 0.125*sin(x)*(0.5*cos(x) + 1)^2*(0.5*cos(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 0.33333333333333333333333333333333*sin(2*x)*sin(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(pi*y) + 1) + (0.5*sin(x)*(25.0*cos(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 16.666666666666666666666666666667*sin(2*x)*sin(x)*(cos(pi*y) + 1))*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/(pi*(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1))) + 0.25*sin(x)^2*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi + 0.5)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 0.5*cos(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi + 0.5)*(0.5*cos(x) + 1)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) + 0.33333333333333333333333333333333*sin(2*x)*sin(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi + 0.5)*(0.5*cos(x) + 1)*(cos(pi*y) + 1) - (0.5*sin(x)*(25.0*cos(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 16.666666666666666666666666666667*sin(2*x)*sin(x)*(cos(pi*y) + 1))*(0.5*cos(x) + 1)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/(pi*(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1))
