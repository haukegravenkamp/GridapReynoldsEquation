f = @(x,y)sin(x).*-5.0e-1+((cos(x).*5.0e-1+1.0).^3.*(pi.^2.*cos(y.*pi).*sin(x).*(atan(sin(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1)./pi-5.0e-1).*(cos(x.*2.0)./3.0-1.0./3.0).*5.0e-1-(pi.*sin(y.*pi).^2.*sin(x).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^2.*2.5e+1)./(sin(x).^2.*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^2.*6.25e+2+1.0)+(pi.*cos(y.*pi).*sin(x).^2.*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).^2.*1.25e+1)./(sin(x).^2.*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^2.*6.25e+2+1.0)+pi.*sin(y.*pi).^2.*sin(x).^4.*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^4.*1.0./(sin(x).^2.*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^2.*6.25e+2+1.0).^2.*1.5625e+4))./1.2e+1+((cos(x).*5.0e-1+1.0).^3.*(cos(x.*2.0).*sin(x).*(atan(sin(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1)./pi-5.0e-1).*(cos(y.*pi)+1.0).*6.666666666666667e-1+sin(x.*2.0).*cos(x).*(atan(sin(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1)./pi-5.0e-1).*(cos(y.*pi)+1.0).*6.666666666666667e-1+sin(x).*(atan(sin(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1)./pi-5.0e-1).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*5.0e-1+(sin(x.*2.0).*sin(x).*(cos(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1-sin(x.*2.0).*sin(x).*(cos(y.*pi)+1.0).*1.666666666666667e+1).*(cos(y.*pi)+1.0).*6.666666666666667e-1)./(pi.*(sin(x).^2.*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^2.*6.25e+2+1.0))+(sin(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*(sin(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1+cos(x.*2.0).*sin(x).*(cos(y.*pi)+1.0).*3.333333333333333e+1+sin(x.*2.0).*cos(x).*(cos(y.*pi)+1.0).*3.333333333333333e+1).*5.0e-1)./(pi.*(sin(x).^2.*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^2.*6.25e+2+1.0))-(cos(x).*(cos(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1-sin(x.*2.0).*sin(x).*(cos(y.*pi)+1.0).*1.666666666666667e+1).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*1.0)./(pi.*(sin(x).^2.*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^2.*6.25e+2+1.0))-(sin(x).*(sin(x.*2.0).*sin(x).^2.*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).*8.333333333333333e+2-cos(x).*sin(x).*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^2.*1.25e+3).*(cos(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1-sin(x.*2.0).*sin(x).*(cos(y.*pi)+1.0).*1.666666666666667e+1).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*1.0./(sin(x).^2.*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^2.*6.25e+2+1.0).^2.*5.0e-1)./pi))./1.2e+1+sin(x).*(cos(x).*5.0e-1+1.0).^2.*(sin(x.*2.0).*sin(x).*(atan(sin(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1)./pi-5.0e-1).*(cos(y.*pi)+1.0).*-3.333333333333333e-1+cos(x).*(atan(sin(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1)./pi-5.0e-1).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*5.0e-1+(sin(x).*(cos(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1-sin(x.*2.0).*sin(x).*(cos(y.*pi)+1.0).*1.666666666666667e+1).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*5.0e-1)./(pi.*(sin(x).^2.*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^2.*6.25e+2+1.0))).*1.25e-1+sin(x).^2.*(atan(sin(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1)./pi+5.0e-1).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e-1-cos(x).*(atan(sin(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1)./pi+5.0e-1).*(cos(x).*5.0e-1+1.0).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*5.0e-1+sin(x.*2.0).*sin(x).*(atan(sin(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1)./pi+5.0e-1).*(cos(x).*5.0e-1+1.0).*(cos(y.*pi)+1.0).*3.333333333333333e-1-(sin(x).*(cos(x).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*2.5e+1-sin(x.*2.0).*sin(x).*(cos(y.*pi)+1.0).*1.666666666666667e+1).*(cos(x).*5.0e-1+1.0).*(cos(y.*pi)+1.0).*(cos(x.*2.0)./3.0-1.0./3.0).*5.0e-1)./(pi.*(sin(x).^2.*(cos(y.*pi)+1.0).^2.*(cos(x.*2.0)./3.0-1.0./3.0).^2.*6.25e+2+1.0));

f = @(x,y)sin(x).*-5.0e-1+((cos(x)+2.0).^3.*(sin(x).*(atan(sin(x).*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*1.666666666666667e+1).*3.183098861837907e-1-5.0e-1).*(cos(y.*pi)+1.0).*(cos(x).^2.*2.0-1.0).*6.666666666666667e-1+sin(x).*(atan(sin(x).*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*1.666666666666667e+1).*3.183098861837907e-1-5.0e-1).*(cos(y.*pi)+1.0).*(cos(x).^2.*(2.0./3.0)-2.0./3.0).*5.0e-1+cos(x).^2.*sin(x).*(atan(sin(x).*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*1.666666666666667e+1).*3.183098861837907e-1-5.0e-1).*(cos(y.*pi)+1.0).*1.333333333333333+(sin(x).^2.*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2.*3.0-1.0).*(cos(y.*pi)+1.0).*(cos(x).^2.*(2.0./3.0)-2.0./3.0).*2.25e+2)./(pi.*(sin(x).^2.*(cos(x).^2-1.0).^2.*(cos(y.*pi)+1.0).^2.*2.5e+3+9.0))-(cos(x).^2.*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*(cos(y.*pi)+1.0).*(cos(x).^2.*(2.0./3.0)-2.0./3.0).*4.5e+2)./(pi.*(sin(x).^2.*(cos(x).^2-1.0).^2.*(cos(y.*pi)+1.0).^2.*2.5e+3+9.0))+(cos(x).^2.*sin(x).^2.*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*(cos(y.*pi)+1.0).*6.0e+2)./(pi.*(sin(x).^2.*(cos(x).^2-1.0).^2.*(cos(y.*pi)+1.0).^2.*2.5e+3+9.0))-(cos(x).^2.*sin(x).^4.*(cos(y.*3.141592653589793)+1.0).^3.*1.0./(sin(x).^2.*(cos(x).^2-1.0).^2.*(cos(y.*pi)+1.0).^2.*2.5e+3+9.0).^2.*(cos(x).^2-1.0).^2.*(cos(y.*pi)+1.0).*(cos(x).^2.*(2.0./3.0)-2.0./3.0).*3.375e+6)./pi))./9.6e+1+((cos(x)+2.0).^3.*(pi.^2.*cos(y.*pi).*sin(x).*(atan(sin(x).*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*1.666666666666667e+1).*3.183098861837907e-1-5.0e-1).*(cos(x).^2.*(2.0./3.0)-2.0./3.0).*5.0e-1-(pi.*sin(y.*pi).^2.*sin(x).^2.*(cos(x).^2-1.0).^2.*1.0e+2)./(sin(x).^2.*(cos(x).^2-1.0).^2.*(cos(y.*pi)+1.0).^2.*2.5e+3+9.0)+pi.*sin(y.*pi).^2.*sin(x).^4.*(cos(x).^2-1.0).^4.*1.0./(sin(x).^2.*(cos(x).^2-1.0).^2.*(cos(y.*pi)+1.0).^2.*2.5e+3+9.0).^2.*(cos(y.*pi)+1.0).^2.*2.5e+5+(pi.*cos(y.*pi).*sin(x).^2.*(cos(x).^2-1.0).^2.*(cos(y.*pi)+1.0).*5.0e+1)./(sin(x).^2.*(cos(x).^2-1.0).^2.*(cos(y.*pi)+1.0).^2.*2.5e+3+9.0)))./9.6e+1+sin(x).*(cos(x)+2.0).^2.*(cos(x).*(atan(sin(x).*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*1.666666666666667e+1).*3.183098861837907e-1-5.0e-1).*(cos(y.*pi)+1.0).*(cos(x).^2.*(2.0./3.0)-2.0./3.0).*5.0e-1-cos(x).*sin(x).^2.*(atan(sin(x).*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*1.666666666666667e+1).*3.183098861837907e-1-5.0e-1).*(cos(y.*pi)+1.0).*6.666666666666667e-1+(cos(x).*sin(x).*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*(cos(y.*pi)+1.0).*(cos(x).^2.*(2.0./3.0)-2.0./3.0).*2.25e+2)./(pi.*(sin(x).^2.*(cos(x).^2-1.0).^2.*(cos(y.*pi)+1.0).^2.*2.5e+3+9.0))).*3.125e-2+sin(x).^2.*(atan(sin(x).*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*1.666666666666667e+1).*3.183098861837907e-1+5.0e-1).*(cos(y.*pi)+1.0).*(cos(x).^2.*(2.0./3.0)-2.0./3.0).*2.5e-1-cos(x).*(atan(sin(x).*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*1.666666666666667e+1).*3.183098861837907e-1+5.0e-1).*(cos(x)./2.0+1.0).*(cos(y.*pi)+1.0).*(cos(x).^2.*(2.0./3.0)-2.0./3.0).*5.0e-1+cos(x).*sin(x).^2.*(atan(sin(x).*(cos(y.*3.141592653589793)+1.0).*(cos(x).^2-1.0).*1.666666666666667e+1).*3.183098861837907e-1+5.0e-1).*(cos(x)./2.0+1.0).*(cos(y.*pi)+1.0).*6.666666666666667e-1-(cos(x).*sin(x).*(cos(y.*3.141592653589793)+1.0).*(cos(x)./2.0+1.0).*(cos(x).^2-1.0).*(cos(y.*pi)+1.0).*(cos(x).^2.*(2.0./3.0)-2.0./3.0).*2.25e+2)./(pi.*(sin(x).^2.*(cos(x).^2-1.0).^2.*(cos(y.*pi)+1.0).^2.*2.5e+3+9.0))

% symbolic expression from Matlab
f = ((0.5*cos(x) + 1)^3*(0.5*pi^2*cos(pi*y)*sin(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(2*x)/3 - 1/3) - (25.0*pi*sin(pi*y)^2*sin(x)^2*(cos(2*x)/3 - 1/3)^2)/(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1) + (12.5*pi*cos(pi*y)*sin(x)^2*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3)^2)/(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1) + (15625.0*pi*sin(pi*y)^2*sin(x)^4*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^4)/(625.0*sin(x)^2*(cos(y*pi) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1)^2))/12 - 0.5*sin(x) + ((0.5*cos(x) + 1)^3*(0.66666666666666666666666666666667*cos(2*x)*sin(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(pi*y) + 1) + 0.66666666666666666666666666666667*sin(2*x)*cos(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(pi*y) + 1) + 0.5*sin(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) + (0.66666666666666666666666666666667*sin(2*x)*sin(x)*(25.0*cos(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 16.666666666666666666666666666667*sin(2*x)*sin(x)*(cos(pi*y) + 1))*(cos(pi*y) + 1))/(pi*(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1)) + (0.5*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3)*(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) + 33.333333333333333333333333333333*cos(2*x)*sin(x)*(cos(pi*y) + 1) + 33.333333333333333333333333333333*sin(2*x)*cos(x)*(cos(pi*y) + 1)))/(pi*(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1)) - (1.0*cos(x)*(25.0*cos(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 16.666666666666666666666666666667*sin(2*x)*sin(x)*(cos(pi*y) + 1))*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/(pi*(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1)) - (0.5*sin(x)*(833.33333333333333333333333333333*sin(2*x)*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3) - 1250.0*cos(x)*sin(x)*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2)*(25.0*cos(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 16.666666666666666666666666666667*sin(2*x)*sin(x)*(cos(pi*y) + 1))*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/(pi*(625.0*sin(x)^2*(cos(y*pi) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1)^2)))/12 + 0.125*sin(x)*(0.5*cos(x) + 1)^2*(0.5*cos(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 0.33333333333333333333333333333333*sin(2*x)*sin(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi - 0.5)*(cos(pi*y) + 1) + (0.5*sin(x)*(25.0*cos(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 16.666666666666666666666666666667*sin(2*x)*sin(x)*(cos(pi*y) + 1))*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/(pi*(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1))) + 0.25*sin(x)^2*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi + 0.5)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 0.5*cos(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi + 0.5)*(0.5*cos(x) + 1)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) + 0.33333333333333333333333333333333*sin(2*x)*sin(x)*(atan(25.0*sin(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/pi + 0.5)*(0.5*cos(x) + 1)*(cos(pi*y) + 1) - (0.5*sin(x)*(25.0*cos(x)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3) - 16.666666666666666666666666666667*sin(2*x)*sin(x)*(cos(pi*y) + 1))*(0.5*cos(x) + 1)*(cos(pi*y) + 1)*(cos(2*x)/3 - 1/3))/(pi*(625.0*sin(x)^2*(cos(pi*y) + 1)^2*(cos(2*x)/3 - 1/3)^2 + 1))