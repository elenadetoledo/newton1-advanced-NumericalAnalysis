%%%PROJECT 1 - PART 1 - ELENA DE TOLEDO HERNÁNDEZ
%SECTION 3.2 - PROBLEMS 3, 14(b)

%FIRST PART - SECTION 3.2 PROBLEM 3
%Find the positive minimum point of the function f(x)=(x)^-2 * tanx by
%computing the zeros of f' using Newton's method

clear 
clc

f = @(x) x^-2 *tan(x);
% df1 = f' computed by hand to rearrange the function and avoid difference of
%two close numbers
df1 = @(x) -(2*tan(x)-x*(sec(x))^2)/(x^3);
%Similarly, df2 = f''
df2 = @(x) ((2*x^2*(sec(x))^2+6)*tan(x)-4*x*(sec(x))^2)/x^4;

%fprintf('\nInitial guess \tNumber of its \tValue xn1 \tValue f(xn1)');
%fprintf('\nxn \tNumber of its \tValue of x \tImage of f(x)');
minf = 0;
minx = 0;
for j= -1000: pi : 1000
    xn = j;
    count = 0;
    while (count < 50) && (df2(xn) ~= 0)

        
        count = count + 1;
        xn1 = xn - (df1(xn)/df2(xn));
        %fprintf('\n \t%f \t%d \t%f \t%f',xn,count,xn1,f(xn1));%Each zero in each range

        xn = xn1;
        

        if abs(df1(xn) - 0) < 10^-6
            count = 1000;
        end


    end
    
    if (j == -1000) || (f(xn) < minf)
        minx = xn;
        minf = f(xn);
    end
    %fprintf('\n%f \t%d \t%f \t%f',xn,count,minx,minf);
    
        
 
end
fprintf('\nSEC 3.2 PROBLEM 3: The minimum of f(x) is %f,which is obtained at x = %f',minf,minx)

%SECOND PART - SECTION 3.2 PROBLEM 14(b)
%Using Newton's method, find the roots of the nonlinear systems. One
%solution will suffice
%x + e^(-1x)+y^3 = 0, x^2 + 2xy - y^(2) + tan(x) = 0

f1 = @(x,y) x+exp(1)^(-x) + y^(3);
f2 = @(x,y) x^(2) +2*x*y -y^(2) +tan(x);
df1x = @(x,y) 1 - exp(1)^(-x);
df1y = @(x,y) 3*y^(2);
df2x = @(x,y) 2*x + 2*y + (sec(x))^(2);
df2y = @(x,y) 2*x - 2*y;

%Initial guesses
x0 = 1.5;
y0 = -2;
pvect = [x0 ; y0];
N = 0;
tolerance = 10^(-6);
error = 1000; %Initialize to a random first value to enter the while loop
fprintf('\nInitial guess \tNumber of its \tValue xn+1 \tValue yn+1 \tError');

while (error > tolerance) && (N < 50)
    JacobianMatrix = [df1x(x0,y0) df1y(x0,y0); df2x(x0,y0) df2y(x0,y0)];
    N = N + 1;
    vfunct = [f1(x0,y0);f2(x0,y0)];
    %Newton function
    h = -inv(JacobianMatrix)*vfunct;
    next = pvect + h;
    error = norm(pvect -next);
    fprintf('\nx=%f,y=%f \t%d \t%f \t%f \t%f',x0,y0,N,next(1),next(2),error);
    %Update initial values
    x0 = next(1);
    y0 = next(2);
    pvect = [x0 ; y0];
    
end

fprintf('\nSEC 3.2 PROBLEM 14(b): The solution is x=%f and y=%f ',pvect(1),pvect(2));

