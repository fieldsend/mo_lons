function y = gecco_workshop_2019_problem3(x,S)

if (sum(x<=0)>0)
    error('input variable cannot be negative');
end
   
if (sum(x>=100)>0)
    error('input variable cannot be larger than 5');
end

x = ceil(x);
y = zeros(1,5);

y(1) = S.F1(x(1),x(2));
y(2) = S.F2(x(1),x(2));
y(3) = S.F3(x(1),x(2));
y(4) = S.F4(x(1),x(2));
y(5) = S.F5(x(1),x(2)); 
 
    


