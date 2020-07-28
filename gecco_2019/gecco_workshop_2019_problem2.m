function y = gecco_workshop_2019_problem1(x,S)

% y = gecco_workshop_2019_problem1(x,S)
%
% Second simple illustrative problem used in GECCO 2019 paper.
%
% Jonathan Fieldsend, University of Exeter, 2019
% See license information in package, available at 
% https://github.com/fieldsend/mo_lons


if (sum(x<=0)>0)
    error('input variable cannot be negative');
end
   
if (sum(x>=5)>0)
    error('input variable cannot be larger than 5');
end

x = ceil(x);
y = zeros(1,2);

A = [10 9 10 11 2;
     9 1 9 11 13;
     9 5 5 10 9;
     9 5 4 5 6;
     9 6 2 3 6]';
 
B = [10 9 10 11 3;
     9 9 9 12 15;
     9 3 2 7 8;
     9 6 2 3 7;
     9 5 4 5 6]';
 

y(1) = A(x(1),x(2));
y(2) = B(x(1),x(2));
 
 
    


