function y = gecco_workshop_2019_problem1(x,S)

% y = gecco_workshop_2019_problem1(x,S)
%
% First simple illustrative problem used in GECCO 2019 paper.
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

A = [9 9 8 8 4;
     9 8 9 8 5;
     9 7 8 7 3;
     8 7 6 4 3;
     7 5 6 5 7]';
 
B = [9 9 8 8 4;
     9 7 8 7 4;
     9 7 8 7 7;
     8 6 7 5 6;
     7 7 3 4 2]';
 

y(1) = A(x(1),x(2));
y(2) = B(x(1),x(2));
 
 
    


