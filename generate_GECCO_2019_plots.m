%generate_GECCO_2019_plots
%
% Script regenerates all the network plots from the GECCO 2019 paper
%
% Jonathan E. Fieldsend and Khulood Alyahya. 2019. 
% Visualising the Landscape of Multi-Objective Problems using Local Optima 
% Networks. In Genetic and Evolutionary Computation Conference Companion 
% (GECCO ?19 Companion), July 13?17, 2019, Prague, Czech Republic. ACM, 
% New York, NY, USA.
%
%
% Jonathan Fieldsend, University of Exeter, 2019
% See license information in package, available at 
% https://github.com/fieldsend/mo_lons

fprintf('Generating PLOS-net, DNON and PLON plots for the first small bi-objective problem\n');

[X,Y,state,neighbours,w,YY] = exaustive_generate_lon('gecco_workshop_2019_problem1', [], false, [4.5 4.5], [0.5 0.5], 5, 2);
figure; DNO=plot_PLOS_net(X,Y,state,neighbours,5,10);

[A,B,EE,C,PO] = process_d_lon(X,Y,w,YY,state, neighbours);
figure; generate_lon(B, C,EE, [],10,50);

[X,Y,state,neighbours] = exaustive_generate_lon('gecco_workshop_2019_problem1', [], true, [4.5 4.5], [0.5 0.5], 5, 2);
[V,B,Adj,EE,C] = process_p_lon(X,Y,state,neighbours);
figure; generate_lon(B, C,EE, [],10,50,Adj);


fprintf('Generating PLOS-net, DNON and PLON plots for the second small bi-objective problem\n');

[X,Y,state,neighbours,w,YY] = exaustive_generate_lon('gecco_workshop_2019_problem2', [], false, [4.5 4.5], [0.5 0.5], 5, 2);
figure; DNO=plot_PLOS_net(X,Y,state,neighbours,5,10);

[A,B,EE,C,PO] = process_d_lon(X,Y,w,YY,state, neighbours);
figure; generate_lon(B, C,EE, [],10,50);

[X,Y,state,neighbours] = exaustive_generate_lon('gecco_workshop_2019_problem2', [], true, [4.5 4.5], [0.5 0.5], 5, 2);
[V,B,Adj,EE,C] = process_p_lon(X,Y,state,neighbours);
figure; generate_lon(B, C,EE, [],10,50,Adj);

fprintf('Generating PLOS-net, DNON and PLON plots for the first five-objective problem\n');
load gecco_2019_mo_lon_many_problem_parameters1

[X,Y,state,neighbours,w,YY] = exaustive_generate_lon('gecco_workshop_2019_problem3', Meta, false, [39.5 39.5], [0.5 0.5], 40, 5);
figure; DNO=plot_PLOS_net(X,Y,state,neighbours,1.5,3);
 
[A,B,EE,C,PO] = process_d_lon(X,Y,w,YY,state, neighbours);
figure; generate_lon(B, C,EE, [],10,50);
 
[X,Y,state,neighbours] = exaustive_generate_lon('gecco_workshop_2019_problem3', Meta, true, [39.5 39.5], [0.5 0.5], 40, 5);
[V,B,Adj,EE,C] = process_p_lon(X,Y,state,neighbours);
figure; [s, t, weights] = generate_lon(B, C,EE, [],2,100,Adj);

fprintf('Generating PLOS-net, DNON and PLON plots for the second five-objective problem\n');
load gecco_2019_mo_lon_many_problem_parameters2

[X,Y,state,neighbours,w,YY] = exaustive_generate_lon('gecco_workshop_2019_problem3', Meta, false, [39.5 39.5], [0.5 0.5], 40, 5);
figure; DNO=plot_PLOS_net(X,Y,state,neighbours,1.5,3);
 
[A,B,EE,C,PO] = process_d_lon(X,Y,w,YY,state, neighbours);
figure; generate_lon(B, C,EE, [],10,50);
 
[X,Y,state,neighbours] = exaustive_generate_lon('gecco_workshop_2019_problem3', Meta, true, [39.5 39.5], [0.5 0.5], 40, 5);
[V,B,Adj,EE,C] = process_p_lon(X,Y,state,neighbours);
figure; [s, t, weights] = generate_lon(B, C,EE, [],2,100,Adj);

