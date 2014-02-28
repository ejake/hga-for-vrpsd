% Read all files an run exact expected distance (only for small instances)
clear all;
% -- Create output file
outputPath = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Outcomes/';
outputFile = 'exdist_outcome_20131001.csv';
outputFullPath = [outputPath outputFile];
oFile = fopen(outputFullPath, 'a');
currentTime = clock;
fprintf(oFile,'Outcomes expected distance algorithms (%u-%u-%u, %u:%u):\n',currentTime(3), currentTime(2), currentTime(1), currentTime(4), currentTime(5));
fprintf(oFile,'instance;n;time exhaustive alg;time recursive alg; time dp alg; tour;Q;average distance;expected distance; dp expected distance\r\n');
% ---
% -- Specify location of input data (instances) 
%Windows:
%path_wildchar = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\*.dat';
%path = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\';
%Linux:
path_wildchar = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/*.dat';
path = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/';
listing = dir(path_wildchar);
for i=1:length(listing)    
    % -- Read instance file
    disp(listing(i).name);
    path_file = [path listing(i).name];
    fid=fopen(path_file, 'r');
    fn = str2num(fgetl(fid));
    data = [];
    for j=1:fn
        txt = fgetl(fid);
        data = [data ;str2num(txt)];
    end
    fclose(fid);
    % ---
    % -- Initialize instance values
    %Number of customers
    n = fn;
    %Demand Probability Distribution for each customer
    DD=[data(:,4) data(:,5)];
    %Location for each customers
    LL=[data(:,2) data(:,3)];
    %Factor for Q
    f=1;
    %load instance object
    instance = InstanceVrpsd(n, f, DD, LL);
    % ---
    % Asses expected distance
    pi = [];
    for k=1: instance.n
        pi = [pi Control(k,0)];
    end    
    tic;
    if instance.n == 5
        [frecD disD allDis] = distanceDistribution5n(pi, instance, 0);
    end
    if instance.n == 8
        [frecD disD allDis] = distanceDistribution8n(pi, instance, 1);
    end
    timeSpent = toc;    
    avgDist = mean(allDis);
    %avgDist = sum(frecD.*disD)/sum(frecD);
    % ---
    % Asses expected distance with recursive algorithm
    tic;
    expd = gammaBackwardEd(0,instance.Q,instance);
    timeSpent2 = toc;
    % ---
    % Asses dp algorithm to expected distance with memorization
    tour = 0:instance.n;
    tic;
    bed = backwardExpectedDistance(tour, instance);
    timeSpent3 = toc;
    % -- Write results in a output file
    fprintf(oFile, '%s;',listing(i).name);
    fprintf(oFile, ' %i',instance.n);
    fprintf(oFile, '%10.2f;',timeSpent);
    fprintf(oFile, '%10.2f;',timeSpent2);
    fprintf(oFile, '%10.2f;',timeSpent3);
    for k=1: (length(pi))
        fprintf(oFile, ' %i',pi(k).m);
    end
    fprintf(oFile, ';');
    fprintf(oFile, ' %i;',instance.Q);
    fprintf(oFile, ' %6.6f;',avgDist);
    fprintf(oFile, ' %6.6f;',expd);
    fprintf(oFile, ' %6.6f;',bed);
    fprintf(oFile,'\r\n');
    % ---
end
% -- Close outputFile
fclose(oFile);
% ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Read all files an run backward Expected Distance for each one
clear all;
% -- Create output file
%outputPath = '/home/undavid/Documents/MATLAB/VRPSD/outcomes/';
outputPath = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Outcomes/';
outputFile = 'bedcyclic_outcome.csv';
outputFullPath = [outputPath outputFile];
oFile = fopen(outputFullPath, 'a');
currentTime = clock;
fprintf(oFile,'Results RA (%u-%u-%u, %u:%u):\n',currentTime(3), currentTime(2), currentTime(1), currentTime(4), currentTime(5));
fprintf(oFile,'instance;n;time;tour;expected_distance_tour;avg_travel_distance_tour;var_travel_distance;avg_nodes_visited\r\n');
    
% ---
% -- Specify location of input data (instances) 
%Windows:
%path_wildchar = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\*.dat';
%path = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\';
%Linux:
path_wildchar = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/*.dat';
path = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/';
listing = dir(path_wildchar);
for i=1:length(listing)
    % -- Read instance file
    path_file = [path listing(i).name];
    fid=fopen(path_file, 'r');
    fn = str2num(fgetl(fid));
    data = [];
    for j=1:fn
        txt = fgetl(fid);
        data = [data ;str2num(txt)];
    end
    fclose(fid);
    % ---
    % -- Initialize instance values
    %Number of customers
    n = fn;
    %Demand Probability Distribution for each customer
    DD=[data(:,4) data(:,5)];
    %Location for each customers
    LL=[data(:,2) data(:,3)];
    %Factor for Q
    f=1;
    %load instance object
    instance = InstanceVrpsd(n, f, DD, LL);
    % ---
    % -- Backward expected distance
    tau = [0 1:n 0];
    tic;
    expectd2 = backwardExpectedDistance(instance, tau, 0, instance.Q);
    timeSpent = toc;
    % ---
    % -- Simulating tour distance
    pi = [];
    for k=1: (length(tau)-2)
        pi = [pi Control(tau(k+1),0)];
    end
    num_ite = 1000;
    rsim = zeros(1,num_ite);
    for j=1:num_ite
        [c trip] = simCyclicTripDistance(pi,instance,0);
        rsim(j) = c;
    end
    avg = mean(rsim);
    variance = var(rsim);
    stg_trip = length(trip);
    % ---
    % -- Write results in a output file
    fprintf(oFile, '%s;',listing(i).name);
    fprintf(oFile, ' %i',instance.n);
    fprintf(oFile, '%10.2f;',timeSpent);
    for k=2: (length(tau)-1)
        fprintf(oFile, ' %i',tau(k));
    end
    fprintf(oFile, ';');
    fprintf(oFile, ' %6.6f;',expectd2);
    fprintf(oFile, ' %6.6f;',avg);
    fprintf(oFile, ' %6.6f;',variance);
    fprintf(oFile, ' %i;',stg_trip);
    fprintf(oFile,'\r\n');
    % ---
end
% -- Close outputFile
fclose(oFile);
% ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Explore all possible TSP routes (exhaustive algorithm)
%It section should be used for small instances
clear all;
% -- Create output file
%outputPath = '/home/undavid/Documents/MATLAB/VRPSD/outcomes/';
outputPath = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Outcomes/';
outputFile = 'aped_outcome_20131002.csv';
outputFullPath = [outputPath outputFile];
oFile = fopen(outputFullPath, 'w');
currentTime = clock;
fprintf(oFile,'Results RA (%u-%u-%u, %u:%u):\n',currentTime(3), currentTime(2), currentTime(1), currentTime(4), currentTime(5));
fprintf(oFile,'instance;n;time;tour;expected_distance_tour\r\n');

% ---
% -- Specify location of input data (instances) 
%Windows:
%path_wildchar = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\*.dat';
%path = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\';
%Linux:
%path_wildchar = '/home/undavid/Documents/MATLAB/VRPSD/test_instance/*.dat';
%path = '/home/undavid/Documents/MATLAB/VRPSD/test_instance/';
path_wildchar = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/*.dat';
path = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/';
listing = dir(path_wildchar);
for i=1:length(listing)
    % -- Read instance file
    path_file = [path listing(i).name];
    fid=fopen(path_file, 'r');
    fn = str2num(fgetl(fid));
    data = [];
    for j=1:fn
        txt = fgetl(fid);
        data = [data ;str2num(txt)];
    end
    fclose(fid);
    % ---
    % -- Initialize instance values
    %Number of customers
    n = fn;
    %Demand Probability Distribution for each customer
    DD=[data(:,4) data(:,5)];
    %Location for each customers
    LL=[data(:,2) data(:,3)];
    %Factor for Q
    f=1;
    %load instance object
    instance = InstanceVrpsd(n, f, DD, LL);
    disp(listing(i).name);
    % ---
    %Exhaustive exploration of expected distance for all tau permutations
    sub_tau = 1:n;
    all_perms_tau = perms(sub_tau);
    disp(size(all_perms_tau,1));
    for j=1:size(all_perms_tau,1)
        tau = [0 all_perms_tau(j,:) 0];
        tic;
        expectd2 = backwardExpectedDistance(tau, instance);
        timeSpent = toc;
        % ---
    % -- Write results in a output file    
    fprintf(oFile, '%s;',listing(i).name);
    fprintf(oFile, '%i;',n);
    fprintf(oFile, '%10.2f;',timeSpent);
    for k=1: length(tau)
        fprintf(oFile, ' %i',tau(k));
    end
    fprintf(oFile, ';');
    fprintf(oFile, ' %10.2f;',timeSpent);
    fprintf(oFile, ' %6.6f;',expectd2);
    fprintf(oFile,'\r\n');
    % ---
    end    
    % ---
end
% -- Close outputFile
fclose(oFile);
% ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read all files an run RA for each one
clear all;
% -- Create output file
%outputPath = '/home/undavid/Documents/MATLAB/VRPSD/outcomes/';
outputPath = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Outcomes/';
outputFile = 'ra_outcome20140225.csv';
outputFullPath = [outputPath outputFile];
oFile = fopen(outputFullPath, 'a');
currentTime = clock;
fprintf(oFile,'Results RA (%u-%u-%u, %u:%u):\n',currentTime(3), currentTime(2), currentTime(1), currentTime(4), currentTime(5));
fprintf(oFile,'instance;n;time;policy_m;policy_a;time_expected_distance;expected_distance\r\n');
%fprintf(oFile, '%10.2f; %i; %i; %10.2f; %6.6f; %6.6f; %i; %6.6f; %i \r\n',timeSpent,policy_m,policy_a,timeExSpent,expectd,simCost0, simTour0, simCost1, simTour1);
    
% ---
% -- Specify location of input data (instances) 
%Windows:
%path_wildchar = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\*.dat';
%path = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\';
%Linux:
path_wildchar = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/*.dat';
path = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/';
listing = dir(path_wildchar);
for i=1:length(listing)
    % -- Read instance file
    path_file = [path listing(i).name];
    fid=fopen(path_file, 'r');
    fn = str2num(fgetl(fid));
    data = [];
    for j=1:fn
        txt = fgetl(fid);
        data = [data ;str2num(txt)];
    end
    fclose(fid);
    disp(listing(i).name);
    % ---
    % -- Initialize instance values
    %Number of customers
    n = fn;
    %Demand Probability Distribution for each customer
    DD=[data(:,4) data(:,5)];
    %Location for each customers
    LL=[data(:,2) data(:,3)];
    %Factor for Q
    f=1;
    %load instance object
    instance = InstanceVrpsd(n, f, DD, LL);
    %Initial state
    st = zeros(1, instance.n+2); %state = (l, q_l, j_1, ... , j_n)
    st(1)=0; st(2) = instance.Q;
    st(3:instance.n+2) = -99;
    % ---
    % -- Rollout algorithm
    % base tour (a priori solution) - Cyclic heuristic
    l = 1;
    tau_l = [l:instance.n 1:l-1];    
    tic; % start timer
    pi = rollout (instance, State(instance.n, instance.Q), tau_l);
    timeSpent = toc; % stop timer
    % ---
    % -- Prepare outcomes data
    policy_m = zeros(1,length(pi));
    policy_a = zeros(1,length(pi));
    for k=1: length(pi)
        policy_m(k) = pi(k).m;
        policy_a(k) = pi(k).a;
    end
    %Expected distance    
    tic % start timer
    expectd = -1;
    if length(pi) == n
        expectd = backwardExpectedDistance([0 policy_m], instance);
    end
    timeExSpent =toc; % stop timer
    
    % ---
    % -- Write results in a output file    
    fprintf(oFile, '%s;',listing(i).name);%instance name
    fprintf(oFile, '%i;',instance.n);%instance size
    fprintf(oFile, '%10.2f;',timeSpent);%running time of RA
    fprintf(oFile, ' %s;', num2str(policy_m));%tour constructed by RA
    fprintf(oFile, ' %s;', num2str(policy_a));% proactive replanishment of policy contructed by RA
    fprintf(oFile, ' %10.2f;',timeExSpent);%running time of tour backward expected distance contructed by RA
    fprintf(oFile, ' %6.6f;',expectd);%expected distance    
    fprintf(oFile,'\r\n');
    % ---
    clear instance;
end
% -- Close outputFile
fclose(oFile);
% ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read all files an run GA for each one
clear all;
% -- Create output file
%outputPath = '/home/undavid/Documents/MATLAB/VRPSD/outcomes/';
outputPath = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Outcomes/';
outputFile = 'ra_outcome.csv';
outputFullPath = [outputPath outputFile];
oFile = fopen(outputFullPath, 'w');
currentTime = clock;
fprintf(oFile,'Results GA (%u-%u-%u, %u:%u):\n',currentTime(3), currentTime(2), currentTime(1), currentTime(4), currentTime(5));
fprintf(oFile,'instance;n;time;tour;expected_distance\r\n');
%fprintf(oFile, '%10.2f; %i; %i; %10.2f; %6.6f; %6.6f; %i; %6.6f; %i \r\n',timeSpent,policy_m,policy_a,timeExSpent,expectd,simCost0, simTour0, simCost1, simTour1);
    
% ---
% -- Specify location of input data (instances) 
%Windows:
%path_wildchar = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\*.dat';
%path = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\';
%Linux:
path_wildchar = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/*.dat';
path = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/';
listing = dir(path_wildchar);
for i=1:length(listing)
    % -- Read instance file
    path_file = [path listing(i).name];
    fid=fopen(path_file, 'r');
    fn = str2num(fgetl(fid));
    data = [];
    for j=1:fn
        txt = fgetl(fid);
        data = [data ;str2num(txt)];
    end
    fclose(fid);
    % ---
    % -- Initialize instance values
    %Number of customers
    n = fn;
    %Demand Probability Distribution for each customer
    DD=[data(:,4) data(:,5)];
    %Location for each customers
    LL=[data(:,2) data(:,3)];
    %Factor for Q
    f=1;
    %load instance object
    instance = InstanceVrpsd(n, f, DD, LL);
    % ---
    % -- Genetic algorithm
    pop_size = 
    tic; % start timer
    [optRoute, expectd] = vrpsd_ga(instance, 10, 10, 1, 1);
    timeSpent = toc; % stop timer
    % ---
    % -- Prepare outcomes data
    tour = [0 optRoute 0];
    
    % ---
    % -- Write results in a output file    
    fprintf(oFile, '%s;',listing(i).name);
    fprintf(oFile, ' %i;',instance.n);
    fprintf(oFile, '%10.2f;',timeSpent);
    for k=1: length(optRoute)
        fprintf(oFile, ' %i',tour(k));
    end
    fprintf(oFile, ';');
    fprintf(oFile, ' %6.6f',expectd);
    fprintf(oFile,'\r\n');
    % ---
end
% -- Close outputFile
fclose(oFile);
% ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read only one instance
clear all;
%Windows:
%fid=fopen('D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\i_5r1.dat', 'rt');
%Linux:
fid=fopen('/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/dummy_n5.dat', 'rt');
%fid=fopen('/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/i_5r1.dat', 'rt');
%fid=fopen('/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/i_8r1.dat', 'rt');
%fid=fopen('/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/i_20r2.dat', 'rt');
fn = str2num(fgetl(fid));
data = [];
for i=1:fn
    txt = fgetl(fid);
    data = [data ;str2num(txt)];
end
fclose(fid);

%% Initialize instance values
%Number of customers
n = fn;
%Demand Probability Distribution for each customer
DD=[data(:,4) data(:,5)];
%Location for each customers
LL=[data(:,2) data(:,3)];
%Factor for Q
f=1;
%load instance object
instance = InstanceVrpsd(n, f, DD, LL);
%Initial state
st = zeros(1, instance.n+2); %state = (l, q_l, j_1, ... , j_n)
st(1)=0; st(2) = instance.Q;
st(3:instance.n+2) = -99;

%% Run expected distance algorithm
tau = [0 1:n 0];
%tau = [0 5 3 2 1 4 0];
tic;
expectd = expectedDistance(instance, tau, 0, instance.Q);
disp('Forward expected distance algorithm:');
disp(expectd);
timeSpent = toc;
disp('Time: ');
disp(timeSpent);

%% Run backward expected distance
tau = [0 1:n 0];
%tau = [1:n 0];
%tau = [0 4 5 2 3 1 0];
tic;
expectd2 = backwardExpectedDistance(instance, tau, 0, instance.Q);
timeSpent = toc;
%% Run rollout algorithm
tau = [0 1:n 0];
%tau = [0 1 2 3 4 5 0];
x0 = State(instance.n, instance.Q);
tic;
rapolicy = rollout( tau, instance, x0 );
timeSpent = toc;
%% Run GA
% -- Genetic algorithm
iters = 100;
fp = 2;
pop_size = instance.n*fp;
tic; % start timer
[optRoute, expectd] = vrpsd_ga(instance, pop_size, iters, 0, 1);
timeSpent = toc; % stop timer

%% Run hybrid GA

%% Simulation
% Demand distribution
a = avgDemand();
x = 7:17;
plot(a,x);
p = a./243;
plot(x,p);

%% -- Simulating tour distance

num_ite = 1000;
rsim = zeros(1, num_ite);
instance.Q = 7;

pi = [];
for k=1: instance.n
    pi = [pi Control(k,0)];
end

for j=1:num_ite
    [c trip] = simCyclicTripDistance(pi,instance,0);
    rsim(j) = c;
end
avg = mean(rsim);
variance = var(rsim);
stg_trip = length(trip);

%% Real expected distance

pi = [];
for k=1: instance.n
    pi = [pi Control(k,0)];
end

[frecD disD allDis] = distanceDistribution5n(pi, instance, 0);
avgDist = sum(frecD.*disD)/sum(frecD);
plot(frecD);
hist(allDis);
fprintf('Real expected distance: %6.4f\n', mean(allDis));
%disp(avgDist);

%% Review real backward expected distance
J= zeros(instance.n + 1, instance.Q +1);
J(6,:) = instance.d(5+1,1);%d(n,0) q_l is irrelevant

%Demand of customer 5
J(5,5:7) = instance.d(4+1,4+2);
J(5,4) = (instance.d(4+1,4+2)+2*instance.d(1,4+2))*1/3 + instance.d(4+1,4+2)*(2/3);
J(5,3) = (instance.d(4+1,4+2)+2*instance.d(1,4+2))*2/3 + instance.d(4+1,4+2)*(1/3);
J(5,2) = (instance.d(4+1,4+2)+2*instance.d(1,4+2))*3/3 + instance.d(4+1,4+2)*(0);
J(5,1) = instance.d(4+1,1)+instance.d(1,4+2);

%Demand of customer 4
J(4,4:7) = instance.d(3+1,3+2);
J(4,3) = (instance.d(3+1,3+2)+2*instance.d(1,3+2))*1/3 + instance.d(3+1,3+2)*(2/3);
J(4,2) = (instance.d(3+1,3+2)+2*instance.d(1,3+2))*2/3 + instance.d(3+1,3+2)*(1/3);
J(4,1) = instance.d(3+1,1)+instance.d(1,3+2);

%Demand of customer 3
J(3,4:7) = instance.d(2+1,2+2);
J(3,3) = (instance.d(2+1,2+2)+2*instance.d(1,2+2))*1/3 + instance.d(2+1,2+2)*(2/3);
J(3,2) = (instance.d(2+1,2+2)+2*instance.d(1,2+2))*2/3 + instance.d(2+1,2+2)*(1/3);
J(3,1) = instance.d(2+1,1)+instance.d(1,2+2);

%Demand of customer 2
J(2,4:7) = instance.d(1+1,1+2);
J(2,3) = (instance.d(1+1,1+2)+2*instance.d(1,1+2))*1/3 + instance.d(1+1,1+2)*(2/3);
J(2,2) = (instance.d(1+1,1+2)+2*instance.d(1,1+2))*2/3 + instance.d(1+1,1+2)*(1/3);
J(2,1) = instance.d(1+1,1)+instance.d(1,1+2);

J(1,:) = instance.d(1+1,1);%d(1,0) q_l is irrelevant

%% Review real backward expected distance (2)
% (1/10/2013)

%q_l frecuency to asses probability for dummy instance
PQ = [0 0 1/3 1/3 1/3 0 0 ; 
    2/9 4/9 2/9 1/9 0 0 0 ;
    1/4 1/4 1/12 1/12 1/6 1/6 0 ;
    1/6 1/6 1/6 1/6 1/6 1/6 0 ;
    1/6 2/9 2/9 1/6 1/9 1/9 0];

tic; % start timer
e = gammaEd(instance.Q,0,instance,PQ);
timeSpent = toc; % stop timer
disp('Forward expected algorithm (2): ');
disp(e);
disp('Time: ');
disp(timeSpent);


%% New algorithm to asses expected distance
%------------------
% - In forward approach is to close to real expected distance, applied
% to dummy instance was 0.1 different 
% - This algorithm implementation assume a tour sequential, i.e. \tau = [0 1 2 ... n]
% (1/10/2013)
%------------------
l=0;
ql=instance.Q;
tic;
el = gammaBackwardEd(l,ql,instance);
timeSpent = toc;
fprintf('Expected distance to sequential tour, starting in %d with q_%d = %d: %6.4f (%6.4f sec)\n', l, l, ql, el, timeSpent);
%Expected distance to sequential tour, starting in 0 with q_0 = 6: 6.4944 (0.0672 sec)


%% Assesing expected distance
%------------------
% Implementing dynamic programming backward algorithm with memorization
% (2/10/2013)
%------------------

tic;
tour = 0:instance.n;
bed = backwardExpectedDistance(tour, instance);
timeSpent = toc;
fprintf('Backward Expected distance to sequential tour: %6.4f (%6.4f sec)\n', bed, timeSpent);

%% Computing rollout algorithm
%--------------------
% Rollout algorithm implementation using new (right) algorithm to asses
% expected distance
% (12/2/2014)
%--------------------


% base tour (a priori solution) - Cyclic heuristic
l = 1;
tau_l = [l:instance.n 1:l-1];
%tau_l = [2 3 5 1 4];
tau = zeros(1,instance.n);

pi = rollout (instance, State(instance.n, instance.Q), tau_l);

fprintf('pi = (');
for i=1:length(pi)
    fprintf(' (%i,%i)',pi(i).m,pi(i).a);
    tau(i) = pi(i).m;
end
fprintf(')\n');

tic;
bed = backwardExpectedDistance([0 tau], instance);
timeSpent = toc;
fprintf('Backward Expected distance pi tour: %6.4f (%6.4f sec)\n', bed, timeSpent);
