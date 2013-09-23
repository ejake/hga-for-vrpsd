% Read all files an run exact expected distance (only for small instances)
clear all;
% -- Create output file
%outputPath = '/home/undavid/Documents/MATLAB/VRPSD/outcomes/';
%outputPath = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Outcomes/';
outputPath = '/home/lisi/Documents/vrpsd/outcomes/';
outputFile = 'exdist_cyclic_outcome.csv';
outputFullPath = [outputPath outputFile];
oFile = fopen(outputFullPath, 'a');
currentTime = clock;
fprintf(oFile,'Results RA (%u-%u-%u, %u:%u):\n',currentTime(3), currentTime(2), currentTime(1), currentTime(4), currentTime(5));
fprintf(oFile,'instance;n;time;tour;Q;average distance;total demand;frecuency;distance\r\n');
% ---
% -- Specify location of input data (instances) 
%Windows:
%path_wildchar = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\*.dat';
%path = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\';
%Linux:
path_wildchar = '/home/lisi/Documents/vrpsd/instances/*.dat';
path = '/home/lisi/Documents/vrpsd/instances/';
%path_wildchar = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/*.dat';
%path = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/';
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
        [frecD disD] = distanceDistribution5n(pi, instance, 1);
    end
    if instance.n == 8
        [frecD disD] = distanceDistribution8n(pi, instance, 1);
    end
    timeSpent = toc;
    avgDist = sum(frecD.*disD)/sum(frecD);
    % -- Write results in a output file
    fprintf(oFile, '%s;',listing(i).name);
    fprintf(oFile, ' %i',instance.n);
    fprintf(oFile, '%10.2f;',timeSpent);
    for k=1: (length(pi))
        fprintf(oFile, ' %i',pi(k).m);
    end
    fprintf(oFile, ';');
    fprintf(oFile, ' %i;',instance.Q);
    fprintf(oFile, ' %6.6f;',avgDist);
    for k=1: (length(frecD))
        fprintf(oFile, ' %i',frecD(k));
    end
    fprintf(oFile, ';');
    for k=1: (length(disD))
        fprintf(oFile, ' %i',disD(k));
    end
    fprintf(oFile, ';');
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
outputFile = 'aped_outcome.csv';
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
        expectd2 = backwardExpectedDistance(instance, tau, 0, instance.Q);
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
outputFile = 'ra_outcome.csv';
outputFullPath = [outputPath outputFile];
oFile = fopen(outputFullPath, 'w');
currentTime = clock;
fprintf(oFile,'Results RA (%u-%u-%u, %u:%u):\n',currentTime(3), currentTime(2), currentTime(1), currentTime(4), currentTime(5));
fprintf(oFile,'instance;n;time;policy_m;policy_a;time_expected_distance;expected_distance;avg_travel_distance_tour;tour_travel_distance_sim_tou;avg_travel_distance_simulate_policy;tour_travel_distance_sim_policy\r\n');
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
    %Initial state
    st = zeros(1, instance.n+2); %state = (l, q_l, j_1, ... , j_n)
    st(1)=0; st(2) = instance.Q;
    st(3:instance.n+2) = -99;
    % ---
    % -- Rollout algorithm
    tau = [0 1:n 0];
    x0 = State(instance.n, instance.Q);
    tic; % start timer
    rapolicy = rollout( tau, instance, x0 );
    timeSpent = toc; % stop timer
    % ---
    % -- Prepare outcomes data
    policy_m = zeros(1,length(rapolicy));
    policy_a = zeros(1,length(rapolicy));
    for k=1: length(rapolicy)
        policy_m(k) = rapolicy(k).m;
        policy_a(k) = rapolicy(k).a;
    end
    %Expected distance
    eTau = [0 policy_m 0];
    tic % start timer
    expectd = -1;
    if length(eTau) == n+2
        expectd = backwardExpectedDistance(instance, eTau, 0, instance.Q);
    end
    timeExSpent =toc; % stop timer
    %Simulating tour
    simCost0 = -1;
    simTour0 = [];
    simCost1 = -1;
    simTour1 = [];
    if length(eTau) == n+2
        [simCost0 simTour0] = simTripDistance(rapolicy, instance, 0);
        [simCost1 simTour1] = simTripDistance(rapolicy, instance, 1);
    end 
    
    % ---
    % -- Write results in a output file    
    fprintf(oFile, '%s;',listing(i).name);
    fprintf(oFile, '%10.2f;',timeSpent);
    for k=1: length(policy_m)
        fprintf(oFile, ' %i',policy_m(k));
    end
    fprintf(oFile, ';');
    for k=1: length(policy_a)
        fprintf(oFile, ' %i',policy_a(k));
    end
    fprintf(oFile, ';');
    fprintf(oFile, ' %10.2f;',timeExSpent);
    fprintf(oFile, ' %6.6f;',expectd);
    fprintf(oFile, ' %6.6f;',simCost0); 
    for k=1: length(simTour0)
        fprintf(oFile, ' %i',simTour0(k));
    end
    fprintf(oFile, ';');
    fprintf(oFile, ' %6.6f;',simCost1); 
    for k=1: length(simTour1)
        fprintf(oFile, ' %i',simTour1(k));
    end
    fprintf(oFile,'\r\n');
    % ---
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
timeSpent = toc;

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
disp('Real expected distance');
disp(avgDist);

%% Review real backward expected distance
J= zeros(instance.n + 1, instance.Q +1);
J(6,:) = instance.d(5+1,1);%d(n,0) q_l is irrelevant

%Demand of customer 5
J(5,5:7) = instance.d(4+1,4+2)
J(5,4) = (instance.d(4+1,4+2)+2*instance.d(1,4+2))*1/3 + instance.d(4+1,4+2)*(2/3)
J(5,3) = (instance.d(4+1,4+2)+2*instance.d(1,4+2))*2/3 + instance.d(4+1,4+2)*(1/3)
J(5,2) = (instance.d(4+1,4+2)+2*instance.d(1,4+2))*3/3 + instance.d(4+1,4+2)*(0)
J(5,1) = instance.d(4+1,1)+instance.d(1,4+2)

%Demand of customer 4
J(4,4:7) = instance.d(3+1,3+2)
J(4,3) = (instance.d(3+1,3+2)+2*instance.d(1,3+2))*1/3 + instance.d(3+1,3+2)*(2/3)
J(4,2) = (instance.d(3+1,3+2)+2*instance.d(1,3+2))*2/3 + instance.d(3+1,3+2)*(1/3)
J(4,1) = instance.d(3+1,1)+instance.d(1,3+2)

%Demand of customer 3
J(3,4:7) = instance.d(2+1,2+2)
J(3,3) = (instance.d(2+1,2+2)+2*instance.d(1,2+2))*1/3 + instance.d(2+1,2+2)*(2/3)
J(3,2) = (instance.d(2+1,2+2)+2*instance.d(1,2+2))*2/3 + instance.d(2+1,2+2)*(1/3)
J(3,1) = instance.d(2+1,1)+instance.d(1,2+2)

%Demand of customer 2
J(2,4:7) = instance.d(1+1,1+2)
J(2,3) = (instance.d(1+1,1+2)+2*instance.d(1,1+2))*1/3 + instance.d(1+1,1+2)*(2/3)
J(2,2) = (instance.d(1+1,1+2)+2*instance.d(1,1+2))*2/3 + instance.d(1+1,1+2)*(1/3)
J(2,1) = instance.d(1+1,1)+instance.d(1,1+2)

J(1,:) = instance.d(1+1,1);%d(1,0) q_l is irrelevant




