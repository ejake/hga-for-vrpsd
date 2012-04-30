%% Read all files an run backward Expected Distance for each one
clear all;
% -- Create output file
%outputPath = '/home/undavid/Documents/MATLAB/VRPSD/outcomes/';
outputPath = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Outcomes/';
outputFile = 'bed_outcome.csv';
outputFullPath = [outputPath outputFile];
oFile = fopen(outputFullPath, 'a');
currentTime = clock;
fprintf(oFile,'Results RA (%u-%u-%u, %u:%u):\n',currentTime(3), currentTime(2), currentTime(1), currentTime(4), currentTime(5));
fprintf(oFile,'instance;n;time;tour;expected_distance_tour;avg_travel_distance_tour;avg_nodes_visited\r\n');
    
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
    for k=1: (length(tau)-1)
        pi = [pi Control(tau(k+1),0)];
    end
    num_ite = 1000;
    rsim = zeros(1,num_ite);
    for j=1:nume_ite
        [c trip] = simTripDistance(pi,instance,0);
        rsim(j) = c;
    end
    avg = mean(rsim);
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
    fprintf(oFile, ' %i;',stg_trip);
    fprintf(oFile,'\r\n');
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
        policy_m = rapolicy(i).m;
        policy_a = rapolicy(i).a;
    end
    %Expected distance
    eTau = [0 policy_m 0];
    tic % start timer
    expectd = -1;
    if length(eTau) == n+2
        expectd = expectedDistance(instance, eTau, 0, instance.Q);
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

%% Read only one instance
clear all;
%Windows:
%fid=fopen('D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\i_5r1.dat', 'rt');
%Linux:
fid=fopen('/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/i_5r1.dat', 'rt');
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

