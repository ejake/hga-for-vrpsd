%% Read all files an run RA for each one
clear all;
% -- Create output file
%outputPath = '/home/undavid/Documents/MATLAB/VRPSD/outcomes/';
outputPath = '/home/andres/Documents/hga_vrpsd/outcomes/';
outputFile = 'ra_outcome20140722.csv';
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
path_wildchar = '/home/andres/Documents/hga_vrpsd/instances/large/*.dat';
path = '/home/andres/Documents/hga_vrpsd/instances/large/';
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
    clear n;
end
% -- Close outputFile
fclose(oFile);
% ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
