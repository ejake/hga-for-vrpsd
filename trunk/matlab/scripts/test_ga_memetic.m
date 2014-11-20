%% Read all files an run GA for each one
clear all;
% -- Create output file
%outputPath = '/home/undavid/Documents/MATLAB/VRPSD/outcomes/';
outputPath = '/home/andres/Documents/hga_vrpsd/outcomes/';
outputFile = 'ga_memetic_20141119.csv';
outputFullPath = [outputPath outputFile];
oFile = fopen(outputFullPath, 'a');
currentTime = clock;
fprintf(oFile,'Results GA - memetic (%u-%u-%u, %u:%u):\n',currentTime(3), currentTime(2), currentTime(1), currentTime(4), currentTime(5));
fprintf(oFile,'instance;n;time;policy_m;policy_a;expected_distance;pop_size;num_iter;epsilon;m_iter_without_change;prob_mutation;alpha\r\n');
%fprintf(oFile, '%10.2f; %i; %i; %10.2f; %6.6f; %6.6f; %i; %6.6f; %i \r\n',timeSpent,policy_m,policy_a,timeExSpent,expectd,simCost0, simTour0, simCost1, simTour1);
    
% ---
% -- Specify location of input data (instances) 
%Windows:
%path_wildchar = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\*.dat';
%path = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\';
%Linux:
path_wildchar = '/home/andres/Documents/hga_vrpsd/instances/small/*.dat';
path = '/home/andres/Documents/hga_vrpsd/instances/small/';
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
    % -- Genetic algorithm    
    tic; % start timer
    
    pop_size = instance.n;
    num_iter = 60;
    epsilon = 0.01;
    m_without_change = 15;
    prob_mutation = 0.1;
    alpha = 0.5;
    [pi_ga ed_ga fig] = vrpsd_ga(instance, pop_size, num_iter, epsilon, m_without_change, prob_mutation, alpha, 1, 0, 1);
    
    timeSpent = toc; % stop timer
    % ---
    % -- Prepare outcomes data
    policy_m = zeros(1,length(pi));
    policy_a = zeros(1,length(pi));
    for k=1: length(pi_ga)
        policy_m(k) = pi_ga(k).m;
        policy_a(k) = pi_ga(k).a;
    end
    
    % ---
    % -- Write results in a output file    
    fprintf(oFile, '%s;',listing(i).name);
    fprintf(oFile, ' %i;',instance.n);
    fprintf(oFile, '%10.2f;',timeSpent);
    fprintf(oFile, ' %s;', num2str(policy_m));%tour constructed by RA
    fprintf(oFile, ' %s;', num2str(policy_a));% proactive replanishment of policy contructed by RA    
    fprintf(oFile, ' %6.6f',ed_ga);
    fprintf(oFile, ' %i',pop_size);
    fprintf(oFile, ' %i',num_iter);
    fprintf(oFile, ' %6.6f',epsilon);
    fprintf(oFile, ' %i',m_without_change);
    fprintf(oFile, ' %6.6f',prob_mutation);
    fprintf(oFile, ' %6.6f',alpha);
    fprintf(oFile,'\r\n');
    % ---
    hgexport(fig,['/home/andres/Documents/hga_vrpsd/outcomes/ga/memetic/' listing(i).name '.eps']);
    close(fig);
end
% -- Close outputFile
fclose(oFile);
% ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
