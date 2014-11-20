% Test expected distance
% Read all files an run exact expected distance (only for small instances)
clear all;
% -- Create output file
outputPath = '/home/andres/Documents/hga_vrpsd/outcomes/';
outputFile = 'exdist_outcome_20140722.csv';
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
path_wildchar = '/home/andres/Documents/hga_vrpsd/instances/*.dat';
path = '/home/andres/Documents/hga_vrpsd/instances/';
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
    %load instance  object
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
