%% Read all files an run everyone
%Windows:
%path_wildchar = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\*.dat';
%path = 'D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\';
%Linux:
path_wildchar = '/media/DATA_/Documents/Seminario\de\Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/*.dat';
path = '/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/';
listing = dir(path_wildchar);
for i=1:length(listing)
    path_file = [path listing(i).name];
    fid=fopen(path_file, 'r');
    fn = str2num(fgetl(fid));
    data = [];
    for j=1:fn
        txt = fgetl(fid);
        data = [data ;str2num(txt)];
    end
    fclose(fid);
    %Initialize instance values
end
%Initialize instance values
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

%Run algorithm
tau = [0 5 3 2 1 4 0];
tic;
expectd = expectedDistance(instance, tau, 0, Q);
timeSpent = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read only one instance
%Windows:
%fid=fopen('D:\Documents\Seminario de Investigacion\VRP\Experiments\Instances\Novoa\data_thesis\i_5r1.dat', 'rt');
%Linux:
fid=fopen('/media/DATA_/Documents/Seminario de Investigacion/VRP/Experiments/Instances/Novoa/data_thesis/i_5r1.dat', 'rt');
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
tau = [0 5 3 2 1 4 0];
tic;
expectd = expectedDistance(instance, tau, 0, instance.Q);
timeSpent = toc;

%% Run rollout algorithm
tau = [0 1 2 3 4 5 0];
tic;
state = [0, instance.Q, -99, -99, -99, -99 , -99];
rapolicy = rollout( tau, instance, state );
timeSpent = toc;
