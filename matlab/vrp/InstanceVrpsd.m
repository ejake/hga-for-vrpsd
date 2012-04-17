classdef InstanceVrpsd %< handle
    %InstanceVrpsd Summary an single VRPSD instance
    %   Detailed explanation goes here
    
    properties
        n = 5; %Number of customers
        f=1; %Factor for vehicles Q
        %The mean of the demand for this instance
        meanDemand=0;
        Q; %Vehicle capacity
        veh;
        Cust = [];
        %Demand Probability Distribution for each customer
        DemPD = [];
        %Location for each customers
        LL = [];
        %Depot location
        depot=[0,0]
        %distance matrix
        d;
    end
    
    methods
        function obj = InstanceVrpsd(numCustomers, factorF, DemandsPD, Locations)
            obj.n = numCustomers;
            obj.f = factorF;
            obj.meanDemand = sum((DemandsPD(:,1)+ DemandsPD(:,2))/2)/obj.n;
            obj.Q = (obj.meanDemand*obj.n)/(1+obj.f);
            obj.veh = Vehicle(obj.Q,0,obj.Q);
            obj.DemPD = DemandsPD;
            obj.LL = Locations;
            for i =1:obj.n
                auxCust = Customer(Locations(i,:), DemandsPD(i,:));
                obj.Cust = [obj.Cust; auxCust];
            end
            %assesing distance
            obj.d = zeros(obj.n+1,obj.n+1);
            for i =1:obj.n%distance to depot
                obj.d(1,i+1) = Distance(Locations(i,:),obj.depot);
                obj.d(i+1,1) = obj.d(1,i+1);
            end            
            for i = 1:obj.n%distance between customers
                for j=i+1:obj.n
                    obj.d(i+1,j+1) = Distance(Locations(i,:),Locations(j,:));
                    obj.d(j+1,i+1) = obj.d(i+1,j+1);
                end                
            end
            %end assesing distances
        end
    end
    
end

