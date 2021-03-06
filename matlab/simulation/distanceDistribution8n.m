function [frecDis, disDis, distances] = distanceDistribution8n(policy, instance, cyclic)

    %cyclic = logical(cyclic);
    %numDemSum = 0;
    minDem = 0;
    maxDem = 0;
    for i = 1: instance.n        
        minDem = minDem + instance.Cust(i).PD(1);
        maxDem = maxDem + instance.Cust(i).PD(2);        
    end
    numDemSum = maxDem - minDem + 1;
    frecDis = zeros(1,numDemSum);
    disDis = zeros(1,numDemSum);
    %cases = 0;
    distot = 0;
    iter = 1;
    distances = zeros(1,(instance.Cust(1).PD(2)-instance.Cust(1).PD(1) + 1)*(instance.Cust(2).PD(2)-instance.Cust(2).PD(1) + 1)*(instance.Cust(3).PD(2)-instance.Cust(3).PD(1) + 1)*(instance.Cust(4).PD(2)-instance.Cust(4).PD(1) + 1)*(instance.Cust(5).PD(2)-instance.Cust(5).PD(1) + 1)*(instance.Cust(6).PD(2)-instance.Cust(6).PD(1) + 1)*(instance.Cust(7).PD(2)-instance.Cust(7).PD(1) + 1)*(instance.Cust(8).PD(2)-instance.Cust(8).PD(1) + 1));
    for x1 = instance.Cust(1).PD(1):instance.Cust(1).PD(2) % Demand of customer 1
        for x2 = instance.Cust(2).PD(1):instance.Cust(2).PD(2) % Demand of customer 2
            for x3 = instance.Cust(3).PD(1):instance.Cust(3).PD(2) % Demand of customer 3
                for x4 = instance.Cust(4).PD(1):instance.Cust(4).PD(2) % Demand of customer 4
                    for x5 = instance.Cust(5).PD(1):instance.Cust(5).PD(2) % Demand of customer 5
                        for x6 = instance.Cust(6).PD(1):instance.Cust(6).PD(2) % Demand of customer 6
                            for x7 = instance.Cust(7).PD(1):instance.Cust(7).PD(2) % Demand of customer 7
                                for x8 = instance.Cust(8).PD(1):instance.Cust(8).PD(2) % Demand of customer 8
                                    inst = [x1 x2 x3 x4 x5 x6 x7 x8];
                                    totDem = sum(inst);
                                    pos = 0;
                                    pospol = 1;
                                    dist = 0;
                                    q = instance.Q;
                                    while (sum(inst) ~= 0)%asses distance
                                        if ((q == 0) || (policy(pospol).a == 1))%Proactive restocking when q is equal to 0
                                            dist = dist + instance.d(1,pos+1);
                                            pos = 0;
                                            q = instance.Q;
                                        end

                                        %if pos == 0
                                        %    nextpospol = 1;
                                        %else
                                        %    nextpospol = pospol + 1;
                                        %end
                                        %if (nextpospol > instance.n) 
                                        %    nextpospol = 1;
                                        %end

                                        nextpos = policy(pospol).m;
                                        dist = dist + instance.d(pos+1, nextpos+1);
                                        if ( q >= inst(nextpos)) 
                                            q = q -inst(nextpos);
                                            inst(nextpos) = 0; 
                                            pos = nextpos;
                                            pospol = nextpos + 1;
                                        else
                                            inst(nextpos) = inst(nextpos) - q;
                                            q = 0;
                                            pos = nextpos;
                                            %if(cyclic)
                                            %    pospol = nextpospol; %uncomment if run cyclic tour
                                            %end
                                        end

                                    end
                                    if(pos ~= 0)%Add the distance to return to the depot 
                                        dist = dist + instance.d(1,pos+1);
                                    end
                                    frecDis(totDem - minDem + 1) = frecDis(totDem - minDem + 1) + 1;
                                    disDis(totDem - minDem + 1) = dist;
                                    distot = distot + dist;
                                    distances(iter) = dist;
                                    iter = iter+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end