function [sumD] = avgDemand()

    x1= 2:4;
    x2 = 1:3;
    x3 = 1:3;
    x4 = 1:3;
    x5 = 1:3;
    sumD = zeros(1,11);%7:17    
    for i =1:3
        for j=1:3            
            for k=1:3
                for l = 1:3
                    for m=1:5
                        sum = x1(i)+x2(j)+x3(k)+x4(l)+x5(m);
                        sumD(sum-6) = sumD(sum-6) + 1;
                    end
                end                
            end
        end
    end
end