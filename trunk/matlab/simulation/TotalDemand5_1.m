ds = zeros(1, factorial(6));
s = [1 2 3 4 5 6];
ss = perms(s); 
cas=0; 
for t=1: size(ss,1)
	policy = [0 ss(t,1); 0 ss(t,2); 0 ss(t,3); 0 ss(t,4); 0 ss(t,5); 0 ss(t,6)];
	points = [0 0; 1 0; 1 0.5; 1 1; 0 1; 0 0.5];
	Q = 2.4;
    disDem = zeros(1,11);
	demVals= 7:17;
    disDis = zeros(1,243);
	disVals= 1:104;
	cases = 0;
	tot = 0;
	var = 0;
	distot = 0;
    for x1 =2:4
        for x2=1:3            
            for x3=1:3
                for x4 =1:3
                   for x5=2:4  
					    xt = x1+x2+x3+x4+x5;
						disDem(xt-6) = disDem(xt-6) +1;
                        tot = tot + xt;
						var = var + (xt -12)^2;
						inst = [0 x1 x2 x3 x4 x5];
						cases = cases +1;
						pos = 1;
						pospol = 1;
						dist = 0;
						q = Q;
                        while (sum(inst) ~= 0)
						
							if ((q == 0) || (policy(pospol,1) == 1))                            
								dist = dist + sqrt( points(pos,1)^2 +points(pos,2)^2 );
							    pos = 1;
							    q = Q;
						    end
							nextpospol = pospol + 1;
							if (nextpospol == 7) 
							     nextpospol = 1;
							end
                            nextpos = policy(nextpospol,2);
							dist = dist + sqrt( (points(pos,1)-points(nextpos,1))^2 + (points(pos,2)-points(nextpos,2))^2);
							if ( q >= inst(nextpos)) 
							    q = q -inst(nextpos);
								inst(nextpos) = 0; 
								pos = nextpos;
								pospol = nextpospol;
							else
                                inst(nextpos) = inst(nextpos) - q;
								q = 0;
								pos = nextpos;
								%pospol = nextpospol;
							end
				        
						end
						disDis(cases) = dist;
						distot = distot + dist;
                    end
                end
               
            end
        end
    end
	tot=tot/cases;
	var=var/cases;
	distot=distot/cases;
	cas = cas + 1;
	ds(cas) = distot;
end