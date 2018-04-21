function [k,labda,error] = findWeibull(v,m,maxError)
%FINDWEIBULL Finds the weibull distribution parameters k and labda using
%the variance v and the mean m of the distribution.
%   

correctValue = sqrt(v)/m;
first_interval_k = [0.0001 10];
dk = 1/10;
kvector = first_interval_k(1):first_interval_k(2)*dk:first_interval_k(2);

counter = 0;
error = 10000;
maxIterations = 1000;



while (error > maxError) && (counter < maxIterations);

    values = zeros(1,length(kvector));
    for i=1:1:length(kvector)
        values(i) = sqrt(gamma(1+2/kvector(i))/gamma(1+1/kvector(i))^2-1);
    end
    
    errors = abs(values-correctValue);
    [error,I] = min(errors);  
    dk = dk * 0.5;
    k= kvector(I); 
    
    if I == 1
       kvector = kvector(1):dk:kvector(2); 
    elseif I == length(kvector)
       kvector = kvector(I-1):dk:kvector(end);
    else
        kvector = kvector(I-1):dk:kvector(I+1);
    end

counter = counter + 1;
end

disp(counter);
labda = m / gamma(1/k + 1);

end


