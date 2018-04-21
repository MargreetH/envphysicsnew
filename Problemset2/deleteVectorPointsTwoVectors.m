
%This function removes both values from 2 vectors if one of the vectors has
%an instance of value at that position.
function [y,z] = deleteVectorPointsTwoVectors(a,b, value)

if length(a) ~= length(b)
   return; 
end

counter = 1;

for i=1:1:length(a)
    
      if a(i) == value | b(i) == value %don't use the value
          %do nothing
      else %do use the value
          y(counter) = a(i);
          z(counter) = b(i);
          counter = counter + 1;
      end
end


end