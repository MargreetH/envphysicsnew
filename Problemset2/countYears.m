function y = countYears(v)


y=1;
curYear = v(1);
for i = 1:1:length(v)
    if v(i) ~= curYear
        y = y + 1;
        curYear = v(i);
    else
        
    end
    
end

end