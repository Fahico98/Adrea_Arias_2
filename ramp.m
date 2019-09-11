
function toReturn = ramp(t)
    toReturn = zeros(1, length(t));
    for i = 1 : length(t)
        if t(i) >= 0
            toReturn(i) = t(i);
        else
            toReturn(i) = 0;
        end
    end
end