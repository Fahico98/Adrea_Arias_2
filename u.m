
function toReturn = u(t)
    toReturn = zeros(1, length(t));
    for i = 1 : length(t)
        if t(i) >= 0
            toReturn(i) = 1;
        else
            toReturn(i) = 0;
        end
    end
end