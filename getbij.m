function [bij] = getbij(x, i, j, t)
if(j==1)
    if(t(i) <= x && x < t(i+1))
        bij = 1;
        return;
    else
        bij = 0;
        return;
    end
end

bij = (x - t(i)) / (t(i+j) - t(i)) * getbij(x, i, j-1, t)...
        + (t(i+j+1) - x) / (t(i+j+1) - t(i+1)) * getbij(x, i+1, j-1, t);

end

