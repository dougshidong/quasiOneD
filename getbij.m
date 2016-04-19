function [bij] = getbij(x, i, j, t)
if(j==0)
    if(t(i) <= x && x < t(i+1))
        bij = 1;
        return;
    else
        bij = 0;
        return;
    end
end

h = getbij(x, i, j-1, t);
k = getbij(x, i+1, j-1, t);
bij = 0;
if(h ~= 0)
    bij = bij + (x - t(i)) / (t(i+j) - t(i)) * h;
end
if(k ~= 0)
    bij = bij + (t(i+j+1) - x) / (t(i+j+1) - t(i+1)) * k;
end

return;
end