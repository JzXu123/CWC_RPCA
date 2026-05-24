function x = normalize_minmax(x)
%NORMALIZE_MINMAX Scale numeric input to [0, 1].

x = double(x);
lo = min(x(:));
hi = max(x(:));
if hi > lo
    x = (x - lo) / (hi - lo);
else
    x = zeros(size(x));
end
end
