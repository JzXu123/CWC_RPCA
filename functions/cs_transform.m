function result_transform = cs_transform(result, threshold)
%CS_TRANSFORM Capped-square score stretching used by the original demo.

result_transform = result;
lower = result < threshold;
result_transform(lower) = (result(lower) / threshold).^2;
result_transform(~lower) = 1;
end
