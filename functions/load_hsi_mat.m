function [data, gt_map] = load_hsi_mat(file_name)
%LOAD_HSI_MAT Load a hyperspectral dataset with flexible variable names.

S = load(file_name);

if isfield(S, 'data')
    data = S.data;
else
    data = [];
    names = fieldnames(S);
    for i = 1:numel(names)
        value = S.(names{i});
        if isnumeric(value) && ndims(value) == 3
            data = value;
            break;
        end
    end
    if isempty(data)
        error('No H x W x B data cube found in %s.', file_name);
    end
end

if isfield(S, 'map')
    gt_map = S.map;
elseif isfield(S, 'mask')
    gt_map = S.mask;
elseif isfield(S, 'gt')
    gt_map = S.gt;
elseif isfield(S, 'GT')
    gt_map = S.GT;
else
    gt_map = [];
end
end
