function vec = ind2subPlus(vecSize,index)
vec = cell(size(vecSize));
[vec{:}] = ind2sub(vecSize,index);
vec = cell2mat(vec);
end