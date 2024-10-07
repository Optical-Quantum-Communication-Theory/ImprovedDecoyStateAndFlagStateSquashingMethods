function index = sub2indPlus(vecSize,vec)
vec = num2cell(vec);
index = sub2ind(vecSize,vec{:});
end

