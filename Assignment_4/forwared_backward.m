function [smoothed, var] = forwared_backward(true, noise, alpha)

%Step 1. Forward raws

raw_forw = zeros(length(true), length(true));
raw_forw(1, :) = noise(1, :);

for i = 2:length(true)
    raw_forw(i, :) = raw_forw(i - 1, :) + alpha*(noise(i, :) - raw_forw(i - 1, :));
end
%Step 2. Backward raws
raw_back = zeros(length(true), length(true));
raw_back(end, :) = raw_forw(end, :);

for i = (length(true) - 1):-1:1
    raw_back(i, :) = raw_back(i + 1, :) + alpha*(raw_forw(i, :) - raw_back(i + 1, :));
end
%Step 3. Forward columns
col_forw = zeros(length(true), length(true));
col_forw(:,1) = raw_back(:,1);

for i = 2:length(true)
    col_forw(:, i) = col_forw(:,i - 1) + alpha*(raw_back(:, i) - col_forw(:,i - 1));
end

%Step 4. Backwards columns

col_back = zeros(length(true), length(true));
col_back(end, :) = col_forw(end, :);

for i = (length(true) - 1):-1:1
    col_back(:, i) = col_back(:, i + 1) + alpha*(col_forw(:, i) - col_back(:,i + 1));
end

smoothed = col_back;

var = sum((reshape((true-smoothed),...
    [1,length(true)*length(true)])).^2)/(length(true)*length(true)-1);
end

