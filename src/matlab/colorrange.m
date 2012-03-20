function colors = colorrange(n)

colors = zeros(n,3);

h = 1/(n-1);

colors(:,3) = transpose(0:h:1);
%colors(:,2) = transpose(0:h:1);
colors(:,1) = transpose(1:(-h):0);

end