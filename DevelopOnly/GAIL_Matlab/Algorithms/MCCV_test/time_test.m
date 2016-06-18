counts = [10 100 1e3 2e3 3e3 5e3 8e3 1e4 2e4 3e4 5e4 8e4 ]
   % 1e5 1e5 2e5 3e5 5e5 8e5 1e6];

n = numel(counts);

fn = {};

for i = 1:n 
    fn{i} = @() gmean(counts(i));
end 


time = ones(1, n);

gmean(1e4);

for i = 1:n
    tic,fn{i}();,time(i) = toc;
end
    
plot(counts, time./counts, 'o')