function p = setup_pool(numWorkers)

p = gcp('nocreate');
if isempty(p)
    tic; p=parpool(numWorkers); fprintf('Pool initiation: '); toc; p
elseif p.NumWorkers ~= numWorkers
    delete(gcp('nocreate'));
    tic; p=parpool(numWorkers); fprintf('Pool initiation: '); toc; p
else
    fprintf('Pool already exists.'); p
end