function [status,ME] = setupParpool(threads)
ME.message = '';
status = 0;
p = gcp('nocreate');
if isempty(p)
    try
        parpool(threads);
    catch ME
        status = 1;
    end
elseif p.NumWorkers ~= threads
    try
        delete(p);
        parpool(threads);
    catch ME
        status = 1;
    end
end
end