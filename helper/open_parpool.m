% open parpool
if isempty(gcp('nocreate'))
    parpool(24,'SpmdEnabled',false);
end;