function FLAG=test_memery(memorySize)
if nargin<1
    memorySize=4e9;
end

if ismac
    % Code to run on Mac platform
    [r,w] = unix('free | grep Mem');
    stats = str2double(regexp(w, '[0-9]*', 'match'));
    %     memsize = stats(1)/1e6;
    %     freemem = (stats(3)+stats(end))/1e6;
    memsize = stats(1);
    freemem = (stats(3)+stats(end));
elseif isunix
    % Code to run on Linux platform
elseif ispc
    % Code to run on Windows platform
    [userview systemview]=memory;
    freemem=userview.MemAvailableAllArrays;
else
    disp('Platform not supported')
end


if freemem>memorySize
    FLAG=true;
else
    FLAG=false;
end
end
