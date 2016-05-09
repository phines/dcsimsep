function randseed(seed)
% seed the random number generator in a way that will probably be unique
% for each call to this function. This is particularly helpful
% when doing parallel calls that might occur simultaneously
% aka MATLAB metaphysics (c)Edu Cotilla-Sanchez
% if seed is provided, just use this number

if nargin==0
    if isunix
        [~,pid] = system('echo $PPID');
        pid = str2num(pid);
        [~,hostname] = system('hostname');
    else
        pid = randi(1e5);
    end
    seed = (now-floor(now))*1e8 + pid + sum(single(hostname(:)));
end

s = RandStream.create('mt19937ar','seed',seed);
RandStream.setGlobalStream(s);
