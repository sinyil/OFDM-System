function start = detectpreamble(signal,preamble,threshold)
% Correlates signal with preamble including a energy
% normalization
% 
% If it detects a peak larger than the threshold it 
% returns the index. Otherwise it returns a vector with
% the correlation sequence.
%

p = preamble(:);
plen = length(p);
slen = length(signal);

c = zeros(slen,1);

for i=1:slen-plen
    
    fragment    = signal(i:i+plen-1);
    normalize   = sum(abs(fragment).^2);
    
    corval      = sum((p.*fragment))^2/normalize; 
    
    c(i)        = corval;
    
    if corval > threshold        
        start = i + plen;
        return
    end
    
end
disp('Frame start not found.')
start = c;
return