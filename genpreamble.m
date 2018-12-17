function preamble = genpreamble(length)
% Task 2.3.1


preamble        = zeros(length,1);

% init state
preamble(1:8)   = ones(8,1);

% cycle through states
for i=9:length
    next        = [1 0 1 1 1 0 0 0 ]*preamble(i-8:i-1); % gen next state
    nextgf2     = mod(next,2); % convert to GF(2)
    preamble(i) = nextgf2;    
end
