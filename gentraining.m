function training = gentraining(length)

training        = zeros(length,1);

% init state
training(1:8)   = ones(8,1);

% cycle through states
for i=9:length
    next        = [1 0 1 0 1 0 0 1 ]*training(i-8:i-1); % gen next state
    nextgf2     = mod(next,2); % convert to GF(2)
    training(i) = nextgf2;    
end
