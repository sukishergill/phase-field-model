% Generate large collection of cos 

num_modes = 100;              
u = zeros(Grid.N(1), Grid.N(2));              
rng(30);                      

for i = 1:num_modes
    kx = (rand() * 20) - 10;  
    ky = (rand() * 20) - 10;

    phase = rand() * 2*pi;    
    
    amplitude = rand();       
    
    % Accumulate the 2D cosine wave
    u = u + amplitude * cos(kx * Grid.xx + ky * Grid.yy + phase);
end

u = (u - min(u(:))) / (max(u(:)) - min(u(:))) - 0.5;
