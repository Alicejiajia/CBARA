function [state] = single_state_init(N, Ts, x_start, y_start, x_end, y_end)

    x1 = linspace(x_start, x_end, N).';
    y1 = linspace(y_start, y_end, N).';
    start1 = [x_start, y_start];
    ed1 = [x_end, y_end];
    [vx1, vy1] = cpu_vxy(start1, ed1, Ts, N); 
    vx1 = linspace(vx1, vx1, N).';
    vy1 = linspace(vy1, vy1, N).';
    
    state_ideal = [x1, y1, vx1, vy1];
    std_pos = 0.5; 
    std_vel = 0.3;
    
    noise_pos = std_pos * randn(N, 2);
    noise_vel = std_vel * randn(N, 2);

    state = state_ideal + [noise_pos, noise_vel];
    
    state(1, :) = state_ideal(1, :); 

    plot(state(:,1), state(:,2), '--d', 'color', [0.4940 0.1840 0.5560], 'linewidth', 1.5, 'MarkerIndices', [1 N], 'markersize', 4);
    hold on;
end
