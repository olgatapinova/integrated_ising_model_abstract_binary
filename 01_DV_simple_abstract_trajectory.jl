function simulate_v_x_error_rate(r0, N, η, T, N1, N2, L, ϵ1, ϵ2, N_group1, N_group2, timemax, dLcrossed)
    # Parameters:
    # T - temperature
    # η - inhibition
    # ϵ1,2 - biases per target
    # L - threshold
    # N_group1,2 - number of spins in each group of spins (per target)
    # N1,2 - initial number of the active spins
    # r0 - units of spin-flipping rate
    # timemax - maximum time of the simulation
    # dLcrossed - distance to cut if finish far away from the target (% of L)

    ## Initialization

    N1_vec = Array{Int16,1}(); # active neurons (spins ON) in subgroup I
    N2_vec = Array{Int16,1}(); # active neurons (spins ON) in subgroup II

    t_vec = Array{Float64,1}(); # time points
    v_vec = Array{Float64,1}(); # difference in number of spins from different subgroups
    x_vec = Array{Float64,1}(); # position in time

    v = (N1 - N2) / N; # define the initial difference

    t = 0; # first time point
    x = 0; # initial position

    v_vec = push!(v_vec, v); 
    t_vec = push!(t_vec, t); # Initialization
    x_vec = push!(x_vec, x);
    N1_vec = push!(N1_vec, N1); 
    N2_vec = push!(N2_vec, N2); # Initialization

    ## Simulation
    
    while abs(x_vec[end]) < L && t < timemax

        # v is given already
        # calculate the rates
        r1_on  = r0*(N_group1 - N1) / (1 + exp((-2v + η - ϵ1) / T)); # group 1: ++1
        r1_off = r0*N1 / (1 + exp((2v-η + ϵ1) / T)); # group 1: --1

        r2_on  = r0*(N_group2 - N2) / (1 + exp((2v+ η - ϵ2) / T)); # group 2: ++1
        r2_off = r0*N2 / (1 + exp((-2v-η + ϵ2) / T)); # group 2: --1
        
        rate_tot = r1_on + r1_off + r2_on + r2_off; # total rate

        rand_time = 1 - rand(); # get a random number for gillespie
        dt = (1/rate_tot) * log(1/rand_time); # calculate next time interval dt

        rand_flip = 1 - rand(); # get a second random number

        # calculate next step: one spin flip occurs
        if rand_flip <= (r1_on / rate_tot)
            N1 += 1;
        elseif rand_flip <= ((r1_on + r1_off) / rate_tot)
            N1 -= 1;
        elseif rand_flip <= ((r1_on + r1_off + r2_on) / rate_tot)
            N2 += 1;
        else
            N2 -= 1;
        end

        t += dt; # update time point
        v = (N1 - N2) / N; # calculate new difference

        # update vectors
        v_vec = push!(v_vec, v); 
        t_vec = push!(t_vec, t); 
        x_vec = push!(x_vec, x_vec[end] + v*dt);
        N1_vec = push!(N1_vec, N1); 
        N2_vec = push!(N2_vec, N2);
        
    end
    
    #cut if finish far away from the threshold
    if abs(x_vec[end]) - L > L*dLcrossed
        
        if sign(x_vec[end]) > 0
            x_vec[end] = L;
            dx = L - x_vec[end-1];
        else
            x_vec[end] = -L;
            dx = -L - x_vec[end-1];
        end
        t_vec[end] = t_vec[end-1] + dx / v_vec[end];
    end
    
    return t_vec, v_vec, x_vec
    
end
