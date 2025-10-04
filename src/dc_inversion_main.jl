using DelimitedFiles
using LinearAlgebra
using Plots

# Include the function files
include("dc1d.jl")
include("sens_1d_dc.jl")

function main_dc_inversion()
    # Load synthetic data (replace with actual field data for real applications)
    current_dir = @__DIR__
    data_path = joinpath(current_dir, "..", "examples", "model1.dat")
    indata = readdlm(data_path)

    abbytwo = indata[:, 1]  # AB/2 values (electrode spacing)
    apro = log.(indata[:, 2])  # Log of observed apparent resistivity

    # Starting model (initial guess)
    rho = [500.0, 100.0]  # Resistivities of two layers (Ω·m)
    h = [50.0]  # Thickness of first layer (m)
    par = vcat(rho, h)  # Combined parameter vector

    # Noise level assumption (50% relative error)
    noise = 0.5
    ntr = 10  # Number of iterations

    # Forward calculation for initial model
    rhoa = dc1d(rho, h, abbytwo)
    aprc = log.(rhoa)  # Log of calculated apparent resistivity

    # Calculate sensitivity matrix
    sens_apr = sens_1d_dc(rho, h, abbytwo)

    # Prepare data for inversion
    dobs = apro  # Observed data (log apparent resistivity)
    dcomp = aprc  # Computed data (log apparent resistivity)
    A = sens_apr  # Sensitivity matrix

    # Discrepancy vector (difference between observed and computed)
    de = dobs - dcomp
    len_rhoa = length(rhoa)
    cnph = dobs  # Normalization factor

    # Weighting based on noise assumption
    weight = ones(length(dobs))
    de = de ./ (noise .* weight)
    A = A ./ (noise .* weight)

    # Initialize inversion parameters
    lem = 10.0  # Initial damping parameter
    minlm = 0.0  # Minimum damping parameter
    rms = sqrt(sum(((dobs - dcomp) .^ 2) ./ (cnph .^ 2)) / length(dobs))
    err = Float64[]  # To store RMS error history

    # Arrays to store inversion results
    data = zeros(length(rhoa), ntr)
    model = zeros(length(par), ntr)

    # Iterative inversion process
    rms1 = 1e20  # Initialize previous RMS error

    for iter in 1:ntr
        if rms <= rms1  # If solution is converging
            # Solve the regularized least-squares problem
            I_mat = Matrix{Float64}(I, size(A, 2), size(A, 2))
            dm = (A' * A + lem * I_mat) \ (A' * de)

            println("Iteration $iter: λ=$lem, RMS=$rms")

            # Update model parameters
            cdm = par .* exp.(dm)
            par = cdm
            rho = par[1:length(rho)]
            h = par[length(rho)+1:end]

            # Forward calculation with updated model
            rhoa = dc1d(rho, h, abbytwo)
            aprc = log.(rhoa)
            sens_apr = sens_1d_dc(rho, h, abbytwo)

            # Update data and sensitivity matrices
            dobs = apro
            dcomp = aprc
            A = sens_apr
            de = dobs - dcomp

            # Apply weighting - FIXED VERSION
            de = de ./ (noise .* weight)
            # Correct way to broadcast division for matrix A
            A = A ./ (noise .* weight)  # This automatically broadcasts correctly

            # Update RMS error
            rms1 = rms
            rms = sqrt(sum(((dobs - dcomp) .^ 2) ./ (cnph .^ 2)) / length(dobs))
            push!(err, rms)

            # Reduce damping parameter for next iteration
            lem *= 0.5

            if lem < minlm
                lem = minlm
            end
        else
            # If solution is diverging, increase damping
            iii = 0
            while rms > rms1 && iii <= 10
                lem /= 0.5
                println("Divergence correction $iii: λ=$lem, RMS=$rms")

                # Solve with increased damping
                I_mat = Matrix{Float64}(I, size(A, 2), size(A, 2))
                dm = (A' * A + lem * I_mat) \ (A' * de)

                # Update model
                cdm = par .* exp.(dm)
                par = cdm
                rho = par[1:length(rho)]
                h = par[length(rho)+1:end]

                # Forward calculation
                rhoa = dc1d(rho, h, abbytwo)
                aprc = log.(rhoa)
                sens_apr = sens_1d_dc(rho, h, abbytwo)

                # Update matrices
                dobs = apro
                dcomp = aprc
                A = sens_apr
                de = dobs - dcomp

                # Apply weighting - FIXED VERSION
                de = de ./ (noise .* weight)
                A = A ./ (noise .* weight)  # Simple broadcasting works

                # Update RMS error
                rms1 = rms
                rms = sqrt(sum(((dobs - dcomp) .^ 2) ./ (cnph .^ 2)) / length(dobs))
                iii += 1
            end
        end

        # Store results for this iteration
        data[:, iter] = rhoa
        model[:, iter] = par
    end

    # Store final RMS error
    push!(err, rms)

    return rho, h, abbytwo, apro, aprc, err, data, model, par
end

function plot_results(rho, h, abbytwo, apro, aprc, err, data, model, par)
    # Create convergence plot
    p1 = plot(err, yaxis=:log, marker=:circle, color=:black, linewidth=2, markersize=6,
        title="Convergence History", xlabel="Iteration Number", ylabel="RMS Error",
        grid=true, gridstyle=:dash, minorgrid=true)

    # Create data fitting plot
    p2 = plot(abbytwo, exp.(apro), seriestype=:scatter, marker=:circle, color=:blue,
        markersize=6, label="Observed",
        xaxis=:log, yaxis=:log,
        title="Data Fit: Observed vs. Calculated Apparent Resistivity",
        xlabel="AB/2 (m)", ylabel="ρa (Ω·m)",
        grid=true, gridstyle=:dash, minorgrid=true)
    plot!(abbytwo, exp.(aprc), color=:red, linewidth=2, label="Calculated")

    # Model visualization - FIXED: Use the final parameters directly
    NN = length(rho)

    # Use the final parameters that were returned, not from the model matrix
    rho_final = par[1:NN]
    h_final = par[NN+1:end]

    # Prepare model for plotting
    resis = Float64[]
    depth1 = Float64[]
    depth_sum = 0.0

    # Handle the case where we have multiple layers
    if length(h_final) > 0
        h_temp = vcat(h_final, h_final[end])  # Add a final layer thickness (use last thickness)
    else
        h_temp = [25.0]  # Default if no thickness parameters
    end

    for kk in 1:length(rho_final)
        temp1 = [rho_final[kk], rho_final[kk]]
        if kk <= length(h_temp)
            depth_sum += h_temp[kk]
        else
            depth_sum += 10.0  # Default thickness if not specified
        end
        temp2 = [depth_sum, depth_sum]
        append!(resis, temp1)
        append!(depth1, temp2)
    end

    depth_plot = vcat(0.0, depth1[1:end-1])

    # True model (for comparison - adjust based on your synthetic model)
    # If you don't have a true model, you can comment this out
    resisTrue = [100.0, 100.0, 10.0, 10.0]  # Adjust based on your true model
    depthTrue = [0.0, 20.0, 20.0, 40.0]     # Adjust based on your true model

    # Create the model comparison plot
    p3 = plot(resisTrue, -depthTrue, color=:blue, linewidth=3, label="True Model",
        seriestype=:steppost, xlabel="Resistivity (Ω·m)", ylabel="Depth (m)",
        title="Resistivity Model Comparison", grid=true, gridstyle=:dash)
    plot!(resis, -depth_plot, color=:red, linewidth=3, linestyle=:dash,
        label="Inverted Model", seriestype=:steppost,
        ylims=(-maximum(depth_plot) * 1.1, 0))

    # Subplot for data fitting
    p4 = plot(abbytwo, exp.(apro), seriestype=:scatter, marker=:circle, color=:blue,
        markersize=6, label="Observed",
        xaxis=:log, yaxis=:log,
        title="Data Fit",
        xlabel="AB/2 (m)", ylabel="ρa (Ω·m)",
        grid=true, gridstyle=:dash, minorgrid=true)
    plot!(abbytwo, exp.(aprc), color=:red, linewidth=2, label="Calculated")

    # Display individual plots
    display(p1)
    display(p2)

    # Create combined plot layout
    p_combined = plot(p3, p4, layout=(1, 2), size=(1200, 500))
    display(p_combined)
    # Explicitly display the combined plot
    # display(p_combined)

    # Keep the window open (choose one method)
    println("Press Enter to exit...")
    readline()

    return p1, p2, p_combined
end

# Run the inversion and plotting
if abspath(PROGRAM_FILE) == @__FILE__
    rho, h, abbytwo, apro, aprc, err, data, model, par = main_dc_inversion()
    plot_results(rho, h, abbytwo, apro, aprc, err, model, par)
end