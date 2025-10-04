# Simple runner script
include("../src/dc1d.jl")
include("../src/sens_1d_dc.jl")
include("../src/dc_inversion_main.jl")

# Execute the inversion
println("Starting DC resistivity inversion...")
rho, h, abbytwo, apro, aprc, err, data, model, par = main_dc_inversion()
println("Inversion completed!")

# Plot results - FIXED: Added the missing 'data' parameter
println("Generating plots...")
plot_results(rho, h, abbytwo, apro, aprc, err, data, model, par)
println("Done!")

