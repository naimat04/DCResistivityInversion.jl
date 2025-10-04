function sens_1d_dc(rho, h, abbytwo)
    """
    Sensitivity Calculation by Brute force method for 1D DC data
    
    Parameters:
    rho : Vector of resistivity values for each layer (Ω·m)
    h : Vector of layer thicknesses (m)
    abbytwo : Vector of electrode spacing (AB/2 values in meters)
    
    Returns:
    sens_apr : Sensitivity matrix of log(apparent resistivity) w.r.t model parameters
    """
    
    # Forward computation for first assumed model
    rhoa = dc1d(rho, h, abbytwo)
    
    # Sensitivity calculation
    nm = length(rho) + length(h)  # no. of model parameters
    n_abbytwo = length(abbytwo)
    sens_apr = zeros(n_abbytwo, nm)
    
    for ii in 1:nm
        if ii <= length(rho)
            # Sensitivity w.r.t resistivity parameters
            res = rho[ii]
            rho[ii] = 1.05 * res
            
            rhoa1 = dc1d(rho, h, abbytwo)
            
            # Sensitivity w.r.t log(resistivity)
            sens_apr[:, ii] = (log.(rhoa1) .- log.(rhoa)) / 0.05
            
            # Restore original value
            rho[ii] = res
        else
            # Sensitivity w.r.t thickness parameters
            jj = ii - length(rho)  # index for thickness array
            thick = h[jj]
            h[jj] = 1.05 * thick
            
            rhoa1 = dc1d(rho, h, abbytwo)
            
            # Sensitivity w.r.t log(thickness)
            sens_apr[:, ii] = (log.(rhoa1) .- log.(rhoa)) / 0.05
            
            # Restore original value
            h[jj] = thick
        end
    end
    
    return sens_apr
end