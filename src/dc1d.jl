function dc1d(rho, z, c1c2)
    rho1 = rho[1]
    rho2 = rho[2]
    p1p2 = 2

    # Calculate geometric factors
    r1 = (c1c2 .- p1p2) ./ 2
    r2 = (c1c2 .+ p1p2) ./ 2
    r3 = r2
    r4 = r1

    # Reflection coefficient (scalar)
    k = (rho2 - rho1) / (rho1 + rho2)

    # Geometric factor
    p = (1 ./ r1 .- 1 ./ r2 .- 1 ./ r3 .+ 1 ./ r4) .^ -1

    # Summation of infinite series - FIXED: use z.^2 for element-wise power
    sumterm = zeros(length(c1c2))
    for m in 1:45
        km = k ^ m  # This is fine - both scalars
        # Use broadcasting with z.^2 instead of z^2
        term = (1 ./ sqrt.(r1 .^ 2 .+ 4 * m ^ 2 * z .^ 2) .-
                1 ./ sqrt.(r2 .^ 2 .+ 4 * m ^ 2 * z .^ 2) .-
                1 ./ sqrt.(r3 .^ 2 .+ 4 * m ^ 2 * z .^ 2) .+
                1 ./ sqrt.(r4 .^ 2 .+ 4 * m ^ 2 * z .^ 2))
        sumterm .+= km .* term
    end

    # Calculate apparent resistivity
    Rhom = rho1 .* (1 .+ 2 .* p .* sumterm)
    return Rhom
end