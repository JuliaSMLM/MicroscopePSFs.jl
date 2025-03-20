using Pkg
Pkg.activate("dev")
using Revise
using MicroscopePSFs.Zernike
using CairoMakie

# Check noll2nl and nl2noll functions
function check_noll_conversion(max_j::Integer=21)
    for j in 1:max_j
        n, l = noll2nl(j)
        j_check = nl2noll(n, l)
        if j != j_check
            println("Mismatch: j=$j, (n=$n, l=$l) -> j_check=$j_check")
        end
    end
end



# Function to check RMS normalization of Zernike polynomials
function check_zernike_normalization(max_j::Integer=15, grid_size::Integer=201)
    println("\nChecking RMS normalization of Zernike polynomials:")
    println("-----------------------------------------------------")
    println("j\tn\tl\tRMS\t\tDiff from 1.0")
    println("-----------------------------------------------------")
    
    for j in 1:max_j
        n, l = noll2nl(j)
        
        # Create ZernikeCoefficients with a single coefficient set to 1
        coeffs = ZernikeCoefficients(j)
        coeffs.mag[1] = 0.0
        
        coeffs.mag[j] = 1.0
        
        # Generate the polynomial field
        field = evaluate_pupil(coeffs, grid_size)
        
        # Calculate RMS value over valid points (inside unit circle)
        total_sum = 0.0
        count = 0
        
        # Create normalized coordinate grid
        xs = ys = range(-1, 1, length=grid_size)
        
        for i in 1:grid_size, j_idx in 1:grid_size
            x, y = xs[i], ys[j_idx]
            ρ = sqrt(x^2 + y^2)
            if ρ <= 1.0
                total_sum += abs2(real(field[j_idx, i]))
                count += 1
            end
        end
        
        # Calculate RMS
        rms = sqrt(total_sum / count)
        diff = rms - 1.0
        
        # Print results with formatting
        println("$j\t$n\t$l\t$(round(rms, digits=6))\t$(round(diff, digits=6))")
    end
end

# Function to evaluate a 2D grid of polynomial values for specific n, l using evaluate_pupil
function evaluate_polynomial_grid_nl(n::Integer, l::Integer, grid_size::Integer=201)
    # Convert n,l to Noll index j
    j = nl2noll(n, l)

    # Create ZernikeCoefficients with a single coefficient set to 1
    coeffs = ZernikeCoefficients(j)

    # Set the specific coefficient to 1
    coeffs.mag[1] = 0.0 # was set to 1 by default
    coeffs.mag[j] = 1.0

    # Use evaluate_pupil to generate the field
    field = evaluate_pupil(coeffs, grid_size)

    # Convert complex field to real values for visualization
    # Since we only set magnitude coefficients, we can just take the real part
    result = real.(field)
    return result
end

# Function to plot a grid of Zernike polynomials organized by n and l
function plot_zernike_nl_grid(max_n::Integer=4)
    println("Plotting Zernike polynomials organized by n and l...")
    println("Blue is low, Red is high")
    println("Displayed with positive y up, to match wikipedia")
    fig = Figure(size=(1000, 180 * max_n))

    burd = reverse(cgrad(:RdBu))

    for n in 0:max_n
        # For each n, we want to plot l values from -n to n in steps of 2
        # But l can only have the same parity as n
        l_start = -n + (mod(-n, 2) != mod(n, 2) ? 1 : 0)
        l_values = l_start:2:n

        # Create a row for this n value
        row = n + 1  # Row index (1-based)

        for (col_idx, l) in enumerate(l_values)
            # Get j value
            j = nl2noll(n, l)


            # Create axis
            ax = Axis(fig[row, col_idx],
                aspect=DataAspect(),
                title="j=$j, (n=$n, l=$l)",
                xticklabelsize=0,
                yticklabelsize=0,
                yreversed=false)

            # Generate and plot polynomial
            z_data = evaluate_polynomial_grid_nl(n, l, 201)

            # Plot range that is symmetric around zero
            max_val = maximum(abs.(z_data))
            min_val = -max_val


            # Use a diverging colormap centered at zero
            heatmap!(ax, range(-1, 1, length=size(z_data, 1)),
                range(-1, 1, length=size(z_data, 2)),
                z_data',
                colorrange=(min_val, max_val),
                colormap=burd)

            # Hide axis
            hidedecorations!(ax)
        end
    end

    return fig
end

# Check noll2nl and nl2noll functions
check_noll_conversion(21)

# Check Zernike polynomial normalization
check_zernike_normalization(15)

# Generate plot of Zernike polynomials organized by n and l
plot_zernike_nl_grid(5)




