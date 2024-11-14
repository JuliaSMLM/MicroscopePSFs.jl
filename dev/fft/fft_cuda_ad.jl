using CUDA
using Zygote

function fftfreq(n::Integer)
    if iseven(n)
        k = vcat(0:n÷2-1, -n÷2:-1)
    else
        k = vcat(0:(n-1)÷2, -(n-1)÷2:-1)
    end
    return float(k) ./ n
end

function fourier_shift(shift, input)
    n = length(input)
    # Create frequencies and ensure they match input type
    freq = fftfreq(n)
    freq = convert(Vector{eltype(float(input))}, freq)
    if input isa CUDA.CuArray
        freq = CuArray(freq)
    end
    
    F = my_fft(input)
    # Ensure shift is converted to match input type
    shift_t = convert(eltype(float(input)), shift)
    phase_shift = exp.(-2π * im * freq * shift_t)
    F_shifted = F .* phase_shift
    return real.(my_ifft(F_shifted))
end

function my_fft(x)
    n = length(x)
    if n <= 1
        return x
    end
    
    # Split into even and odd indices
    even = x[1:2:end]
    odd = x[2:2:end]
    
    even = my_fft(even)
    odd = my_fft(odd)
    
    half_n = n ÷ 2
    # Convert angle to match input type
    angle = convert(eltype(float(x)), -2π / n)
    t = (0:half_n-1)
    if x isa CuArray
        t = CuArray(t)
    end
    factor = exp.(angle * im * t)
    
    a = even .+ factor .* odd
    b = even .- factor .* odd
    
    return vcat(a, b)
end

function my_ifft(x)
    n = length(x)
    if n <= 1
        return x
    end
    y = conj(my_fft(conj(x)))
    return y ./ n
end

# Test the implementation
sz = 64
a = rand(sz)
n_shift = 3.0  # Use float for consistency

# CPU version
b = fourier_shift(n_shift, a)
grads, = Zygote.jacobian(s -> fourier_shift(s, a), n_shift)

# GPU version
a_d = CuArray(a)
b_d = fourier_shift(n_shift, a_d)
grads_d, = Zygote.jacobian(s -> fourier_shift(s, a_d), n_shift)

# Compare gradients
println("CPU gradient: ", grads)
println("GPU gradient: ", Array(grads_d))