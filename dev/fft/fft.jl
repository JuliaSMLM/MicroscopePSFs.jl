module PureFFT

"""
    split_even_odd(x::AbstractVector{T}) where T
Split sequence into even and odd parts using reshape.
"""
function split_even_odd(x::AbstractVector{T}) where T
    N = length(x)
    reshaped = reshape(x, 2, N÷2)
    return view(reshaped, 1, :), view(reshaped, 2, :)
end

"""
    butterfly_combine(evens::AbstractVector, odds::AbstractVector, twiddles::AbstractVector)
Perform butterfly operations on vectors.
"""
function butterfly_combine(evens::AbstractVector, odds::AbstractVector, twiddles::AbstractVector)
    odd_terms = odds .* twiddles
    lower = evens .+ odd_terms
    upper = evens .- odd_terms
    return vcat(lower, upper)
end

"""
    generate_twiddles(N::Integer)
Generate twiddle factors for N-point FFT.
"""
function generate_twiddles(N::Integer)
    k = 0:(N÷2-1)
    cis.(-2π * k / N)
end

"""
    base_case_fft(x::AbstractVector{<:Number})
Handle the base case of length-1 vectors without scalar indexing.
"""
function base_case_fft(x::AbstractVector{<:Number})
    reshape(complex.(x), :)
end

"""
    fft(x::AbstractVector{<:Number}) -> Vector{Complex{Float64}}

Compute the Discrete Fourier Transform using a pure functional radix-2 FFT.
This implementation uses only GPU-compatible operations.
Input length must be a power of 2.
"""
function fft(x::AbstractVector{<:Number})
    N = length(x)
    
    # Verify input length is power of 2
    if !ispow2(N)
        throw(ArgumentError("Input length must be a power of 2"))
    end
    
    # Base case
    N == 1 && return base_case_fft(x)
    
    # Split into even and odd parts using reshape
    evens, odds = split_even_odd(x)
    
    # Recursively compute FFT of even and odd subsequences
    even_fft = fft(evens)
    odd_fft = fft(odds)
    
    # Generate twiddle factors 
    twiddles = generate_twiddles(N)
    
    # Combine using vectorized butterfly operations
    return butterfly_combine(even_fft, odd_fft, twiddles)
end

"""
    ifft(X::AbstractVector{<:Number}) -> Vector{Complex{Float64}}

Compute the Inverse Discrete Fourier Transform.
Uses only GPU-compatible operations.
"""
function ifft(X::AbstractVector{<:Number})
    N = length(X)
    conj.(fft(conj.(X))) ./ N
end

# Convenience method for real inputs
fft(x::AbstractVector{<:Real}) = fft(complex.(x))

export fft, ifft

end # module

module FFTShift

using CUDA
using Zygote
using ..PureFFT

"""
    fft_shift(x::AbstractVector{T}, shift::Real) where T
Shift a signal by a real-valued amount using FFT phase multiplication.
The shift parameter can be fractional.
"""
function fft_shift(x::AbstractVector{T}, shift::Real) where T
    N = length(x)
    k = vcat(0:N÷2-1, -N÷2:-1)
    
    # Create phase factors
    phase = cis.(-2π * k * shift / N)
    phase = T <: Complex ? phase : complex.(phase)
    
    # Move phase to GPU if input is on GPU
    phase_device = x isa CuArray ? CuArray(phase) : phase
    
    # Perform FFT, multiply by phase, then IFFT
    X = fft(x)
    X_shifted = X .* phase_device
    return real.(ifft(X_shifted))
end

export fft_shift

end # module

# Example usage and gradient testing:
using CUDA
using Zygote
using Plots
using .FFTShift: fft_shift

# Create a test signal (Gaussian pulse)
N = 64
x = collect(1:N)
signal = exp.(-(x .- N/4).^2 ./ 50)

# Test the shift function
shift = 10.0
shifted = fft_shift(signal, shift)

# Loss function
loss(s, x) = sum(abs2, fft_shift(x, s))

# Compute gradient more cleanly
∇shift = Zygote.gradient(s -> loss(s, signal), shift)[1]

# Basic visualization
p = plot(signal, label="Original", title="Signal Shifting")
plot!(shifted, label="Shifted")
println("CPU Gradient: ", ∇shift)

# Test GPU version
signal_gpu = CuArray(signal)
shifted_gpu = fft_shift(signal_gpu, shift)
∇shift_gpu = Zygote.gradient(s -> loss(s, signal_gpu), shift)[1]
println("GPU Gradient: ", ∇shift_gpu)

display(p)