using Revise
using Pkg; Pkg.activate(".")

using FFTW
using Plots

"derivative by FFT, assuming sampling rate of 1."
function fft_derivative(f::AbstractVector{Real})
    N = length(f) 
    k = fftfreq( N ) * 2π  # compute the angular frequencies (wavenumbers)
    # for even N exclude Nyquist frequency according to SG Johnson
    # https://math.mit.edu/~stevenj/fft-deriv.pdf
    if iseven(N)
        k[N/2] = 0
    end
    F_deriv = 1im * k .* fft(f)    # multiply FFT by i*k in the frequency domain
    return real( ifft(F_deriv) )     # inverse FFT to get the derivative in time domain
end

# Example usage
x = (1:1024) / 1024
f = sin.(9.3*2π * x)  # Sample time series (e.g., sine wave)
# taper the ends with an n-point quarter sine^2 wave taper
taper(n) = @. sin( π/2 * (0:n-1)/n )^2
f_taper(f,n) = [f[1:n].*taper(n); f[n+1:end-n]; f[end-n+1:end].*reverse(taper(n))]
f_tapered = f_taper(f, 24)

f_exact_deriv = 9.3*2π*cos.(9.3*2π * x)
# 
dx = x[2] - x[1]
f_deriv = fft_derivative(f) / dx

# Plot the original function and its derivative
plot(x, f, label="Original function")
plot!(x, f_deriv, label="First derivative", marker=".", xlim=(0.25, 0.75))
plot!(x, f_exact_deriv, label="Exact derivative")

plot(x, f_deriv - f_exact_deriv, label="Exact derivative")
# works well even near edges when N/2 wavenumber is properly zeroed out.
