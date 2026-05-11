module DisplacementUncertainty

using Statistics

export displacement_variance_pitch_roll, propagate_rho_uncertainty,
       combine_variances, estimate_pitch_roll_uncertainty

"""
Uncertainty propagation for pitch/roll errors through displacement calculations.
Used with weighted total least squares (TLS) to properly weight structure function fits.
"""

"""
    displacement_variance_pitch_roll(rng1, rng2, pitch1, pitch2, roll1, roll2, σ_pitch, σ_roll)

Compute displacement variance contribution from pitch/roll measurement uncertainty.

Uses first-order error propagation: Var(f) = (∂f/∂x)² Var(x) + (∂f/∂y)² Var(y)

# Arguments
- `rng1, rng2`: range gate distances (m)
- `pitch1, pitch2`: pitch angles at times t1, t2 (rad)
- `roll1, roll2`: roll angles at times t1, t2 (rad)
- `σ_pitch`: pitch uncertainty standard deviation (rad)
- `σ_roll`: roll uncertainty standard deviation (rad)

# Returns
- `(δdr2, δdz2)`: variance contributions to dr² and dz²

# Physics
Displacements in beam geometry:
- dz = rng₂·cos(θ₂)cos(φ₂) - rng₁·cos(θ₁)cos(φ₁)  (vertical)
- dx = X + rng₂·(-sin(θ₂)) - rng₁·(-sin(θ₁))        (along-wind)
- dy = Y + rng₂·cos(θ₂)sin(φ₂) - rng₁·cos(θ₁)sin(φ₁) (cross-wind)

where θ=pitch, φ=roll, X,Y are advection by mean wind.
"""
function displacement_variance_pitch_roll(rng1, rng2, pitch1, pitch2, roll1, roll2, σ_pitch, σ_roll)
    # Trig functions
    cp1, sp1 = cos(pitch1), sin(pitch1)
    cr1, sr1 = cos(roll1), sin(roll1)
    cp2, sp2 = cos(pitch2), sin(pitch2)
    cr2, sr2 = cos(roll2), sin(roll2)

    # Partial derivatives for dz = rng·cos(pitch)·cos(roll)
    # ∂dz/∂pitch = -rng·sin(pitch)·cos(roll)
    # ∂dz/∂roll = -rng·cos(pitch)·sin(roll)
    ∂dz_∂pitch1 = -rng1 * sp1 * cr1
    ∂dz_∂pitch2 = -rng2 * sp2 * cr2
    ∂dz_∂roll1 = -rng1 * cp1 * sr1
    ∂dz_∂roll2 = -rng2 * cp2 * sr2

    # Partial derivatives for dx component: rng·(-sin(pitch))
    # ∂dx/∂pitch = -rng·cos(pitch)
    ∂dx_∂pitch1 = rng1 * cp1
    ∂dx_∂pitch2 = -rng2 * cp2

    # Partial derivatives for dy = rng·cos(pitch)·sin(roll)
    # ∂dy/∂pitch = -rng·sin(pitch)·sin(roll)
    # ∂dy/∂roll = rng·cos(pitch)·cos(roll)
    ∂dy_∂pitch1 = rng1 * (-sp1) * sr1
    ∂dy_∂pitch2 = rng2 * (-sp2) * sr2
    ∂dy_∂roll1 = rng1 * cp1 * cr1
    ∂dy_∂roll2 = rng2 * cp2 * cr2

    # Variance propagation (assuming independent pitch and roll errors at each time)
    var_dz = (∂dz_∂pitch1^2 + ∂dz_∂pitch2^2) * σ_pitch^2 +
             (∂dz_∂roll1^2 + ∂dz_∂roll2^2) * σ_roll^2

    var_dx = (∂dx_∂pitch1^2 + ∂dx_∂pitch2^2) * σ_pitch^2

    var_dy = (∂dy_∂pitch1^2 + ∂dy_∂pitch2^2) * σ_pitch^2 +
             (∂dy_∂roll1^2 + ∂dy_∂roll2^2) * σ_roll^2

    # Total displacement dr² = dx² + dy² + dz²
    # For small uncertainties, variance adds
    δdr2 = var_dx + var_dy + var_dz
    δdz2 = var_dz

    return (δdr2, δdz2)
end


"""
    propagate_rho_uncertainty(dr2, dz2, var_dr2, var_dz2)

Propagate dr² and dz² uncertainties to ρ = r^(2/3)·(1 - (dz/dr)²/4)

# Arguments
- `dr2, dz2`: squared displacement values
- `var_dr2, var_dz2`: their variances from uncertainty propagation

# Returns
- `var_ρ`: variance of ρ (missing if dr2 invalid)

# Physics
Structure function fitting uses the coordinate:
    ρ = r^(2/3)·(1 - (dz/dr)²/4)

where r = √(dr²). This coordinate accounts for anisotropy in the velocity structure
function due to the vertical component being suppressed relative to horizontal.
"""
function propagate_rho_uncertainty(dr2, dz2, var_dr2, var_dz2)
    if dr2 <= 0 || !isfinite(dr2) || !isfinite(dz2)
        return missing
    end

    # ρ = dr2^(1/3) · (1 - dz2/(4·dr2))
    factor = dr2^(1/3)
    correction = 1 - dz2/(4*dr2)

    # Partial derivatives via chain rule
    # ∂ρ/∂dr2 = (1/3)·dr2^(-2/3)·(1 - dz2/(4·dr2)) + dr2^(1/3)·dz2/(4·dr2²)
    ∂ρ_∂dr2 = (1/3)*dr2^(-2/3)*correction + factor*dz2/(4*dr2^2)

    # ∂ρ/∂dz2 = dr2^(1/3)·(-1/(4·dr2)) = -dr2^(-2/3)/4
    ∂ρ_∂dz2 = -dr2^(-2/3)/4

    # Variance propagation
    var_ρ = ∂ρ_∂dr2^2 * var_dr2 + ∂ρ_∂dz2^2 * var_dz2

    return var_ρ
end


"""
    combine_variances(propagated_var, sample_var, n_samples)

Combine propagated instrument variance with empirical binned sample variance.

# Arguments
- `propagated_var`: variance from uncertainty propagation (instrument error)
- `sample_var`: empirical variance within bin (sampling + turbulence)
- `n_samples`: number of samples averaged in bin

# Returns
- `total_var`: combined variance estimate for binned mean

# Statistics
For a binned mean of n independent samples:
    Var(mean) = Var_instrument + Var_sampling/n

where Var_sampling includes both measurement scatter and real turbulent variability.
"""
function combine_variances(propagated_var, sample_var, n_samples)
    # Instrument variance affects every sample
    # Sample variance is reduced by averaging (central limit theorem)
    if n_samples > 1
        return propagated_var + sample_var / n_samples
    else
        # Single sample: both variances contribute fully
        return propagated_var + sample_var
    end
end


"""
    estimate_pitch_roll_uncertainty(pitch, roll)

Estimate pitch/roll measurement uncertainty from short-term data variability.

Uses standard deviation of the data as a proxy for measurement uncertainty.
This is appropriate when the true pitch/roll should be slowly varying (ship motion
on scales >> 1 second), so high-frequency variability indicates sensor noise.

# Arguments
- `pitch`: pitch angle time series (rad)
- `roll`: roll angle time series (rad)

# Returns
- `(σ_pitch, σ_roll)`: uncertainty estimates in radians

# Notes
Falls back to conservative 2° uncertainty if no valid data (MDV-only case).
Floors at 0.1° (VectorNav spec) to avoid underestimating uncertainty.
"""
function estimate_pitch_roll_uncertainty(pitch::AbstractVector, roll::AbstractVector)
    # Filter to finite values
    pitch_valid = filter(isfinite, pitch)
    roll_valid = filter(isfinite, roll)

    if !isempty(pitch_valid) && !isempty(roll_valid)
        # Use std dev as uncertainty estimate
        σ_pitch_deg = max(std(pitch_valid) * 180/π, 0.1)  # floor at VN spec
        σ_roll_deg = max(std(roll_valid) * 180/π, 0.1)

        return (deg2rad(σ_pitch_deg), deg2rad(σ_roll_deg))
    else
        # No valid VN data - use conservative default
        # (acknowledges missing motion correction)
        return (deg2rad(2.0), deg2rad(2.0))
    end
end

end # module DisplacementUncertainty
