c = 3*10^8 # Speed of light
f = 1e9 # Frequency
λ = c/f # Wavelength 

# Radar observation model (noisy power, phase, delay, doppler)

# decrease noise_std with certain action factors


# Example radar return model
function radar_return_power(λ,R; Pt=1, Gt=1, Gr=1,rcs=1, L=1, noise_std=0)
    # Pt: Transmit power
    # Gt: Transmit gain
    # Gr: Receive gain
    # rcs: Radar cross section
    # L: Losses
    return Pt * Gt * Gr * rcs * λ^2 / (16 * π^3 * R^4 * L) + noise_std*randn()
end

function radar_return_delay(R; noise_std=0)
    # R: Range to target (m)
    # Modify to scpecify a range gate with τ and Tipp
    return 2*R/c + noise_std*randn()
end

function radar_range(T_delay; noise_std=0)
    # T_delay: Delay (s)
    return c*T_delay/2 + noise_std*randn()
end

function radar_view_angle(x,y,radar_location; noise_std=0)
    return atand(x-radar_location[1], y-radar_location[2]) + noise_std*randn()
end

function radar_doppler(state; noise_std=0)
    # state: (x, y, vx, vy)    

    x, y, vx, vy = state
    r_hat = [x,y]
    v_hat = [vx, vy]
    Ψ = dot(r_hat *r_hat) / (norm(r_hat) * norm(v_hat)) # Angle between r and v
    V = sqrt(vx^2 + vy^2)                               # Velocity magnitude
    
    return 2*V*cos(Ψ)/λ + noise_std*randn()

end

function radar_range_resolution(τ; noise_std=0)
    # τ: Pulse width (s)
    return c*τ/2 + noise_std*randn()
end

function radar_unambiguous_range(T_ipp; noise_std=0)
    # T_ipp: interpulse period (s)
    return c*T_ipp/2 + noise_std*randn()
end

function radar_meas_2_obs_state(measurements; noise_std=0)
    # measurements: (range, delay, doppler, power)
    
    power, view_angle, delay, doppler_freq = measurements

    R_obs = radar_range(delay)

    


    return (px, py, vx, vy)
end