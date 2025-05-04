using StaticArrays
using LinearAlgebra
using Distributions
using Plots

# Abstract transmitter type
abstract type AbstractTransmitter end

# PointTransmitter
struct PointTransmitter <: AbstractTransmitter
    position::SVector{2, Float64}
    frequency::Float64
    power::Float64
    tx_time::Float64
    pulse_width::Float64
    steering_angle::Float64
    beamwidth::Float64
    is_isotropic::Bool
end

# Reflector
mutable struct Reflector
    position::SVector{2, Float64}
    pulse_width::Float64
    reflection_gain::Float64
    has_reflected::Bool
end

# Receiver
mutable struct Receiver
    position::SVector{2, Float64}
    received_power::Vector{Float64}
end

# Radar environment
mutable struct RadarEnvironment
    grid_x::Vector{Float64}
    grid_y::Vector{Float64}
    wavefield::Matrix{Float64}
    transmitters::Vector{PointTransmitter}
    reflectors::Vector{Reflector}
    receivers::Vector{Receiver}
    current_time::Float64
end

function RadarEnvironment(grid_x::Vector{Float64}, grid_y::Vector{Float64})
    wavefield = zeros(length(grid_x), length(grid_y))
    transmitters = PointTransmitter[]
    reflectors = Reflector[]
    receivers = Receiver[]
    current_time = 0.0
    return RadarEnvironment(grid_x, grid_y, wavefield, transmitters, reflectors, receivers, current_time)
end

function add_transmitter!(env::RadarEnvironment, tx::PointTransmitter)
    push!(env.transmitters, tx)
end

function rm_transmitter!(env::RadarEnvironment, tx::PointTransmitter)
    env.transmitters = filter(x -> x != tx, env.transmitters)
end

function add_reflector!(env::RadarEnvironment, refl::Reflector)
    push!(env.reflectors, refl)
end

function add_receiver!(env::RadarEnvironment, rx::Receiver)
    push!(env.receivers, rx)
end

function pulse_amplitude(tx_time::Float64, pulse_width::Float64, t::Float64)
    return pdf(Normal(tx_time, pulse_width), t)
end

function findnearest(grid::Vector{Float64}, val::Float64)
    return findmin(abs.(grid .- val))[2]
end

function snap_to_grid(grid::Vector{Float64}, pos::SVector{2, Float64})
    ix = findnearest(grid, pos[1])
    iy = findnearest(grid, pos[2])
    return SVector(grid[ix], grid[iy])
end

function propagate!(env::RadarEnvironment)
    env.wavefield .= 0.0

    new_transmitters = PointTransmitter[]

    for tx in env.transmitters
        if env.current_time < tx.tx_time
            continue
        end

        for (i, xi) in enumerate(env.grid_x)
            for (j, yj) in enumerate(env.grid_y)
                pos = SVector(xi, yj)
                r = norm(pos - tx.position)

                # --- Beam steering logic ---
                delta = pos - tx.position
                if r == 0
                    angle_diff = 0.0
                else
                    angle = atan(delta[2], delta[1]) * 180 / π
                    if angle < 0
                        angle += 360
                    end

                    angle_diff = abs(angle - tx.steering_angle)
                    if angle_diff > 180
                        angle_diff = 360 - angle_diff
                    end
                end

                # Check if point is inside beam
                if !tx.is_isotropic
                    if angle_diff > tx.beamwidth / 2
                        continue  # skip → outside beam
                    end
                end

                # --- Pulse propagation ---
                t_arrival = r / 300
                amplitude = pulse_amplitude(tx.tx_time, tx.pulse_width, env.current_time - t_arrival)

                if amplitude <= 1e-12
                    continue
                end

                propagated_power = tx.power * amplitude / (r == 0 ? 1.0 : r^2)
                env.wavefield[i, j] += propagated_power

                # Reflector check
                for refl in env.reflectors
                    if refl.has_reflected
                        continue
                    end

                    if pos == refl.position
                        if propagated_power > 1e-12
                            reflected_power = propagated_power * refl.reflection_gain
                            reflect_time = env.current_time

                            reflected_tx = PointTransmitter(
                                refl.position,
                                tx.frequency,
                                reflected_power,
                                reflect_time,
                                refl.pulse_width,
                                0.0,          # steering angle irrelevant
                                360.0,        # full beam
                                true          # isotropic flag
                            )
                            push!(new_transmitters, reflected_tx)

                            println("Reflection triggered at ", refl.position, " power = ", reflected_power)

                            refl.has_reflected = true
                        end
                    end
                end
            end
        end
    end

    append!(env.transmitters, new_transmitters)

    # Final receiver sampling (after propagation)
    for rx in env.receivers
        xi = findnearest(env.grid_x, rx.position[1])
        yi = findnearest(env.grid_y, rx.position[2])
        observed = env.wavefield[xi, yi]
        push!(rx.received_power, observed)
    end
end


function step!(env::RadarEnvironment, dt::Float64)
    env.current_time += dt
    propagate!(env)
end

function plot_wavefield(env::RadarEnvironment)
    log_wavefield = 10 .* log10.(env.wavefield' .+ 1e-12)

    heatmap(env.grid_x, env.grid_y, log_wavefield, 
        c=:viridis, xlabel="x (m)", ylabel="y (m)", 
        title="Radar Wavefield at t = $(round(env.current_time, digits=2)) s",
        clims=(-120, 0))
end