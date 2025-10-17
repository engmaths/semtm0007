# Julia v1.10.10
# Status `~/ddpm-data-gen/Project.toml`
#   [15e1cf62] NPZ v0.4.3
#   [91a5bcdd] Plots v1.41.1
#   [90137ffa] StaticArrays v1.9.15
#   [ed894a53] WaterLily v1.5.2
#   [9a3f8284] Random
#   [10745b16] Statistics v1.10.0

using WaterLily, StaticArrays, Plots, Statistics, NPZ, Random, CUDA

const USE_CUDA = false  # CUDA.functional()  # CUDA doesn't seem to be working unfortunately
const MEM = USE_CUDA ? CuArray : Array

if USE_CUDA
    println("Using CUDA")
else
    println("Not using CUDA")
end

T = Float32
default_centers = SVector{2, T}.([(0.0, 5.2), (-3.7, 1.9), (-0.5, 0.8), (0.5, -2.2), (-3.0, -2.3), (-1.5, -4.5)])  # [10rand(SVector{2, T}) .- 5 for _ in 1:6]
default_radii = [1.5, 0.9, 1.2, 1.3, 0.8, 1.2]  # (3rand(T, length(centers)) .+ 3) ./ 4

function circle(n, m; Re = 550, U = 1, mem = MEM, T = Float32, centers, radii)
    R, x0 = m ÷ 18, m ÷ 2
    # 6 circles with random centers x,y ∈ [5,5] and radius r ∈ [0.75,1.5]
    body = sum(1:length(centers)) do i
        center = centers[i]
        radius = radii[i]
        AutoBody((x, t) -> √sum(abs2, x .- x0 - center * R) - radius * R)
    end
    # make a simulation
    return Simulation((n, m), (U, 0), R; ν = U * R / Re, body, mem, T)
end

# make a simulation and run it
function run_simulation(; Re = 550, U = 1, centers = default_centers, radii = default_radii)
    sim = circle(3 * 2^7, 2^8; Re, U, centers, radii)

    duration = 30

    # sim_gif!(
    #     sim, duration = duration, clims = (-5, 5), remeasure = false, plotbody = true, axis = ([], false),
    #     cfill = :seismic, legend = false, border = :none
    # )

    t₀ = round(WaterLily.sim_time(sim))

    # Compute the initial vorticity
    WaterLily.@inside sim.flow.σ[I] = WaterLily.curl(3, I, sim.flow.u) * sim.L / sim.U
    σ = [copy(sim.flow.σ)]

    @time for tᵢ in range(t₀, t₀ + duration; step = 0.1)
        # Step the simulation
        WaterLily.sim_step!(sim, tᵢ; remeasure = false)
        # Compute the vorticity
        WaterLily.@inside sim.flow.σ[I] = WaterLily.curl(3, I, sim.flow.u) * sim.L / sim.U
        # Store a copy of the vorticity field at this time step
        push!(σ, copy(sim.flow.σ))
        # Logging
        # println("tU/L=", round(tᵢ, digits = 4), ", Δt=", round(sim.flow.Δt[end], digits = 3))
    end

    n, m = size(σ[1]) .- 2  # exclude ghost cells
    bins_x = 24
    bins_y = 16
    bin_size = 16
    x_offset = 1
    y_offset = 1
    sigma = stack(σ)
    sigma_mean = [mean(σ[k][x_offset .+ (i - 1) * bin_size .+ (1:bin_size), y_offset .+ (j - 1) * bin_size .+ (1:bin_size)]) for i in 1:bins_x, j in 1:bins_y, k in eachindex(σ)]
    return (sigma[2:(end - 1), 2:(end - 1), :], sigma_mean)
end

function run_simulation_gif(; Re = 550, U = 1, centers = default_centers, radii = default_radii)
    sim = circle(3 * 2^7, 2^8; Re, U, centers, radii)

    duration = 30

    return sim_gif!(
        sim, duration = duration, clims = (-5, 5), remeasure = false, plotbody = true, axis = ([], false),
        cfill = :seismic, legend = false, border = :none
    )
end

function generate_data_base_case()
    Re = 550.0
    U = 1.0
    radii = default_radii
    (sigma, sigma_mean) = run_simulation(; Re, U, radii)
    npzwrite(joinpath(@__DIR__, "vorticity_base_case.npz"); sigma, sigma_mean, Re, U, radii)
    return
end

function generate_data_vary_Re_U(N = 500)
    rng = Random.Xoshiro(123)
    Re = 50 .+ rand(rng, N) .* (2000 - 50)  # Random Re in [50, 2000]
    U = 0.5 .+ rand(rng, N) .* (2.5 - 0.5)  # Random U in [0.5, 2.5]
    radii = [default_radii for _ in 1:N]
    data = Array{Float32, 3}[]
    for (i, (Re_i, U_i, radii_i)) in enumerate(zip(Re, U, radii))
        _, sigma_mean = run_simulation(; Re = Re_i, U = U_i, radii = radii_i)
        push!(data, sigma_mean)
        println("Completed simulation $i / $N with Re=$Re_i, U=$U_i (vary_Re_U)")
    end
    sigma_mean = stack(data)
    radii_stacked = stack(radii)
    npzwrite(joinpath(@__DIR__, "vorticity_vary_Re_U.npz"); sigma_mean, Re, U, radii = radii_stacked)
    return
end

function generate_data_vary_radii(N = 500)
    rng = Random.Xoshiro(321)
    Re = 50 .+ rand(rng, N) .* (2000 - 50)  # Random Re in [50, 2000]
    U = 0.5 .+ rand(rng, N) .* (2.5 - 0.5)  # Random U in [0.5, 2.5]
    radii = [default_radii .- 0.2 .+ rand(rng, length(default_radii)) .* 0.4 for _ in 1:N]
    data = Array{Float32, 3}[]
    for (i, (Re_i, U_i, radii_i)) in enumerate(zip(Re, U, radii))
        _, sigma_mean = run_simulation(; Re = Re_i, U = U_i, radii = radii_i)
        push!(data, sigma_mean)
        println("Completed simulation $i / $N with Re=$Re_i, U=$U_i (vary_radii)")
    end
    sigma_mean = stack(data)
    radii_stacked = stack(radii)
    npzwrite(joinpath(@__DIR__, "vorticity_vary_radii.npz"); sigma_mean, Re, U, radii = radii_stacked)
    return
end
