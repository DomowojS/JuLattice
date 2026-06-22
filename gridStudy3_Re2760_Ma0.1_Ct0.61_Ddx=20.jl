############################
## Main file for JuLattice #
############################
include("src/Plotter.jl")
include("src/Logger.jl")
include("src/TurbulenceModel.jl")
include("src/Kernel.jl")


using JLD2          # for saving last plot (cross-version compatible)
using MeshGrid, GLMakie
using .Plotter, .Logger
using .TurbulenceModel 
using .Kernel

function run_JuLattice()
    ####################################  Initialize  ####################################
    ##-------- User Settings --------##
    # Disc parameters
    D = 0.05        # diameter
    C_T = 0.61 #0.65 #S67 Rotor      # thrust coefficient

    # convert global C_T to local C_T from 1D momentum theory
    # C_T = 4a(1-a) => a = (1- sqrt(1-C_T)) / 2
    a = (1.0 - sqrt(1.0 - C_T)) / 2.0
    C_T_local = C_T / (1.0 - a)^2
    println("C_T = $(C_T) | C_T_local = $(C_T_local)")

     # Domainsize from cylinder validation
    length_X = 17.5 * D   
    length_Y = 10 * 0.023       
    length_Z = 10 * 0.023       

    # Grid spacing (physical units per lattice unit)
    # value from grid independence study
    delta_x         = 0.0025  #0.00092 

    # Fluid Settings 
    Kinematic_Viscosity = 1e-6                                       # m^2/s 
    reynoldsNumber =   2760 #2760                                    # Target Reynolds number
    Mach_Number = 0.1 # 0.05                                         # Target Mach number (Ma = U_lattice/c_s)
                                                                     # Keep Ma < 0.1 for incompressible flow!

    # Simulation Settings
    Simulation_Time = 60                                             # s
    
                              
   
    # Smagorinsky constant CS
    CS              = 1/3                  # CS ↑ = eddy viscosity ↑

    # Plot Requests (Flags)
    Plotvx = false;
    Plotmag = true;
    Plotdebug = false;
    Plotvorticity = false;
    vorticity_mode = :component # :component (ω_z / ω_y)   or   :magnitude (|ω|)

    ##-------- Compute LBM Parameters from Mach Number --------##
    lattice_speedOfSound    = 1.0 / sqrt(3)
    Inflow_Velocity         = reynoldsNumber * Kinematic_Viscosity / D
    speedOfSound            = Inflow_Velocity / Mach_Number
    delta_t                 = delta_x * lattice_speedOfSound / speedOfSound
    lattice_viscosity       = Kinematic_Viscosity * delta_t / (delta_x)^2
    lattice_inflow_velocity = Mach_Number * lattice_speedOfSound

    # nu_lattice = c_s² * (tau - 0.5) => tau = nu_lattice / c_s² + 0.5
    τ       = lattice_viscosity / (lattice_speedOfSound * lattice_speedOfSound) + 0.5
    omega   = 1.0 / τ

    

    # force Fx in stream direction
    # F_x_lat = -0.5 * C_T * (lattice_inflow_velocity^2)

    fluiddensity = 1.0 # lattice units
    simulationTime = ceil(Int, Simulation_Time / delta_t);  #lattice units
    
    ##-------- Convert user settings to lattice units --------##
    # Domain
    gridlengthX = ceil(Int, length_X / delta_x);
    gridlengthY = ceil(Int, length_Y / delta_x);
    gridlengthZ = ceil(Int, length_Z / delta_x);

    # Disc
    disc_radius = (D / 2) / delta_x
    # disc_x = Int(round((5.5 * 0.023) / delta_x)) + 1
    disc_x = Int(round((5.5 * D) / delta_x)) + 1
    disc_y = Int(round((length_Y / 2) / delta_x)) + 1
    disc_z = Int(round((length_Z / 2) / delta_x)) + 1
    # disc thickness in x-direction
    disc_thickness = 1


    # Reynolds Check:
    # Re_lattice = U*R/v -> should match Re_phys since quantities are scaled
    Re_lattice = floor(Int, ((lattice_inflow_velocity .* 2 .* disc_radius)/lattice_viscosity)) 
    Re_phys = Inflow_Velocity * D / Kinematic_Viscosity

    # Define Slice indices for plotting
    midY = 2 + Int(round((gridlengthY-2)/2))
    midZ = 2 + Int(round((gridlengthZ-2)/2))

    frontY = 2
    backY = gridlengthY-1

    ##-------- Probe Setup --------## 
    D_lat   = Int(round(D / delta_x))
    probe_labels    = ["2D", "4D", "6D", "8D", "10D"]
    probe_xs        = [disc_x + 2*D_lat, disc_x + 4*D_lat, disc_x + 6*D_lat, disc_x + 8*D_lat, disc_x + 10*D_lat]
    n_probes        = length(probe_xs)
    probe_z         = midZ
    probe_ys        = collect(2:gridlengthY-1)
    n_probe_y       = length(probe_ys)

    # sample_dt_phys      = 0.1   # sampling rate = 10Hz
    sample_dt_phys      = 0.01    # sampling rate = 100Hz
    log_dt_phys         = 1.0   # logging rate for csv-flush
    sample_interval     = max(1, round(Int, sample_dt_phys / delta_t))
    samples_per_flush   = max(1, round(Int, log_dt_phys / sample_dt_phys))

    sample_times        = zeros(samples_per_flush)
    cumulativ_count     = 0
    buf_ptr = 0

    # Sample buffers for logging
    sample_buf_u      = [zeros(n_probe_y, samples_per_flush) for _ in 1:n_probes]
    sample_buf_v      = [zeros(n_probe_y, samples_per_flush) for _ in 1:n_probes]
    sample_buf_w      = [zeros(n_probe_y, samples_per_flush) for _ in 1:n_probes]
    sample_buf_mean_u = [zeros(n_probe_y, samples_per_flush) for _ in 1:n_probes]
    sample_buf_mean_v = [zeros(n_probe_y, samples_per_flush) for _ in 1:n_probes]
    sample_buf_mean_w = [zeros(n_probe_y, samples_per_flush) for _ in 1:n_probes]
    sample_buf_std_u  = [zeros(n_probe_y, samples_per_flush) for _ in 1:n_probes]
    sample_buf_std_v  = [zeros(n_probe_y, samples_per_flush) for _ in 1:n_probes]
    sample_buf_std_w  = [zeros(n_probe_y, samples_per_flush) for _ in 1:n_probes]

    cumulativ_mean_u = [zeros(n_probe_y) for _ in 1:n_probes]
    cumulativ_mean_v = [zeros(n_probe_y) for _ in 1:n_probes]
    cumulativ_mean_w = [zeros(n_probe_y) for _ in 1:n_probes]
    cumulativ_M2_u   = [zeros(n_probe_y) for _ in 1:n_probes]
    cumulativ_M2_v   = [zeros(n_probe_y) for _ in 1:n_probes]
    cumulativ_M2_w   = [zeros(n_probe_y) for _ in 1:n_probes]
    # time estimation
    time_samples = zeros(11)


    # Inlet momentum coefficients for D3Q19 weights (velocity bounceback)
    inlet_add_face = (2.0 / (18.0 * lattice_speedOfSound^2)) * lattice_inflow_velocity
    inlet_add_edge = (2.0 / (36.0 * lattice_speedOfSound^2)) * lattice_inflow_velocity

    ##--------  classify nodes --------##
    ## create solid node mask
    # Array{Bool} instead of BitArray: single byte load in kernel loop vs bit-unpack
    is_solid  = fill(false, gridlengthX, gridlengthY, gridlengthZ)
    is_fluid  = fill(false, gridlengthX, gridlengthY, gridlengthZ)
    is_disc   = fill(false, gridlengthX, gridlengthY, gridlengthZ)
    
    # Create masks
    for x in 1:gridlengthX, y in 1:gridlengthY, z in 1:gridlengthZ
        # walls
        if y == 1 || y == gridlengthY || z == 1 || z == gridlengthZ
            is_solid[x,y,z] = true
            continue
        end

        # disc at disc_x, circular in y,z-plane with disc_radius
        dx = x - disc_x
        dy = y - disc_y
        dz = z - disc_z
        if abs(dx) < disc_thickness && sqrt(dy^2 + dz^2) <= disc_radius
            is_disc[x,y,z] = true
        end
    end
    # pre compute fluid range
    is_fluid[2:gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= true

    disc_nodes = findall(is_disc)

    for idx in disc_nodes
        is_fluid[idx] = false
    end

    # fluid mask
    n_fluid_nodes = sum(is_fluid)
    n_mnups_nodes = n_fluid_nodes

    ##-------- precompute wall BC index lists --------##
    # Per-face filtered (x,z) or (x,y) lists
    wall_front_x = Int[]; wall_front_z = Int[]   # y=1  face → writes to y=2
    wall_back_x  = Int[]; wall_back_z  = Int[]   # y=Ny face → writes to y=Ny-1
    wall_bot_x   = Int[]; wall_bot_y   = Int[]   # z=1  face → writes to z=2
    wall_top_x   = Int[]; wall_top_y   = Int[]   # z=Nz face → writes to z=Nz-1

    sizehint!(wall_front_x, gridlengthX * gridlengthZ)
    sizehint!(wall_front_z, gridlengthX * gridlengthZ)
    sizehint!(wall_back_x,  gridlengthX * gridlengthZ)
    sizehint!(wall_back_z,  gridlengthX * gridlengthZ)
    sizehint!(wall_bot_x,   gridlengthX * (gridlengthY - 2))
    sizehint!(wall_bot_y,   gridlengthX * (gridlengthY - 2))
    sizehint!(wall_top_x,   gridlengthX * (gridlengthY - 2))
    sizehint!(wall_top_y,   gridlengthX * (gridlengthY - 2))

    # y-faces: full z range — is_solid[x,2,z] is true at z=1/Nz corners → auto-skipped
    for z in 1:gridlengthZ, x in 1:gridlengthX
        if !is_solid[x, 2, z]
            push!(wall_front_x, x); push!(wall_front_z, z)
        end
        if !is_solid[x, gridlengthY-1, z]
            push!(wall_back_x, x); push!(wall_back_z, z)
        end
    end

    # z-faces: interior y only — y-wall corners are in the y-face lists above
    for y in 2:gridlengthY-1, x in 1:gridlengthX
        if !is_solid[x, y, 2]
            push!(wall_bot_x, x); push!(wall_bot_y, y)
        end
        if !is_solid[x, y, gridlengthZ-1]
            push!(wall_top_x, x); push!(wall_top_y, y)
        end
    end
    
    ##-------- Array Allocation --------##
    #Q = 19; #D3Q19
    # D3Q19 — f and fS are 4D arrays: f[q, x, y, z]
    # first dimension q is indexed with direction constants Q000..Q0PP (exported by Kernel):
    # f[Q000] = rest (0,0,0)
    # f[QM00], f[QP00] = x-axis (±1,0,0)
    # f[Q0M0], f[Q0P0] = y-axis (0,±1,0)
    # f[Q00M], f[Q00P] = z-axis (0,0,±1)
    # f[QMM0], f[QMP0], f[QPM0], f[QPP0] = xy-plane edges
    # f[QM0M], f[QM0P], f[QP0M], f[QP0P] = xz-plane edges
    # f[Q0MM], f[Q0MP], f[Q0PM], f[Q0PP] = yz-plane edges

    # f  = pre-collision populations; fS = post-collision push target
    # undef + parallel first-touch: each thread touches its own z-slice so Linux
    # places those pages on the local NUMA node, matching the collision_stream! access pattern.
    f  = Array{Float64}(undef, NQ, gridlengthX, gridlengthY, gridlengthZ)
    fS = Array{Float64}(undef, NQ, gridlengthX, gridlengthY, gridlengthZ)
    Threads.@threads for z in 1:gridlengthZ
        for y in 1:gridlengthY, x in 1:gridlengthX
            @inbounds for q in 1:NQ
                f[q,x,y,z]  = 0.0
                fS[q,x,y,z] = 0.0
            end
        end
    end

    # Initialise macroscopic variables — same first-touch pattern
    # rho        = Array{Float64}(undef, gridlengthX, gridlengthY, gridlengthZ)
    u          = Array{Float64}(undef, gridlengthX, gridlengthY, gridlengthZ)
    v          = Array{Float64}(undef, gridlengthX, gridlengthY, gridlengthZ)
    w          = Array{Float64}(undef, gridlengthX, gridlengthY, gridlengthZ)
    velocityX   = Array{Float64}(undef, gridlengthX, gridlengthY, gridlengthZ)
    velocityY   = Array{Float64}(undef, gridlengthX, gridlengthY, gridlengthZ)
    velocityZ   = Array{Float64}(undef, gridlengthX, gridlengthY, gridlengthZ)
    velocityMag = Array{Float64}(undef, gridlengthX, gridlengthY, gridlengthZ)
    vortZ       = Array{Float64}(undef, gridlengthX, gridlengthY, gridlengthZ)
    vortY       = Array{Float64}(undef, gridlengthX, gridlengthY, gridlengthZ)
    Threads.@threads for z in 1:gridlengthZ
        for y in 1:gridlengthY, x in 1:gridlengthX
            @inbounds begin
                # rho[x,y,z]        = 1.0
                u[x,y,z]          = 0.0
                v[x,y,z]          = 0.0
                w[x,y,z]          = 0.0
                velocityX[x,y,z]  = 0.0
                velocityY[x,y,z]  = 0.0
                velocityZ[x,y,z]  = 0.0
                velocityMag[x,y,z] = 0.0
                vortZ[x,y,z]      = 0.0
                vortY[x,y,z]      = 0.0
            end
        end
    end

    ##--------  Initialize distribution functions FLUID NODES and SOLID NODES  --------##
    for z in 1:gridlengthZ
        for y in 1:gridlengthY
            for x in 1:gridlengthX

                # Equilibirum Initialisation
                # ux = is_solid[x, y, z] ? 0.0 : lattice_inflow_velocity
                # uy = 0.0
                # uz = 0.0

                # disturbed Initialisation rand()-> random output[0,1] -0.5 -> shifts to [-0.5,0.5]
                ux = is_solid[x, y, z] ? 0.0 : lattice_inflow_velocity
                uy = is_solid[x, y, z] ? 0.0 : (rand() - 0.5) * 0.02 * lattice_inflow_velocity
                uz = is_solid[x, y, z] ? 0.0 : (rand() - 0.5) * 0.02 * lattice_inflow_velocity

                rho_init = fluiddensity

                # Pre-compute polynomial factors
                ux2 = ux * ux
                uy2 = uy * uy
                uz2 = uz * uz

                Pm_u = 1 - 3*ux + 3*ux2
                P0_u = 1 - 1.5*ux2
                Pp_u = 1 + 3*ux + 3*ux2

                Pm_v = 1 - 3*uy + 3*uy2
                P0_v = 1 - 1.5*uy2
                Pp_v = 1 + 3*uy + 3*uy2

                Pm_w = 1 - 3*uz + 3*uz2
                P0_w = 1 - 1.5*uz2
                Pp_w = 1 + 3*uz + 3*uz2

                # Push scheme: Rest particle (0,0,0) - weight 1/3
                f[Q000,x,y,z] = rho_init * P0_u * P0_v * P0_w / 3.0

                # Push scheme: Face neighbors - weight 1/18
                f[QM00,x,y,z] = rho_init * Pm_u * P0_v * P0_w / 18.0
                f[QP00,x,y,z] = rho_init * Pp_u * P0_v * P0_w / 18.0

                f[Q0M0,x,y,z] = rho_init * P0_u * Pm_v * P0_w / 18.0
                f[Q0P0,x,y,z] = rho_init * P0_u * Pp_v * P0_w / 18.0

                f[Q00M,x,y,z] = rho_init * P0_u * P0_v * Pm_w / 18.0
                f[Q00P,x,y,z] = rho_init * P0_u * P0_v * Pp_w / 18.0

                # Push scheme: Edge neighbors - weight 1/36
                # XY-plane edges
                f[QMM0,x,y,z] = rho_init * Pm_u * Pm_v * P0_w / 36.0
                f[QMP0,x,y,z] = rho_init * Pm_u * Pp_v * P0_w / 36.0
                f[QPM0,x,y,z] = rho_init * Pp_u * Pm_v * P0_w / 36.0
                f[QPP0,x,y,z] = rho_init * Pp_u * Pp_v * P0_w / 36.0

                # XZ-plane edges
                f[QM0M,x,y,z] = rho_init * Pm_u * P0_v * Pm_w / 36.0
                f[QM0P,x,y,z] = rho_init * Pm_u * P0_v * Pp_w / 36.0
                f[QP0M,x,y,z] = rho_init * Pp_u * P0_v * Pm_w / 36.0
                f[QP0P,x,y,z] = rho_init * Pp_u * P0_v * Pp_w / 36.0

                # YZ-plane edges
                f[Q0MM,x,y,z] = rho_init * P0_u * Pm_v * Pm_w / 36.0
                f[Q0MP,x,y,z] = rho_init * P0_u * Pm_v * Pp_w / 36.0
                f[Q0PM,x,y,z] = rho_init * P0_u * Pp_v * Pm_w / 36.0
                f[Q0PP,x,y,z] = rho_init * P0_u * Pp_v * Pp_w / 36.0

            end#x
        end#y
    end#z
   
    ##-------- Initialise fS's --------##
    fS .= f

    # Force garbage collection to free unused memory
    GC.gc()

    ##--------  Logging  --------##
    Log_Simulation_Header()
    Log_Grid_Dimensions(gridlengthX, gridlengthY, gridlengthZ, delta_x, delta_t)
    Log_Fluid_Parameters(τ, omega, Inflow_Velocity, lattice_inflow_velocity, Re_phys, Re_lattice)
    Log_Simulation_Start()
    
    ##-------- Logging into .CSV --------##
    run_tag = "Re$(reynoldsNumber)_Ma$(Mach_Number)_dx$(delta_x)"

    wake_csv_paths = ["simulation_data/wake_profil_$(lbl)_$(run_tag).csv" for lbl in probe_labels]

    mkpath("simulation_data")
    mkpath("visualization")

    for path in wake_csv_paths
        open(path, "w") do io
            println(io, "t_phys, y_phys, u, v, w, mean_u, mean_v, mean_w, std_u, std_v, std_w")
        end
    end

    ##--------  Plot calls  --------## 
    # Initialise Observables
    vx_xy_obs = nothing; step_text_vx_xy = nothing
    vx_xz_obs = nothing; step_text_vx_xz = nothing
    mag_xy_obs = nothing; step_text_mag_xy = nothing
    mag_xz_obs = nothing; step_text_mag_xz = nothing
    vx_xz_front_obs = nothing; step_text_vx_xz_front = nothing
    vx_xz_back_obs = nothing; step_text_vx_xz_back = nothing
    vort_xy_obs = nothing; step_text_vort_xy = nothing
    vort_xz_obs = nothing; step_text_vort_xz = nothing

    
    if Plotvx
        vx_xy_obs, step_text_vx_xy, vx_xz_obs, step_text_vx_xz = 
            setup_vx_plot(gridlengthX, gridlengthY, gridlengthZ, velocityX, midY, midZ)
    end

    if Plotmag
        mag_xy_obs, step_text_mag_xy, mag_xz_obs, step_text_mag_xz = 
            setup_mag_plot(gridlengthX, gridlengthY, gridlengthZ, velocityMag, midY, midZ;
                   colorrange=(0.0, 1.2 * Inflow_Velocity),
                   velocity_scale=delta_x / delta_t)
    end

    if Plotdebug
        vx_xz_front_obs, step_text_vx_xz_front, vx_xz_back_obs, step_text_vx_xz_back =
            setup_debug_plots(gridlengthX, gridlengthZ, velocityX, frontY, backY)
    end
   
    if Plotvorticity
    vort_xy_obs, step_text_vort_xy, vort_xz_obs, step_text_vort_xz =
        setup_vorticity_plot(gridlengthX, gridlengthY, gridlengthZ, vortZ, vortY, midY, midZ; mode=vorticity_mode)
    end
    
    ##--------  MAIN  LOOP --------## 
    # Run Simulation Loop
    for i in 1:simulationTime


        # mnups tracking start + estimation start
        t0 = time_ns()
            
        # with rho
            # collision_stream!(
            #     gridlengthX, gridlengthY, gridlengthZ, τ, CS, is_fluid,
            #     rho, u, v, w,
            #     f, fS, disc_nodes, F_x_lat
            # )        
        # # without rho global u
        #     collision_stream!(
        #         gridlengthX, gridlengthY, gridlengthZ, τ, CS, is_fluid,
        #         u, v, w,
        #         f, fS, disc_nodes, F_x_lat
        #     )        
    
        # local u
        collision_stream!(
            gridlengthX, gridlengthY, gridlengthZ, τ, CS, is_fluid,
            u, v, w,
            f, fS, disc_nodes, C_T_local
        )        
        
        ##-------- free slip walls --------##
        # y-faces (front+back) in one barrier, z-faces (bot+top) in another.
        # y-faces and z-faces stay sequential to avoid corner node conflicts.
        let nf = length(wall_front_x), nb = length(wall_back_x)
            @inbounds Threads.@threads :static for i in 1:(nf + nb)
                if i <= nf
                    x = wall_front_x[i]; z = wall_front_z[i]
                    fS[Q0P0, x, 2, z] = fS[Q0M0, x, 1, z]
                    fS[QMP0, x, 2, z] = fS[QMM0, x, 1, z]
                    fS[QPP0, x, 2, z] = fS[QPM0, x, 1, z]
                    fS[Q0PM, x, 2, z] = fS[Q0MM, x, 1, z]
                    fS[Q0PP, x, 2, z] = fS[Q0MP, x, 1, z]
                else
                    j = i - nf
                    x = wall_back_x[j]; z = wall_back_z[j]
                    fS[Q0M0, x, gridlengthY-1, z] = fS[Q0P0, x, gridlengthY, z]
                    fS[QMM0, x, gridlengthY-1, z] = fS[QMP0, x, gridlengthY, z]
                    fS[QPM0, x, gridlengthY-1, z] = fS[QPP0, x, gridlengthY, z]
                    fS[Q0MM, x, gridlengthY-1, z] = fS[Q0PM, x, gridlengthY, z]
                    fS[Q0MP, x, gridlengthY-1, z] = fS[Q0PP, x, gridlengthY, z]
                end
            end
        end
        let nbot = length(wall_bot_x), ntop = length(wall_top_x)
            @inbounds Threads.@threads :static for i in 1:(nbot + ntop)
                if i <= nbot
                    x = wall_bot_x[i]; y = wall_bot_y[i]
                    fS[Q00P, x, y, 2] = fS[Q00M, x, y, 1]
                    fS[QP0P, x, y, 2] = fS[QP0M, x, y, 1]
                    fS[QM0P, x, y, 2] = fS[QM0M, x, y, 1]
                    fS[Q0PP, x, y, 2] = fS[Q0PM, x, y, 1]
                    fS[Q0MP, x, y, 2] = fS[Q0MM, x, y, 1]
                else
                    j = i - nbot
                    x = wall_top_x[j]; y = wall_top_y[j]
                    fS[Q00M, x, y, gridlengthZ-1] = fS[Q00P, x, y, gridlengthZ]
                    fS[QP0M, x, y, gridlengthZ-1] = fS[QP0P, x, y, gridlengthZ]
                    fS[QM0M, x, y, gridlengthZ-1] = fS[QM0P, x, y, gridlengthZ]
                    fS[Q0PM, x, y, gridlengthZ-1] = fS[Q0PP, x, y, gridlengthZ]
                    fS[Q0MM, x, y, gridlengthZ-1] = fS[Q0MP, x, y, gridlengthZ]
                end
            end
        end

        # Corner handling: free slip corners get noslip bounceback values
        # Corner handling after "normal" freeslip logic to rewrite corner populations
        @inbounds for x in 1:gridlengthX
            # front/bot edge
            if !is_solid[x, 2, 2]
                fS[Q0PP, x, 2, 2] = fS[Q0MM, x, 1, 1]
            end

            # front/top edge
            if !is_solid[x, 2, gridlengthZ-1]
                fS[Q0PM, x, 2, gridlengthZ-1] = fS[Q0MP, x, 1, gridlengthZ]
            end

            # back/bot edge
            if !is_solid[x, gridlengthY-1, 2]
                fS[Q0MP, x, gridlengthY-1, 2] = fS[Q0PM, x, gridlengthY, 1]
            end

            # back/top edge
            if !is_solid[x, gridlengthY-1, gridlengthZ-1]
                fS[Q0MM, x, gridlengthY-1, gridlengthZ-1] = fS[Q0PP, x, gridlengthY, gridlengthZ]
            end
        end
        ##-------- free slip walls  end --------##


        # INLET: moving wall bounceback with momentum addition
        # # compute inflow populations fp00S, fpp0S, fpm0S, fp0pS, fp0mS
        # # momentum coefficients for D3Q19 weights
        # # Compute new populations
        @views fS[QP00, 2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fS[QM00, 1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_face
        @views fS[QPP0, 2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fS[QMP0, 1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        @views fS[QPM0, 2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fS[QMM0, 1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        @views fS[QP0P, 2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fS[QM0P, 1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        @views fS[QP0M, 2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fS[QM0M, 1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        
        # OUTLET: interpolation (Non reflective Geier et al. 2015)
        # f_new(x_b, t) = cs * f(x_{b-1}, t-dt) + (1 - cs) * f(x_b, t-dt)
        @views fS[QM00, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= lattice_speedOfSound * f[QM00, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ (1 - lattice_speedOfSound) * f[QM00, gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        @views fS[QMM0, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= lattice_speedOfSound * f[QMM0, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ (1 - lattice_speedOfSound) * f[QMM0, gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        @views fS[QMP0, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= lattice_speedOfSound * f[QMP0, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ (1 - lattice_speedOfSound) * f[QMP0, gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        @views fS[QM0M, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= lattice_speedOfSound * f[QM0M, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ (1 - lattice_speedOfSound) * f[QM0M, gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        @views fS[QM0P, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= lattice_speedOfSound * f[QM0P, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ (1 - lattice_speedOfSound) * f[QM0P, gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]


        # Swap: SWAP POINTERS new distribution to "old"
        f, fS = fS, f     

        # mnups tracking end
        t_mnups_s = (time_ns() - t0) * 1e-9
        mnups = n_mnups_nodes / (t_mnups_s * 1e6)

        # simulation time estimation
        if i >= 5 && i<= 15
            time_samples[i-4] = t_mnups_s
            println("Step $i mainloop: $(round(t_mnups_s * 1000, digits=1))ms")
        end
        
        if i == 15
            avg_time = sum(time_samples) / 11
            est_total_s = avg_time * simulationTime
            est_hours = floor(Int, est_total_s / 3600)
            est_minutes = floor(Int, (est_total_s % 3600) / 60)
            println("---> Estimated total simulation time: ~$(est_hours)h $(est_minutes)min ($simulationTime) steps x $(round(avg_time*1000, digits=1))ms")    
        end


        ##-------- Probe Sampling (cumulativ mean) --------##
        if i % sample_interval == 0
            buf_ptr         += 1
            cumulativ_count += 1
            sample_times[buf_ptr] = i * delta_t

            fac = cumulativ_count > 1 ? 1.0 / (cumulativ_count -1) : 0.0
            
            for p in 1:n_probes
                px = probe_xs[p]
                @inbounds for j in eachindex(probe_ys)
                    y = probe_ys[j]
                    up = u[px, y , probe_z]
                    vp = v[px, y , probe_z]
                    wp = w[px, y , probe_z]

                    sample_buf_u[p][j, buf_ptr] = up
                    sample_buf_v[p][j, buf_ptr] = vp
                    sample_buf_w[p][j, buf_ptr] = wp

                    du = up - cumulativ_mean_u[p][j]
                    cumulativ_mean_u[p][j] += du / cumulativ_count
                    cumulativ_M2_u[p][j]    += du * (up - cumulativ_mean_u[p][j])

                    dv = vp - cumulativ_mean_v[p][j]
                    cumulativ_mean_v[p][j] += dv / cumulativ_count
                    cumulativ_M2_v[p][j]   += dv * (vp - cumulativ_mean_v[p][j])

                    dw = wp - cumulativ_mean_w[p][j]
                    cumulativ_mean_w[p][j] += dw / cumulativ_count
                    cumulativ_M2_w[p][j]   += dw * (wp - cumulativ_mean_w[p][j])

                    sample_buf_mean_u[p][j, buf_ptr] = cumulativ_mean_u[p][j]
                    sample_buf_mean_v[p][j, buf_ptr] = cumulativ_mean_v[p][j]
                    sample_buf_mean_w[p][j, buf_ptr] = cumulativ_mean_w[p][j]
                    sample_buf_std_u[p][j, buf_ptr]  = sqrt(cumulativ_M2_u[p][j] * fac)
                    sample_buf_std_v[p][j, buf_ptr]  = sqrt(cumulativ_M2_v[p][j] * fac)
                    sample_buf_std_w[p][j, buf_ptr]  = sqrt(cumulativ_M2_w[p][j] * fac)
                end
            end

            if buf_ptr == samples_per_flush
                for p in 1:n_probes
                    open(wake_csv_paths[p], "a") do io
                        for s in 1:samples_per_flush-1
                            for j in eachindex(probe_ys)
                                y_phys = (probe_ys[j]-1) * delta_x
                                println(io, "$(sample_times[s]), $y_phys, $(sample_buf_u[p][j,s]), $(sample_buf_v[p][j,s]), $(sample_buf_w[p][j,s]), $(sample_buf_mean_u[p][j,s]), $(sample_buf_mean_v[p][j,s]), $(sample_buf_mean_w[p][j,s]), $(sample_buf_std_u[p][j,s]), $(sample_buf_std_v[p][j,s]), $(sample_buf_std_w[p][j,s])")
                            end
                        end
                    end
                    for buf in (sample_buf_u[p], sample_buf_v[p], sample_buf_w[p],
                                sample_buf_mean_u[p], sample_buf_mean_v[p], sample_buf_mean_w[p],
                                sample_buf_std_u[p], sample_buf_std_v[p], sample_buf_std_w[p])
                        buf[:,1] .= buf[:, samples_per_flush]
                    end
                end
                sample_times[1] = sample_times[samples_per_flush]
                buf_ptr = 1
            end
        end

        # Logging timestep and mnups in console
        if (i % 100 == 0) || (i == simulationTime)
            Log_Simulation_Runtime(i, simulationTime)
            println("MNUPS: $(round(mnups, digits=2))")

            F_total = 0.0
            for idx in disc_nodes
                u_disc = u[idx]
                F_total += -0.5 * C_T_local * (u_disc * u_disc)
            end
            A_disc = Float64(length(disc_nodes))
            C_T_check = abs(F_total) / (0.5 * A_disc * lattice_inflow_velocity^2)
            println("C_T target: $(round(C_T, digits=4)) | C_T_local: $(round(C_T_local, digits=4)) | C_T_eff: $(round(C_T_check, digits=4)) | ratio C_T_eff/C_T: $(round(C_T_check/C_T, digits=3))")
        
        end

        

        # Plot of the field
        if any((Plotvx, Plotdebug, Plotmag, Plotvorticity)) && ((i % 100 == 0) || (i == simulationTime))

            velocityX .= u 
            # velocityY .= v
            # velocityZ .= w
            @. velocityMag = sqrt(u^2 + v^2 + w^2)

            @inbounds for idx in disc_nodes
                velocityMag[idx] = NaN
                velocityX[idx] = NaN
            end

            # Vorticity calculation
            if Plotvorticity
                if vorticity_mode == :magnitude
                    @inbounds for z in 2:gridlengthZ-1, y in 2:gridlengthY-1, x in 2:gridlengthX-1
                    if is_solid[x,y,z]
                        vortZ[x,y,z] = NaN; vortY[x,y,z] = NaN
                    else
                        wx = (w[x,y+1,z] - w[x,y-1,z]) * 0.5 - (v[x,y,z+1] - v[x,y,z-1]) * 0.5
                        wy = (u[x,y,z+1] - u[x,y,z-1]) * 0.5 - (w[x+1,y,z] - w[x-1,y,z]) * 0.5
                        wz = (v[x+1,y,z] - v[x-1,y,z]) * 0.5 - (u[x,y+1,z] - u[x,y-1,z]) * 0.5
                        mag = sqrt(wx*wx + wy*wy + wz*wz)
                        vortZ[x,y,z] = mag; vortY[x,y,z] = mag
                    end
                end
                else  # :component
                    @inbounds for z in 2:gridlengthZ-1, y in 2:gridlengthY-1, x in 2:gridlengthX-1
                        if is_solid[x,y,z]
                            vortZ[x,y,z] = NaN; vortY[x,y,z] = NaN
                        else
                            vortZ[x,y,z] = (v[x+1,y,z] - v[x-1,y,z]) * 0.5 - (u[x,y+1,z] - u[x,y-1,z]) * 0.5
                            vortY[x,y,z] = (u[x,y,z+1] - u[x,y,z-1]) * 0.5 - (w[x+1,y,z] - w[x-1,y,z]) * 0.5
                        end
                    end
                end
            end
                    
            update_plots!(Plotmag, Plotvx, Plotdebug,
                          velocityMag, velocityX,
                          gridlengthX, gridlengthY, gridlengthZ, midY, midZ, frontY, backY,
                          i, simulationTime, delta_t,
                          mag_xy_obs, step_text_mag_xy, mag_xz_obs, step_text_mag_xz,
                          vx_xy_obs, step_text_vx_xy, vx_xz_obs, step_text_vx_xz,
                          vx_xz_front_obs, step_text_vx_xz_front, vx_xz_back_obs, step_text_vx_xz_back;
                          Plotvorticity=Plotvorticity,
                          vortZ=vortZ, vortY=vortY,
                          vort_xy_obs=vort_xy_obs, step_text_vort_xy=step_text_vort_xy,
                          vort_xz_obs=vort_xz_obs, step_text_vort_xz=step_text_vort_xz,
                          mag_velocity_scale=delta_x/delta_t)

            yield()
            sleep(0.01)
        end#any((Plotvx, Plotdebug, Plotmag)) && ((i % 200 == 0) || (i == simulationTime))

    end#i in 1:simulationTime

    ##-------- Log final step --------##
    for p in 1:n_probes
        open(wake_csv_paths[p], "a") do io
            for s in 1:buf_ptr
                for j in eachindex(probe_ys)
                    y_phys = (probe_ys[j] - 1) * delta_x
                    println(io, join([sample_times[s], y_phys,
                        sample_buf_u[p][j,s], sample_buf_v[p][j,s], sample_buf_w[p][j,s],
                        sample_buf_mean_u[p][j,s], sample_buf_mean_v[p][j,s], sample_buf_mean_w[p][j,s],
                        sample_buf_std_u[p][j,s], sample_buf_std_v[p][j,s], sample_buf_std_w[p][j,s]
                    ], ", "))
                end
            end
        end
    end

    Log_Simulation_Tail()

    # save last plot for post processing
    snapshot_path = "visualization/snapshot_$(run_tag).jld2"
    jldsave(snapshot_path;
        velocityMag,
        velocityX,
        vortY,
        vortZ,
        gridlengthX,
        gridlengthY,
        gridlengthZ,
        midY,
        midZ,
        delta_x,
        delta_t,
        run_tag
    )

end#run_JuLattice()

run_JuLattice()