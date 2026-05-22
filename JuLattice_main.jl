############################
## Main file for JuLattice #
############################
include("src/Plotter.jl")
include("src/Logger.jl")
include("src/BoundaryConditions.jl")
include("src/TurbulenceModel.jl")
include("src/Kernel.jl")


using Serialization # for saving last plot
using MeshGrid, GLMakie
using .Plotter, .Logger
using .BoundaryConditions
using .TurbulenceModel 
using .Kernel

function run_JuLattice()
    ####################################  Initialize  ####################################
    ##-------- User Settings --------##
    # Cylinder Definition
    Radius   = 0.0115 #0.08 #0.1        # m (D = 0.023m)
    D = 2 * Radius
    
    # #Simulation Domain Settings
    # lateral 13D (both sides y&z)v| outflow 15.5D: (experimental setup) NO-SLIP DOMAIN
    # length_X = 20.5 * D                 # m
    # length_Y = 0.6                      # m
    # length_Z = 0.6                      # m

    # # # lateral 5D (both sides y&z) | outflow 10D: FREE-SLIP DOMAIN
    # length_X = 15.5 * D
    # length_Y = 10 * D
    # length_Z = 10 * D

    # # lateral 10D (both sides y&z) | outflow 15D: FREE-SLIP DOMAIN
    length_X = 20.5 * D   # extended: 20.5D 
    length_Y = 15 * D     # extended: 15D  
    length_Z = 15 * D     # extended: 15D  

    # Fluid Settings 
    Kinematic_Viscosity = 1e-6                                       # m^2/s 
    reynoldsNumber =   2760 #2760                                    # Target Reynolds number
    Mach_Number = 0.1 # 0.05                                         # Target Mach number (Ma = U_lattice/c_s)
                                                                     # Keep Ma < 0.1 for incompressible flow!

    # Simulation Settings
    Simulation_Time = 60 #0.5;                                   # s
    
    # Grid spacing (physical units per lattice unit)
    # 0.0023 => 10 = D/Δx || 0.00115 => 20 = D/Δx || 0.00153 => 15 = D/Δx
    delta_x         = 0.00092                                
   
    # Smagorinsky constant CS
    CS              = 1/3 #0.1 #0.17                 # CS ↑ = eddy viscosity ↑

    # Plot Requests (Flags)
    Plotvx = false;
    Plotmag = true;
    Plotdebug = false;
    Plotvorticity = false;
    vorticity_mode = :component # :component (ω_z / ω_y)   or   :magnitude (|ω|)

    ##-------- Compute LBM Parameters from Mach Number --------##
    lattice_speedOfSound    = 1.0 / sqrt(3)
    Inflow_Velocity         = reynoldsNumber * Kinematic_Viscosity / (2 * Radius)
    speedOfSound            = Inflow_Velocity / Mach_Number
    delta_t                 = delta_x * lattice_speedOfSound / speedOfSound
    lattice_viscosity       = Kinematic_Viscosity * delta_t / (delta_x)^2
    lattice_inflow_velocity = Mach_Number * lattice_speedOfSound

    # nu_lattice = c_s² * (tau - 0.5) => tau = nu_lattice / c_s² + 0.5
    τ       = lattice_viscosity / (lattice_speedOfSound * lattice_speedOfSound) + 0.5
    omega   = 1.0 / τ

    fluiddensity = 1.0 # lattice units
    simulationTime = ceil(Int, Simulation_Time / delta_t);  #lattice units
    
    ##-------- Convert user settings to lattice units --------##
    # Domain
    gridlengthX = ceil(Int, length_X / delta_x);
    gridlengthY = ceil(Int, length_Y / delta_x);
    gridlengthZ = ceil(Int, length_Z / delta_x);

    # Cylinder Position
    cylinder_x      = Int(round((5.5 * D) / delta_x)) + 1
    cylinder_y      = Int(round((length_Y/ 2 ) / delta_x)) + 1
    cylinder_z_top  = length_Z  #length_Z * 0.75 #
    cylinder_z_bot  = 0.0 #length_Z * 0.25 #
    cylinder_radius = Radius/delta_x

    # Grid-idx for is_object (nodes inside of cylinder)
    cylinder_start = 2 #1 #2 + Int(floor(cylinder_z_bot / delta_x))
    cylinder_end   = gridlengthZ - 1 #gridlengthZ #gridlengthZ-1 #2 + Int(ceil(cylinder_z_top / delta_x))  
    
    # Reynolds Check:
    # Re_lattice = U*R/v -> should match Re_phys since quantities are scaled
    Re_lattice = floor(Int, ((lattice_inflow_velocity .* 2 .* cylinder_radius)/lattice_viscosity)) 
    Re_phys = Inflow_Velocity * 2 * Radius / Kinematic_Viscosity

    # Define Slice indices for plotting
    midY = 2 + Int(round((gridlengthY-2)/2))
    midZ = 2 + Int(round((gridlengthZ-2)/2))

    ##-------- Probe Setup --------## 
    D_lat   = Int(round(D / delta_x))
    probe_x_3D = cylinder_x + 3 * D_lat
    probe_x_6D = cylinder_x + 6 * D_lat
    probe_z = midZ
    probe_ys = collect(2:gridlengthY-1)
    n_probe = length(probe_ys) # Vector{Int64}

    # sample_dt_phys      = 0.1   # sampling rate = 10Hz
    sample_dt_phys      = 0.01    # sampling rate = 100Hz
    log_dt_phys         = 1.0   # logging rate for csv-flush
    sample_interval     = max(1, round(Int, sample_dt_phys / delta_t))
    samples_per_flush   = max(1, round(Int, log_dt_phys / sample_dt_phys))

    sample_times        = zeros(samples_per_flush)
    cumulativ_count     = 0
    buf_ptr = 0

    # 3D behind cylinider
    sample_buf_u_3D        = zeros(n_probe, samples_per_flush)
    sample_buf_v_3D        = zeros(n_probe, samples_per_flush)
    sample_buf_w_3D        = zeros(n_probe, samples_per_flush)
    sample_buf_mean_u_3D   = zeros(n_probe, samples_per_flush)
    sample_buf_mean_v_3D   = zeros(n_probe, samples_per_flush)
    sample_buf_mean_w_3D   = zeros(n_probe, samples_per_flush)
    sample_buf_std_u_3D    = zeros(n_probe, samples_per_flush)
    sample_buf_std_v_3D    = zeros(n_probe, samples_per_flush)
    sample_buf_std_w_3D    = zeros(n_probe, samples_per_flush)

    cumulativ_mean_u_3D    = zeros(n_probe); cumulativ_M2_u_3D = zeros(n_probe)
    cumulativ_mean_v_3D    = zeros(n_probe); cumulativ_M2_v_3D = zeros(n_probe)
    cumulativ_mean_w_3D    = zeros(n_probe); cumulativ_M2_w_3D = zeros(n_probe)

    # 6D behind cylinder
    sample_buf_u_6D        = zeros(n_probe, samples_per_flush)
    sample_buf_v_6D        = zeros(n_probe, samples_per_flush)
    sample_buf_w_6D        = zeros(n_probe, samples_per_flush)
    sample_buf_mean_u_6D   = zeros(n_probe, samples_per_flush)
    sample_buf_mean_v_6D   = zeros(n_probe, samples_per_flush)
    sample_buf_mean_w_6D   = zeros(n_probe, samples_per_flush)
    sample_buf_std_u_6D    = zeros(n_probe, samples_per_flush)
    sample_buf_std_v_6D    = zeros(n_probe, samples_per_flush)
    sample_buf_std_w_6D    = zeros(n_probe, samples_per_flush)

    cumulativ_mean_u_6D    = zeros(n_probe); cumulativ_M2_u_6D = zeros(n_probe)
    cumulativ_mean_v_6D    = zeros(n_probe); cumulativ_M2_v_6D = zeros(n_probe)
    cumulativ_mean_w_6D    = zeros(n_probe); cumulativ_M2_w_6D = zeros(n_probe)


    # more slices for debugg plots
    frontY = 2
    backY = gridlengthY-1
    botZ = 2
    topZ = gridlengthZ-1
    nearFrontY = 10
    nearBackY = gridlengthY-10
    nearBotZ = 10
    nearTopZ = gridlengthZ-10

    # Inlet momentum coefficients for D3Q19 weights (velocity bounceback)
    inlet_add_face = (2.0 / (18.0 * lattice_speedOfSound^2)) * lattice_inflow_velocity
    inlet_add_edge = (2.0 / (36.0 * lattice_speedOfSound^2)) * lattice_inflow_velocity

    ##--------  classify nodes --------##
    ## create solid node mask
    # Array{Bool} instead of BitArray: single byte load in kernel loop vs bit-unpack
    is_solid  = fill(false, gridlengthX, gridlengthY, gridlengthZ)
    is_object = fill(false, gridlengthX, gridlengthY, gridlengthZ)
    is_fluid  = fill(false, gridlengthX, gridlengthY, gridlengthZ)

    # pre compute fluid range
    is_fluid[2:gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= true

    # Solid and object mask
    for x in 1:gridlengthX, y in 1:gridlengthY, z in 1:gridlengthZ
        # walls
        if y==1 || y==gridlengthY || z==1 || z==gridlengthZ
            is_solid[x, y, z] = true
            continue
        end

        # cylinder vertically (y-axis)
        dx = x - cylinder_x
        dy = y - cylinder_y
        if (z >= cylinder_start) && (z <= cylinder_end) && (sqrt(dx^2 + dy^2) <= cylinder_radius)
            is_object[x, y, z] = true
            is_solid[x, y, z] = true
        end
    end
    # object_indices = findall(is_object)

    # fluid mask
    is_fluid .&= .!is_object
    n_fluid_nodes = sum(is_fluid)
    n_cylinder_nodes = sum(is_object)
    n_mnups_nodes = n_fluid_nodes + n_cylinder_nodes

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
    

    ##-------- precompute BC --------##
    println("Computing Bouzidi boundary data...")

    boundary_data = compute_object_boundary_data(
        gridlengthX, gridlengthY, gridlengthZ,
        (cylinder_x -2) * delta_x,  (cylinder_y -2) * delta_x,
        cylinder_radius * delta_x,
        cylinder_z_bot, cylinder_z_top,
        is_object, delta_x
    )
    A_lat = (2.0 * cylinder_radius) * Float64(cylinder_end - cylinder_start + 1)
    Cd = 0.0
    Cl = 0.0

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
    rho        = Array{Float64}(undef, gridlengthX, gridlengthY, gridlengthZ)
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
                rho[x,y,z]        = 1.0
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
    D_over_dx = Int(round(D / delta_x))
    run_tag = "Re$(reynoldsNumber)_Ma$(Mach_Number)_DdeltaX$(D_over_dx)"

    wake_csv_path_3D = "simulation_data/wake_profil_3D_$(run_tag).csv"
    wake_csv_path_6D = "simulation_data/wake_profil_6D_$(run_tag).csv"

    forces_csv_path = "simulation_data/forces_$(run_tag).csv"
    mkpath("simulation_data")
    mkpath("visualization")

    # wakevelocities, mean and std
    open(wake_csv_path_3D, "w") do io
        println(io, "t_phys, y_phys, u, v, w, mean_u, mean_v, mean_w, std_u, std_v, std_w")
    end
    open(wake_csv_path_6D, "w") do io
        println(io, "t_phys, y_phys, u, v, w, mean_u, mean_v, mean_w, std_u, std_v, std_w")
    end

    # forces
    forces_io = open(forces_csv_path, "w")
    println(forces_io, "t_phys, Cd, Cl")


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
            setup_mag_plot(gridlengthX, gridlengthY, gridlengthZ, velocityMag, midY, midZ)
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

        #t_estimation = @elapsed begin

        # mnups tracking start + estimation start
        t0 = time_ns()
    
            collision_stream!(
                gridlengthX, gridlengthY, gridlengthZ, τ, CS, is_fluid,
                rho, u, v, w,
                f, fS
            )
        
        #end #end elapsed

        t_estimation = (time_ns() - t0) * 1e-9
        # debug timecheck for mainloop with elapsed
        if i >= 5 && i<= 15
            println("Step $i mainloop: $(round(t_estimation * 1000, digits=1))ms")
        end
        
        if i == 15
            est_total_s = t_estimation * simulationTime
            est_hours = floor(Int, est_total_s / 3600)
            est_minutes = floor(Int, (est_total_s % 3600) / 60)
            println("---> Estimated total simulation time: ~$(est_hours)h $(est_minutes)min ($simulationTime) steps x $(round(t_estimation*1000, digits=1))ms")    
        end
      
        # # bounce-back object | Bouzidi bounceback (IBB)
        # F_x_lat, F_y_lat = apply_bouzidi_bc_3d!(boundary_data,
        #                      fm00S, fp00S, f0m0S, f0p0S, f00mS, f00pS,
        #                      fmm0S, fmp0S, fpm0S, fpp0S,
        #                      fm0mS, fm0pS, fp0mS, fp0pS,
        #                      f0mmS, f0mpS, f0pmS, f0ppS)
        
        # Cd = 2.0 * F_x_lat / (fluiddensity * lattice_inflow_velocity^2 * A_lat)
        # Cl = 2.0 * F_y_lat / (fluiddensity * lattice_inflow_velocity^2 * A_lat)




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


        # # bounce-back walls
        # @inbounds for wall in wall_indices
        #     x, y, z = Tuple(wall)
            
        #     # +x 
        #     if x+1 <= gridlengthX && !is_solid[x+1, y, z]
        #         fp00S[x+1, y, z] = fm00S[x, y, z]
        #     end
        #     # -x 
        #     if x-1 >= 1 && !is_solid[x-1, y, z]
        #         fm00S[x-1, y, z] = fp00S[x, y, z]
        #     end
        #     # +y 
        #     if y+1 <= gridlengthY && !is_solid[x, y+1, z]
        #         f0p0S[x, y+1, z] = f0m0S[x, y, z]
        #     end
        #     # -y 
        #     if y-1 >= 1 && !is_solid[x, y-1, z]
        #         f0m0S[x, y-1, z] = f0p0S[x, y, z]
        #     end
        #     # +z 
        #     if z+1 <= gridlengthZ && !is_solid[x, y, z+1]
        #         f00pS[x, y, z+1] = f00mS[x, y, z]
        #     end
        #     # -z 
        #     if z-1 >= 1 && !is_solid[x, y, z-1]
        #         f00mS[x, y, z-1] = f00pS[x, y, z]
        #     end

        #     # XY
        #     if x+1 <= gridlengthX && y+1 <= gridlengthY && !is_solid[x+1, y+1, z]
        #         fpp0S[x+1, y+1, z] = fmm0S[x, y, z]
        #     end
        #     if x-1 >= 1 && y-1 >= 1 && !is_solid[x-1, y-1, z]
        #         fmm0S[x-1, y-1, z] = fpp0S[x, y, z]
        #     end
        #     if x+1 <= gridlengthX && y-1 >= 1 && !is_solid[x+1, y-1, z]
        #         fpm0S[x+1, y-1, z] = fmp0S[x, y, z]
        #     end
        #     if x-1 >= 1 && y+1 <= gridlengthY && !is_solid[x-1, y+1, z]
        #         fmp0S[x-1, y+1, z] = fpm0S[x, y, z]
        #     end

        #     # XZ
        #     if x+1 <= gridlengthX && z+1 <= gridlengthZ && !is_solid[x+1, y, z+1]
        #         fp0pS[x+1, y, z+1] = fm0mS[x, y, z]
        #     end
        #     if x-1 >= 1 && z-1 >= 1 && !is_solid[x-1, y, z-1]
        #         fm0mS[x-1, y, z-1] = fp0pS[x, y, z]
        #     end
        #     if x+1 <= gridlengthX && z-1 >= 1 && !is_solid[x+1, y, z-1]
        #         fp0mS[x+1, y, z-1] = fm0pS[x, y, z]
        #     end
        #     if x-1 >= 1 && z+1 <= gridlengthZ && !is_solid[x-1, y, z+1]
        #         fm0pS[x-1, y, z+1] = fp0mS[x, y, z]
        #     end

        #     # YZ
        #     if y+1 <= gridlengthY && z+1 <= gridlengthZ && !is_solid[x, y+1, z+1]
        #         f0ppS[x, y+1, z+1] = f0mmS[x, y, z]
        #     end
        #     if y-1 >= 1 && z-1 >= 1 && !is_solid[x, y-1, z-1]
        #         f0mmS[x, y-1, z-1] = f0ppS[x, y, z]
        #     end
        #     if y+1 <= gridlengthY && z-1 >= 1 && !is_solid[x, y+1, z-1]
        #         f0pmS[x, y+1, z-1] = f0mpS[x, y, z]
        #     end
        #     if y-1 >= 1 && z+1 <= gridlengthZ && !is_solid[x, y-1, z+1]
        #         f0mpS[x, y-1, z+1] = f0pmS[x, y, z]
        #     end
            
        # end#for wall in wall_indices



        # bounce-back object | Bouzidi bounceback (IBB)
        F_x_lat, F_y_lat = apply_bouzidi_bc_3d!(boundary_data, fS)
        
        Cd = 2.0 * F_x_lat / (fluiddensity * lattice_inflow_velocity^2 * A_lat)
        Cl = 2.0 * F_y_lat / (fluiddensity * lattice_inflow_velocity^2 * A_lat)


        # INLET: moving wall bounceback with momentum addition
        # # compute inflow populations fp00S, fpp0S, fpm0S, fp0pS, fp0mS
        # # momentum coefficients for D3Q19 weights
        # # Compute new populations
        @views fS[QP00, 2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fS[QM00, 1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_face
        @views fS[QPP0, 2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fS[QMP0, 1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        @views fS[QPM0, 2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fS[QMM0, 1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        @views fS[QP0P, 2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fS[QM0P, 1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        @views fS[QP0M, 2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fS[QM0M, 1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        
        # # OUTLET: no-gradient bounceback 
        # # # all populations that stream in -x direction from previous neighbor
        # # # fm00S, fmm0S, fmp0S, fm0mS, fm0pS
        
        # @views fm00S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm00S[gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        # @views fmm0S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fmm0S[gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        # @views fmp0S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fmp0S[gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        # @views fm0mS[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm0mS[gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        # @views fm0pS[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm0pS[gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        
        
        # OUTLET: interpolation (Non reflective Geier et al. 2015)
        # f_new(x_b, t) = cs * f(x_{b-1}, t-dt) + (1 - cs) * f(x_b, t-dt)
        @views fS[QM00, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= lattice_speedOfSound * f[QM00, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ (1 - lattice_speedOfSound) * f[QM00, gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        @views fS[QMM0, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= lattice_speedOfSound * f[QMM0, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ (1 - lattice_speedOfSound) * f[QMM0, gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        @views fS[QMP0, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= lattice_speedOfSound * f[QMP0, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ (1 - lattice_speedOfSound) * f[QMP0, gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        @views fS[QM0M, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= lattice_speedOfSound * f[QM0M, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ (1 - lattice_speedOfSound) * f[QM0M, gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        @views fS[QM0P, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= lattice_speedOfSound * f[QM0P, gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ (1 - lattice_speedOfSound) * f[QM0P, gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]

        # mnups tracking end
        t_mnups_s = (time_ns() - t0) * 1e-9
        mnups = n_mnups_nodes / (t_mnups_s * 1e6)


        # Swap: SWAP POINTERS new distribution to "old"
        f, fS = fS, f     

        ##-------- Probe Sampling (cumulativ mean) --------##
        if i % sample_interval == 0
            buf_ptr         += 1
            cumulativ_count += 1
            sample_times[buf_ptr] = i * delta_t

            println(forces_io, "$(i * delta_t), $Cd, $Cl")
            flush(forces_io)

            fac = cumulativ_count > 1 ? 1.0 / (cumulativ_count -1) : 0.0
            
            @inbounds for j in eachindex(probe_ys)
                y = probe_ys[j]

                # 3D probe
                u3 = u[probe_x_3D, y, probe_z]
                v3 = v[probe_x_3D, y, probe_z]
                w3 = w[probe_x_3D, y, probe_z]

                sample_buf_u_3D[j, buf_ptr] = u3
                sample_buf_v_3D[j, buf_ptr] = v3
                sample_buf_w_3D[j, buf_ptr] = w3

                du = u3 - cumulativ_mean_u_3D[j];
                cumulativ_mean_u_3D[j] += du / cumulativ_count
                cumulativ_M2_u_3D[j] += du * (u3 - cumulativ_mean_u_3D[j])
                
                dv = v3 - cumulativ_mean_v_3D[j];
                cumulativ_mean_v_3D[j] += dv / cumulativ_count
                cumulativ_M2_v_3D[j] += dv * (v3 - cumulativ_mean_v_3D[j])

                dw = w3 - cumulativ_mean_w_3D[j];
                cumulativ_mean_w_3D[j] += dw / cumulativ_count
                cumulativ_M2_w_3D[j] += dw * (w3 - cumulativ_mean_w_3D[j])

                sample_buf_mean_u_3D[j, buf_ptr] = cumulativ_mean_u_3D[j]
                sample_buf_mean_v_3D[j, buf_ptr] = cumulativ_mean_v_3D[j]
                sample_buf_mean_w_3D[j, buf_ptr] = cumulativ_mean_w_3D[j]
                sample_buf_std_u_3D[j, buf_ptr] = sqrt(cumulativ_M2_u_3D[j] * fac)
                sample_buf_std_v_3D[j, buf_ptr] = sqrt(cumulativ_M2_v_3D[j] * fac)
                sample_buf_std_w_3D[j, buf_ptr] = sqrt(cumulativ_M2_w_3D[j] * fac)

                # 6D probe
                u6 = u[probe_x_6D, y, probe_z]                
                v6 = v[probe_x_6D, y, probe_z]
                w6 = w[probe_x_6D, y, probe_z]

                sample_buf_u_6D[j, buf_ptr] = u6
                sample_buf_v_6D[j, buf_ptr] = v6
                sample_buf_w_6D[j, buf_ptr] = w6


                du = u6 - cumulativ_mean_u_6D[j];
                cumulativ_mean_u_6D[j] += du / cumulativ_count
                cumulativ_M2_u_6D[j] += du * (u6 - cumulativ_mean_u_6D[j])
                
                dv = v6 - cumulativ_mean_v_6D[j];
                cumulativ_mean_v_6D[j] += dv / cumulativ_count
                cumulativ_M2_v_6D[j] += dv * (v6 - cumulativ_mean_v_6D[j])

                dw = w6 - cumulativ_mean_w_6D[j];
                cumulativ_mean_w_6D[j] += dw / cumulativ_count
                cumulativ_M2_w_6D[j] += dw * (w6 - cumulativ_mean_w_6D[j])

                sample_buf_mean_u_6D[j, buf_ptr] = cumulativ_mean_u_6D[j]
                sample_buf_mean_v_6D[j, buf_ptr] = cumulativ_mean_v_6D[j]
                sample_buf_mean_w_6D[j, buf_ptr] = cumulativ_mean_w_6D[j]
                sample_buf_std_u_6D[j, buf_ptr] = sqrt(cumulativ_M2_u_6D[j] * fac)
                sample_buf_std_v_6D[j, buf_ptr] = sqrt(cumulativ_M2_v_6D[j] * fac)
                sample_buf_std_w_6D[j, buf_ptr] = sqrt(cumulativ_M2_w_6D[j] * fac)

            end


            if buf_ptr ==  samples_per_flush
                # 3D flush
                open(wake_csv_path_3D, "a") do io
                    for s in 1:samples_per_flush-1
                        for j in eachindex(probe_ys)
                            y_phys = (probe_ys[j] - 1) * delta_x
                            println(io, "$(sample_times[s]), $y_phys, $(sample_buf_u_3D[j,s]), $(sample_buf_v_3D[j,s]), $(sample_buf_w_3D[j,s]), $(sample_buf_mean_u_3D[j,s]), $(sample_buf_mean_v_3D[j,s]), $(sample_buf_mean_w_3D[j,s]), $(sample_buf_std_u_3D[j,s]), $(sample_buf_std_v_3D[j,s]), $(sample_buf_std_w_3D[j,s])")
                        end
                    end
                end
                
                for buf in (sample_buf_u_3D, sample_buf_v_3D, sample_buf_w_3D,
                            sample_buf_mean_u_3D, sample_buf_mean_v_3D, sample_buf_mean_w_3D,
                            sample_buf_std_u_3D, sample_buf_std_v_3D, sample_buf_std_w_3D)
                    buf[:,1] .= buf[:, samples_per_flush]
                end

                # 6D flush
                open(wake_csv_path_6D, "a") do io
                    for s in 1:samples_per_flush-1
                        for j in eachindex(probe_ys)
                            y_phys = (probe_ys[j] - 1) * delta_x
                            println(io, "$(sample_times[s]), $y_phys, $(sample_buf_u_6D[j,s]), $(sample_buf_v_6D[j,s]), $(sample_buf_w_6D[j,s]), $(sample_buf_mean_u_6D[j,s]), $(sample_buf_mean_v_6D[j,s]), $(sample_buf_mean_w_6D[j,s]), $(sample_buf_std_u_6D[j,s]), $(sample_buf_std_v_6D[j,s]), $(sample_buf_std_w_6D[j,s])")
                        end
                    end
                end

                for buf in (sample_buf_u_6D, sample_buf_v_6D, sample_buf_w_6D,
                            sample_buf_mean_u_6D, sample_buf_mean_v_6D, sample_buf_mean_w_6D,
                            sample_buf_std_u_6D, sample_buf_std_v_6D, sample_buf_std_w_6D)
                    buf[:,1] .= buf[:, samples_per_flush]
                end

                sample_times[1] = sample_times[samples_per_flush]
                buf_ptr = 1
            end

        end

        # Logging timestep and mnups in console
        if (i % 100 == 0) || (i == simulationTime)
            Log_Simulation_Runtime(i, simulationTime)
            println("MNUPS: $(round(mnups, digits=2))")
        end

        # Plot of the field
        if any((Plotvx, Plotdebug, Plotmag, Plotvorticity)) && ((i % 100 == 0) || (i == simulationTime))

            velocityX .= u 
            # velocityY .= v
            # velocityZ .= w
            @. velocityMag = sqrt(u^2 + v^2 + w^2)
            velocityX[is_object] .= NaN
            velocityMag[is_object] .= NaN

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
                          vort_xz_obs=vort_xz_obs, step_text_vort_xz=step_text_vort_xz)

            yield()
            sleep(0.01)
        end#any((Plotvx, Plotdebug, Plotmag)) && ((i % 200 == 0) || (i == simulationTime))

    end#i in 1:simulationTime

    ##-------- Log final step --------##
    # 3D
    open(wake_csv_path_3D, "a") do io
        for s in 1:buf_ptr
            for j in eachindex(probe_ys)
                y_phys = (probe_ys[j] - 1) * delta_x

                println(io, join([sample_times[s], y_phys,
                        sample_buf_u_3D[j,s], sample_buf_v_3D[j,s], sample_buf_w_3D[j,s],
                        sample_buf_mean_u_3D[j,s], sample_buf_mean_v_3D[j,s], sample_buf_mean_w_3D[j,s],
                        sample_buf_std_u_3D[j,s], sample_buf_std_v_3D[j,s], sample_buf_std_w_3D[j,s]
                ], ", "))
            end
        end
    end

    # 6D
    open(wake_csv_path_6D, "a") do io
        for s in 1:buf_ptr
            for j in eachindex(probe_ys)
                y_phys = (probe_ys[j] - 1) * delta_x

                println(io, join([sample_times[s], y_phys,
                sample_buf_u_6D[j,s], sample_buf_v_6D[j,s], sample_buf_w_6D[j,s],
                sample_buf_mean_u_6D[j,s], sample_buf_mean_v_6D[j,s], sample_buf_mean_w_6D[j,s],
                sample_buf_std_u_6D[j,s], sample_buf_std_v_6D[j,s], sample_buf_std_w_6D[j,s],
                ], ", "))
            end
        end
    end


    # close forces.csv
    close(forces_io)

    Log_Simulation_Tail()

    # save last plot for post processing
    snapshot_path = "visualization/snapshot_$(run_tag).jls"
    serialize(snapshot_path, (
        velocityMag     = copy(velocityMag),
        velocityX       = copy(velocityX),
        vortY           = copy(vortY),
        vortZ           = copy(vortZ),
        is_object       = copy(is_object),
        gridlengthX     = gridlengthX,
        gridlengthY     = gridlengthY,
        gridlengthZ     = gridlengthZ,
        midY            = midY,
        midZ            = midZ,
        delta_t         = delta_t,
        run_tag         = run_tag
    ))

end#run_JuLattice()
# run_JuLattice()