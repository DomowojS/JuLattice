############################
## Main file for JuLattice #
############################
include("src/Plotter.jl")
include("src/Logger.jl")
include("src/BoundaryConditions.jl")

using MeshGrid, GLMakie
using .Plotter, .Logger
using .BoundaryConditions

function run_JuLattice()
    ####################################  Initialize  ####################################
    ##-------- User Settings --------##
    # Domain Settings
    length_X = 2              # m
    length_Y = 0.51           # m 
    length_Z = 1              # m

    # Cylinder Definition
    Radius   = 0.1    # m
    # Position = [length_X/4, length_Y/2, length_Z/2] # m [x,y,z]

    # Fluid Settings 
    Fluid_Density = 1000.0; #1000.0;         # kg/m^3
    # Inflow_Velocity = 1.0 #0.4;              # m/s
    Kinematic_Viscosity = 0.0004; #0.001;    # m^2/s 

    # Reynolds and Mach Number Input
    reynoldsNumber = 500                        # Target Reynolds number
    Mach_Number = 0.01;                         # Target Mach number (Ma = U_lattice/c_s)
                                                # Keep Ma < 0.1 for incompressible flow!

    # Simulation Settings
    Simulation_Time = 600;  #8000               # s
    delta_x = 0.01 #0.01;                       # Grid spacing (physical units per lattice unit)
   
    # Plot Requests (Flags)
    Plotvx = true;
    Plotvy = false;
    Plotvz = false;
    Plotvorticity = false;
    Plotdebug = false;


    ####-------- Run Simulation --------#####
    Log_Simulation_Header()

    ##-------- Compute LBM Parameters from Mach Number --------##
    lattice_speedOfSound = 1.0 / sqrt(3)

    Inflow_Velocity = reynoldsNumber * Kinematic_Viscosity / (2 * Radius)

    speedOfSound = Inflow_Velocity / Mach_Number
    # speedOfSound = Kinematic_Viscosity * reynoldsNumber / (Mach_Number * 2 * Radius)
    delta_t = delta_x * lattice_speedOfSound / speedOfSound

    lattice_viscosity = Kinematic_Viscosity * delta_t / (delta_x)^2

    lattice_inflow_velocity = Mach_Number * lattice_speedOfSound

    # Step 4: Relaxation time and omega from lattice viscosity
    # nu_lattice = c_s² * (tau - 0.5) => tau = nu_lattice / c_s² + 0.5
    τ = lattice_viscosity / (lattice_speedOfSound * lattice_speedOfSound) + 0.5
    omega  = 1.0 / τ

    ##-------- Convert user settings to lattice units --------##
    # Domain
    gridlengthX = ceil(Int, length_X / delta_x);
    gridlengthY = ceil(Int, length_Y / delta_x);
    gridlengthZ = ceil(Int, length_Z / delta_x);

    println("nodes in x: $gridlengthX")
    println("nodes in x: $gridlengthY")
    println("nodes in x: $gridlengthZ")


    # Cyliner
    cylinder_x = Int(round(gridlengthX / 3))
    cylinder_y = 2 + Int(round((gridlengthY-2)/2))
    println("cylinder_y = $cylinder_y")
    #cylinder_y = Int(round(gridlengthY / 2))
    cylinder_radius = Radius/delta_x
    cylinder_start = Int(round(gridlengthZ*0.25))
    cylinder_end = Int(round(gridlengthZ*0.75))

    # Fluid
    fluiddensity = Fluid_Density

    # ReynoldsCheck
    lattice_Re = (lattice_inflow_velocity .* 2 .* cylinder_radius)/lattice_viscosity; #Re_lattice = U*R/v -> sollte Re entsprechen weil Größen skaliert wurden
    lattice_Re_Log=floor(Int,lattice_Re)

    # Log 
    Log_Discretization_Settings(delta_x, delta_t, lattice_Re_Log)

    # Print τ value
    println("##### Computed Fluid Values #####")
    println("Computed relaxation time τ  = ", round(τ, digits=10))
    println("Computed omega ω = ", round(omega, digits=10))
    println("Computed Inflow_Velocity (phsical) = ", round(Inflow_Velocity, digits=10))
    println("Reynolds number check:")
    Re_phys = Inflow_Velocity * 2 * Radius / Kinematic_Viscosity
    println("   Re (physical) = ", round(Re_phys, digits=2))
    println("   Re (lattice) = ", round(lattice_Re, digits=2))
    println("#################################")

    # Simulation Settings
    simulationTime = ceil(Int, Simulation_Time / delta_t);

    #Q = 19; #D3Q19
    # D3Q19
    # f000 = rest (0,0,0)
    # fm00, fp00 = x-axis (±1,0,0)
    # f0m0, f0p0 = y-axis (0,±1,0)
    # f00m, f00p = z-axis (0,0,±1)
    # fmm0, fmp0, fpm0, fpp0 = xy-plane edges
    # fm0m, fm0p, fp0m, fp0p = xz-plane edges
    # f0mm, f0mp, f0pm, f0pp = yz-plane edges

    #Define arrays for each direction D3Q9
    # f000 = rest (0,0,0)
    f000 = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # fm00, fp00 = x-axis (±1,0,0)
    fm00 = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fp00 = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # f0m0, f0p0 = y-axis (0,±1,0)
    f0m0 = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    f0p0 = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # f00m, f00p = z-axis (0,0,±1)
    f00m = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    f00p = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # fmm0, fmp0, fpm0, fpp0 = xy-plane edges
    fmm0 = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fmp0 = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fpm0 = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fpp0 = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # fm0m, fm0p, fp0m, fp0p = xz-plane edges
    fm0m = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fm0p = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fp0m = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fp0p = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # f0mm, f0mp, f0pm, f0pp = yz-plane edges
    f0mm = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    f0mp = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    f0pm = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    f0pp = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)

    #Define array for each direction after Collision+stream (S)
    # f000 = rest (0,0,0)
    f000S = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # fm00, fp00 = x-axis (±1,0,0)
    fm00S = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fp00S = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # f0m0, f0p0 = y-axis (0,±1,0)
    f0m0S = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    f0p0S = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # f00m, f00p = z-axis (0,0,±1)
    f00mS = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    f00pS = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # fmm0, fmp0, fpm0, fpp0 = xy-plane edges
    fmm0S = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fmp0S = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fpm0S = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fpp0S = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # fm0m, fm0p, fp0m, fp0p = xz-plane edges
    fm0mS = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fm0pS = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fp0mS = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    fp0pS = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # f0mm, f0mp, f0pm, f0pp = yz-plane edges
    f0mmS = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    f0mpS = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    f0pmS = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    f0ppS = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)

    #Initialise macroscopic variables
    rho = ones(Float32, gridlengthX, gridlengthY, gridlengthZ) .* fluiddensity
    u = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)    #ux
    v = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)    #uy
    w = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)    #uz

    # create grid
    gridX, gridY, gridZ = meshgrid(1:gridlengthX, 1:gridlengthY, 1:gridlengthZ);

    #Swap of Y and X axis: (Y,X,Z) -> (X,Y,Z)
    if size(gridX) == (gridlengthY, gridlengthX, gridlengthZ) #check for format of grids
        gridX = permutedims(gridX, (2,1,3))
        gridY = permutedims(gridY, (2,1,3))
        gridZ = permutedims(gridZ, (2,1,3))
    end    

    # Initialise velocity arrays for plotting
    velocityX = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    velocityY = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    velocityZ = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)

    # omegaX = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # omegaY = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # omegaZ = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)
    # omegaMag = zeros(Float32, gridlengthX, gridlengthY, gridlengthZ)

    #Define mid-Planes for plotting
    midY = 2 + Int(round((gridlengthY-2)/2))
    midZ = 2 + Int(round((gridlengthZ-2)/2))
    println("midY = $midY")
    println("midZ = $midZ")

    # more slices for debugging
    frontY = 2
    backY = gridlengthY-1
    botZ = 2
    topZ = gridlengthZ-1

    nearFrontY = 10
    nearBackY = gridlengthY-10
    nearBotZ = 10
    nearTopZ = gridlengthZ-10


    ## create solid node mask
    is_solid = falses(gridlengthX, gridlengthY, gridlengthZ)
    is_wall = falses(gridlengthX, gridlengthY, gridlengthZ)
    is_object = falses(gridlengthX, gridlengthY, gridlengthZ)

    for x in 1:gridlengthX, y in 1:gridlengthY, z in 1:gridlengthZ
        # walls
        if y==1 || y==gridlengthY || z==1 || z==gridlengthZ
            is_wall[x, y, z] = true
            is_solid[x, y, z] = true
            continue
        end
        
        # cylinder vertically (y-axis)

        dx = x- cylinder_x
        dy = y - cylinder_y
        if (z >= cylinder_start) && (z <= cylinder_end) && (sqrt(dx^2 + dy^2) <= cylinder_radius)
            is_object[x, y, z] = true
            is_solid[x, y, z] = true
        end

    end
    wall_indices = findall(is_wall)
    object_indices = findall(is_object)

    println("Computing Bouzidi boundary data...")
    boundary_data = compute_object_boundary_data(
        gridlengthX, gridlengthY, gridlengthZ,
        cylinder_x * delta_x,  cylinder_y * delta_x,
        cylinder_radius * delta_x,
        cylinder_start, cylinder_end,
        is_object, delta_x
    )

    ## Initialize distribution functions FLUID NODES and SOLID NODES
    for x in 1:gridlengthX
        for y in 1:gridlengthY
            for z in 1:gridlengthZ

                ux = is_solid[x, y, z] ? 0.0 : lattice_inflow_velocity
                #ux = lattice_inflow_velocity
                uy = 0.0
                uz = 0.0
                rho_init = fluiddensity
                
                # Pre-compute polynomial factors
                ux2 = ux * ux
                uy2 = uy * uy
                uz2 = uz * uz
                
                Pm_u = 1 - 3*ux + 3*ux2
                #P0_u = -2 + 3*ux2
                P0_u = 1 - 1.5*ux2
                Pp_u = 1 + 3*ux + 3*ux2
                
                Pm_v = 1 - 3*uy + 3*uy2
                #P0_v = -2 + 3*uy2
                P0_v = 1 - 1.5*uy2
                Pp_v = 1 + 3*uy + 3*uy2
                
                Pm_w = 1 - 3*uz + 3*uz2
                #P0_w = -2 + 3*uz2
                P0_w = 1 - 1.5*uz2
                Pp_w = 1 + 3*uz + 3*uz2
        
                # Push scheme: Rest particle (0,0,0) - weight 1/3
                f000[x,y,z] = rho_init * P0_u * P0_v * P0_w / 3.0
                
                # Push scheme: Face neighbors - weight 1/18
                fm00[x,y,z] = rho_init * Pm_u * P0_v * P0_w / 18.0 
                fp00[x,y,z] = rho_init * Pp_u * P0_v * P0_w / 18.0 
                
                f0m0[x,y,z] = rho_init * P0_u * Pm_v * P0_w / 18.0 
                f0p0[x,y,z] = rho_init * P0_u * Pp_v * P0_w / 18.0
                
                f00m[x,y,z] = rho_init * P0_u * P0_v * Pm_w / 18.0 
                f00p[x,y,z] = rho_init * P0_u * P0_v * Pp_w / 18.0 
                
                # Push scheme: Edge neighbors - weight 1/36
                # XY-plane edges
                fmm0[x,y,z] = rho_init * Pm_u * Pm_v * P0_w / 36.0 
                fmp0[x,y,z] = rho_init * Pm_u * Pp_v * P0_w / 36.0 
                fpm0[x,y,z] = rho_init * Pp_u * Pm_v * P0_w / 36.0 
                fpp0[x,y,z] = rho_init * Pp_u * Pp_v * P0_w / 36.0 
                
                # XZ-plane edges
                fm0m[x,y,z] = rho_init * Pm_u * P0_v * Pm_w / 36.0 
                fm0p[x,y,z] = rho_init * Pm_u * P0_v * Pp_w / 36.0
                fp0m[x,y,z] = rho_init * Pp_u * P0_v * Pm_w / 36.0 
                fp0p[x,y,z] = rho_init * Pp_u * P0_v * Pp_w / 36.0 
                
                # YZ-plane edges
                f0mm[x,y,z] = rho_init * P0_u * Pm_v * Pm_w / 36.0 
                f0mp[x,y,z] = rho_init * P0_u * Pm_v * Pp_w / 36.0 
                f0pm[x,y,z] = rho_init * P0_u * Pp_v * Pm_w / 36.0 
                f0pp[x,y,z] = rho_init * P0_u * Pp_v * Pp_w / 36.0
            
            end
        end
    end

    #Initialise fS-Arrays for the first time as f-Arrays
    f000S .= f000 
    fm00S .= fm00
    fp00S .= fp00
    f0m0S .= f0m0 
    f0p0S .= f0p0 
    f00mS .= f00m 
    f00pS .= f00p 
    fmm0S .= fmm0 
    fmp0S .= fmp0 
    fpm0S .= fpm0 
    fpp0S .= fpp0 
    fm0mS .= fm0m 
    fm0pS .= fm0p 
    fp0mS .= fp0m 
    fp0pS .= fp0p 
    f0mmS .= f0mm 
    f0mpS .= f0mp 
    f0pmS .= f0pm
    f0ppS .= f0pp

    ##  Plot calls
    if any((Plotvorticity, Plotvx, Plotvy, Plotvz, Plotdebug))
                    
        if Plotvx==true
            vx_xy_obs, step_text_vx_xy, fig_vx_xy = Create_Plot_XY(gridlengthX-2, gridlengthY-2, velocityX[2:gridlengthX-1,2:gridlengthY-1,midZ]; title="v_x at z=$(midZ)")
            screen_vx_xy = GLMakie.Screen()
            display(screen_vx_xy, fig_vx_xy)

            vx_xz_obs, step_text_vx_xz, fig_vx_xz = Create_Plot_XZ(gridlengthX-2, gridlengthZ-2, velocityX[2:gridlengthX-1,midY,2:gridlengthZ-1]; title="v_x at y=$(midY)")
            screen_vx_xz = GLMakie.Screen()
            display(screen_vx_xz, fig_vx_xz)
        end

        if Plotdebug == true
            # vx | xz at front wall (y=2)
            vx_xz_front_obs, step_text_vx_xz_front, fig_vx_xz_front = Create_Plot_XZ(gridlengthX, gridlengthZ, velocityX[:,frontY,:]; title="v_x at y=$(frontY) [FRONT]")
            screen_vx_xz_front = GLMakie.Screen()
            display(screen_vx_xz_front, fig_vx_xz_front)

            # # vx | xz at near front (y=10)
            # vx_xz_nearfront_obs, step_text_vx_xz_nearfront, fig_vx_xz_nearfront = Create_Plot_XZ(gridlengthX, gridlengthZ, velocityX[:,nearFrontY,:]; title="v_x at y=$(nearFrontY) [near FRONT]")
            # screen_vx_xz_nearfront = GLMakie.Screen()
            # display(screen_vx_xz_nearfront, fig_vx_xz_nearfront)

            # vx | xz at back wall (y=gridlengthY-1)
            vx_xz_back_obs, step_text_vx_xz_back, fig_vx_xz_back = Create_Plot_XZ(gridlengthX, gridlengthZ, velocityX[:,backY,:]; title="v_x at y=$(backY) [BACK]")
            screen_vx_xz_back = GLMakie.Screen()
            display(screen_vx_xz_back, fig_vx_xz_back)

            # # vx | xz at near back (y=gridlengthY-10)
            # vx_xz_nearback_obs, step_text_vx_xz_nearback, fig_vx_xz_nearback = Create_Plot_XZ(gridlengthX, gridlengthZ, velocityX[:,nearBackY,:]; title="v_x at y=$(nearBackY) [near BACK]")
            # screen_vx_xz_nearback = GLMakie.Screen()
            # display(screen_vx_xz_nearback, fig_vx_xz_nearback)

            # # vy | xz at front wall (y=2)
            # vy_xz_front_obs, step_text_vy_xz_front, fig_vy_xz_front = Create_Plot_XY(gridlengthX, gridlengthZ, velocityY[:,frontY,:]; title="v_y at y=$(frontY) [FRONT]")
            # screen_vy_xz_front = GLMakie.Screen()
            # display(screen_vy_xz_front, fig_vy_xz_front)

            # # vy | xz at back wall (y=gridlengthY-1)
            # vy_xz_back_obs, step_text_vy_xz_back, fig_vy_xz_back = Create_Plot_XY(gridlengthX, gridlengthZ, velocityY[:,backY,:]; title="v_y at y=$(backY) [BACK]")
            # screen_vy_xz_back = GLMakie.Screen()
            # display(screen_vy_xz_back, fig_vy_xz_back)
            
        end

        if Plotvy==true 
            velocityY_obs, text_obj_vy, step_text_vy, fig_vy = Create_Plot(gridlengthX, gridlengthY, velocityY, "Y")
            screen3 = GLMakie.Screen()
            display(screen3, fig_vy)
        end

        if Plotvz==true
            velocityZ_obs, text_obj_vz, step_text_vz, fig_vz = Create_Plot(gridlengthX, gridlengthZ, velocityZ, "Z")
            screen4 = GLMakie.Screen()
            display(screen4, fig_vz)
        end

        if Plotvorticity == true
            #xy slice at z=midZ
            omega_xy_obs, step_text_omega_xy, fig_omega_xy = Create_Vorticity_XY(gridlengthX, gridlengthY, omegaMag[:,:,midZ]; title = "|ω| at z=$(midZ)")
            screen_omega_xy = GLMakie.Screen()
            display(screen_omega_xy, fig_omega_xy)

            #xz slice at y=midY
            omega_xz_obs, step_text_omega_xz, fig_omega_xz = Create_Vorticity_XZ(gridlengthX, gridlengthZ, omegaMag[:,midY,:]; title = "|ω| at y=$(midY)")
            screen_omega_xz = GLMakie.Screen()
            display(screen_omega_xz, fig_omega_xz)    
        end

    end #Plot calls

    println("#################################")
    println("Starting Simulation:")
    # Run Simulation Loop
    for i in 1:simulationTime

        # for x in 2:gridlengthX-1 #Iteration über alle Zellen außer die Randzellen
        #     for y in 2:gridlengthY-1
        #         for z in 2:gridlengthZ-1

        for z in 2:gridlengthZ-1 #Iteration über alle Zellen außer die Randzellen
            for y in 2:gridlengthY-1
                for x in 2:gridlengthX-1


                    if is_solid[x,y,z]
                        continue
                    end
                
                    # Compute macroscopic quantities
                    rho[x,y,z] = f000[x,y,z] + 
                                (fm00[x,y,z] + fp00[x,y,z] + f0m0[x,y,z] + f0p0[x,y,z] + f00m[x,y,z] + f00p[x,y,z]) +
                                (fmm0[x,y,z] + fmp0[x,y,z] + fpm0[x,y,z] + fpp0[x,y,z] + 
                                fm0m[x,y,z] + fm0p[x,y,z] + fp0m[x,y,z] + fp0p[x,y,z] +
                                f0mm[x,y,z] + f0mp[x,y,z] + f0pm[x,y,z] + f0pp[x,y,z])
                    
                    u[x,y,z] = ((-fm00[x,y,z] + fp00[x,y,z]) +
                                (-fmm0[x,y,z] - fmp0[x,y,z] + fpm0[x,y,z] + fpp0[x,y,z]) +
                                (-fm0m[x,y,z] - fm0p[x,y,z] + fp0m[x,y,z] + fp0p[x,y,z])) / rho[x,y,z]
                    
                    v[x,y,z] = ((-f0m0[x,y,z] + f0p0[x,y,z]) +
                                (-fmm0[x,y,z] + fmp0[x,y,z] - fpm0[x,y,z] + fpp0[x,y,z]) +
                                (-f0mm[x,y,z] - f0mp[x,y,z] + f0pm[x,y,z] + f0pp[x,y,z])) / rho[x,y,z]
                    
                    w[x,y,z] = ((-f00m[x,y,z] + f00p[x,y,z]) +
                                (-fm0m[x,y,z] + fm0p[x,y,z] - fp0m[x,y,z] + fp0p[x,y,z]) +
                                (-f0mm[x,y,z] + f0mp[x,y,z] - f0pm[x,y,z] + f0pp[x,y,z])) / rho[x,y,z]
                    
                    # Pre-compute polynomial factors
                    u2 = u[x,y,z] * u[x,y,z]
                    v2 = v[x,y,z] * v[x,y,z]
                    w2 = w[x,y,z] * w[x,y,z]
                    
                    Pm_u = 1 - 3*u[x,y,z] + 3*u2
                    P0_u = 1 - 1.5*u2
                    Pp_u = 1 + 3*u[x,y,z] + 3*u2
                    
                    Pm_v = 1 - 3*v[x,y,z] + 3*v2
                    P0_v = 1 - 1.5*v2
                    Pp_v = 1 + 3*v[x,y,z] + 3*v2
                    
                    Pm_w = 1 - 3*w[x,y,z] + 3*w2
                    P0_w = 1 - 1.5*w2
                    Pp_w = 1 + 3*w[x,y,z] + 3*w2
                    
                    # Push scheme: Rest particle (0,0,0) - weight 1/3
                    f000S[x,y,z] = f000[x,y,z] + omega * (rho[x,y,z] * P0_u * P0_v * P0_w / 3.0 - f000[x,y,z])
                    
                    # Push scheme: Face neighbors - weight 1/18
                    fm00S[x-1,y,z] = fm00[x,y,z] + omega * (rho[x,y,z] * Pm_u * P0_v * P0_w / 18.0 - fm00[x,y,z])
                    fp00S[x+1,y,z] = fp00[x,y,z] + omega * (rho[x,y,z] * Pp_u * P0_v * P0_w / 18.0 - fp00[x,y,z])
                    
                    f0m0S[x,y-1,z] = f0m0[x,y,z] + omega * (rho[x,y,z] * P0_u * Pm_v * P0_w / 18.0 - f0m0[x,y,z])
                    f0p0S[x,y+1,z] = f0p0[x,y,z] + omega * (rho[x,y,z] * P0_u * Pp_v * P0_w / 18.0 - f0p0[x,y,z])
                    
                    f00mS[x,y,z-1] = f00m[x,y,z] + omega * (rho[x,y,z] * P0_u * P0_v * Pm_w / 18.0 - f00m[x,y,z])
                    f00pS[x,y,z+1] = f00p[x,y,z] + omega * (rho[x,y,z] * P0_u * P0_v * Pp_w / 18.0 - f00p[x,y,z])
                    
                    # Push scheme: Edge neighbors - weight 1/36
                    # XY-plane edges
                    fmm0S[x-1,y-1,z] = fmm0[x,y,z] + omega * (rho[x,y,z] * Pm_u * Pm_v * P0_w / 36.0 - fmm0[x,y,z])
                    fmp0S[x-1,y+1,z] = fmp0[x,y,z] + omega * (rho[x,y,z] * Pm_u * Pp_v * P0_w / 36.0 - fmp0[x,y,z])
                    fpm0S[x+1,y-1,z] = fpm0[x,y,z] + omega * (rho[x,y,z] * Pp_u * Pm_v * P0_w / 36.0 - fpm0[x,y,z])
                    fpp0S[x+1,y+1,z] = fpp0[x,y,z] + omega * (rho[x,y,z] * Pp_u * Pp_v * P0_w / 36.0 - fpp0[x,y,z])
                    
                    # XZ-plane edges
                    fm0mS[x-1,y,z-1] = fm0m[x,y,z] + omega * (rho[x,y,z] * Pm_u * P0_v * Pm_w / 36.0 - fm0m[x,y,z])
                    fm0pS[x-1,y,z+1] = fm0p[x,y,z] + omega * (rho[x,y,z] * Pm_u * P0_v * Pp_w / 36.0 - fm0p[x,y,z])
                    fp0mS[x+1,y,z-1] = fp0m[x,y,z] + omega * (rho[x,y,z] * Pp_u * P0_v * Pm_w / 36.0 - fp0m[x,y,z])
                    fp0pS[x+1,y,z+1] = fp0p[x,y,z] + omega * (rho[x,y,z] * Pp_u * P0_v * Pp_w / 36.0 - fp0p[x,y,z])
                    
                    # YZ-plane edges
                    f0mmS[x,y-1,z-1] = f0mm[x,y,z] + omega * (rho[x,y,z] * P0_u * Pm_v * Pm_w / 36.0 - f0mm[x,y,z])
                    f0mpS[x,y-1,z+1] = f0mp[x,y,z] + omega * (rho[x,y,z] * P0_u * Pm_v * Pp_w / 36.0 - f0mp[x,y,z])
                    f0pmS[x,y+1,z-1] = f0pm[x,y,z] + omega * (rho[x,y,z] * P0_u * Pp_v * Pm_w / 36.0 - f0pm[x,y,z])
                    f0ppS[x,y+1,z+1] = f0pp[x,y,z] + omega * (rho[x,y,z] * P0_u * Pp_v * Pp_w / 36.0 - f0pp[x,y,z])

                end
            end
        end

        # bounce-back walls
        for wall in wall_indices
            x, y, z = Tuple(wall)
            # +x 
            if x+1 <= gridlengthX && !is_solid[x+1, y, z]
                fp00S[x+1, y, z] = fm00S[x, y, z]
            end
            # -x 
            if x-1 >= 1 && !is_solid[x-1, y, z]
                fm00S[x-1, y, z] = fp00S[x, y, z]
            end
            # +y 
            if y+1 <= gridlengthY && !is_solid[x, y+1, z]
                f0p0S[x, y+1, z] = f0m0S[x, y, z]
            end
            # -y 
            if y-1 >= 1 && !is_solid[x, y-1, z]
                f0m0S[x, y-1, z] = f0p0S[x, y, z]
            end
            # +z 
            if z+1 <= gridlengthZ && !is_solid[x, y, z+1]
                f00pS[x, y, z+1] = f00mS[x, y, z]
            end
            # -z 
            if z-1 >= 1 && !is_solid[x, y, z-1]
                f00mS[x, y, z-1] = f00pS[x, y, z]
            end

            # XY
            if x+1 <= gridlengthX && y+1 <= gridlengthY && !is_solid[x+1, y+1, z]
                fpp0S[x+1, y+1, z] = fmm0S[x, y, z]
            end
            if x-1 >= 1 && y-1 >= 1 && !is_solid[x-1, y-1, z]
                fmm0S[x-1, y-1, z] = fpp0S[x, y, z]
            end
            if x+1 <= gridlengthX && y-1 >= 1 && !is_solid[x+1, y-1, z]
                fpm0S[x+1, y-1, z] = fmp0S[x, y, z]
            end
            if x-1 >= 1 && y+1 <= gridlengthY && !is_solid[x-1, y+1, z]
                fmp0S[x-1, y+1, z] = fpm0S[x, y, z]
            end

            # XZ
            if x+1 <= gridlengthX && z+1 <= gridlengthZ && !is_solid[x+1, y, z+1]
                fp0pS[x+1, y, z+1] = fm0mS[x, y, z]
            end
            if x-1 >= 1 && z-1 >= 1 && !is_solid[x-1, y, z-1]
                fm0mS[x-1, y, z-1] = fp0pS[x, y, z]
            end
            if x+1 <= gridlengthX && z-1 >= 1 && !is_solid[x+1, y, z-1]
                fp0mS[x+1, y, z-1] = fm0pS[x, y, z]
            end
            if x-1 >= 1 && z+1 <= gridlengthZ && !is_solid[x-1, y, z+1]
                fm0pS[x-1, y, z+1] = fp0mS[x, y, z]
            end

            # YZ
            if y+1 <= gridlengthY && z+1 <= gridlengthZ && !is_solid[x, y+1, z+1]
                f0ppS[x, y+1, z+1] = f0mmS[x, y, z]
            end
            if y-1 >= 1 && z-1 >= 1 && !is_solid[x, y-1, z-1]
                f0mmS[x, y-1, z-1] = f0ppS[x, y, z]
            end
            if y+1 <= gridlengthY && z-1 >= 1 && !is_solid[x, y+1, z-1]
                f0pmS[x, y+1, z-1] = f0mpS[x, y, z]
            end
            if y-1 >= 1 && z+1 <= gridlengthZ && !is_solid[x, y-1, z+1]
                f0mpS[x, y-1, z+1] = f0pmS[x, y, z]
            end
            
        end

        # bounce-back object | Bouzidi bounceback (IBB)
        apply_bouzidi_bc_3d!(boundary_data,
                             fm00S, fp00S, f0m0S, f0p0S, f00mS, f00pS,
                             fmm0S, fmp0S, fpm0S, fpp0S,
                             fm0mS, fm0pS, fp0mS, fp0pS,
                             f0mmS, f0mpS, f0pmS, f0ppS)


        # INLET: moving wall bounceback with momentum addition
        # compute inflow populations fp00S, fpp0S, fpm0S, fp0pS, fp0mS
        # momentum coefficients for D3Q19 weights
        inlet_add_face = (2.0 / (18.0 * lattice_speedOfSound^2)) * lattice_inflow_velocity
        inlet_add_edge = (2.0 / (36.0 * lattice_speedOfSound^2)) * lattice_inflow_velocity

        # Compute new populations
        fp00S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fp00S[1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_face
        fpp0S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fpp0S[1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        fpm0S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fpm0S[1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        fp0pS[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fp0pS[1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        fp0mS[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fp0mS[1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        
        # OUTLET: no-gradient bounceback 
        # all populations that stream in x- direction from previous neighbor
        # fm00S, fmm0S, fmp0S, fm0mS, fm0pS
        
        fm00S[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm00S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1]
        fmm0S[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1] .= fmm0S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1]
        fmp0S[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1] .= fmp0S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1]
        fm0mS[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm0mS[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1]
        fm0pS[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm0pS[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1]
                
        ## Inlet / Outlet BC
        # Inlet Zou-He Velocity
        # unknown: fp00, fpp0, fpm0, fp0p, fp0m (+x direction)
        x_inlet = 1
        ux_inlet = lattice_inflow_velocity

        # Swap: SWAP POINTERS new distribution to "old"
        f000, f000S = f000S, f000
        fm00, fm00S = fm00S, fm00
        fp00, fp00S = fp00S, fp00
        f0m0, f0m0S = f0m0S, f0m0
        f0p0, f0p0S = f0p0S, f0p0
        f00m, f00mS = f00mS, f00m
        f00p, f00pS = f00pS, f00p
        fmm0, fmm0S = fmm0S, fmm0
        fmp0, fmp0S = fmp0S, fmp0
        fpm0, fpm0S = fpm0S, fpm0
        fpp0, fpp0S = fpp0S, fpp0
        fm0m, fm0mS = fm0mS, fm0m
        fm0p, fm0pS = fm0pS, fm0p
        fp0m, fp0mS = fp0mS, fp0m
        fp0p, fp0pS = fp0pS, fp0p
        f0mm, f0mmS = f0mmS, f0mm
        f0mp, f0mpS = f0mpS, f0mp
        f0pm, f0pmS = f0pmS, f0pm
        f0pp, f0ppS = f0ppS, f0pp

        # Plot of the field
        # if ((i % 10 == 0)) || (i == simulationTime)
        # if ((i % 20 == 0)) || (i == simulationTime)
        if ((i % 200 == 0)) || (i == simulationTime)


            #Copy velocities for plotting
            velocityX .= u
            velocityY .= v
            velocityZ .= w

            # Set velocities inside the sphere to zero
            velocityX[is_solid] .= NaN
            velocityY[is_solid] .= NaN
            velocityZ[is_solid] .= NaN

            # # # Compute vorticity 3D
            # #vorticity
            # dv_dx = (circshift(velocityY, (-1,0,0)) .- circshift(velocityY, (1,0,0))) ./ 2
            # dw_dx = (circshift(velocityZ, (-1,0,0)) .- circshift(velocityZ, (1,0,0))) ./ 2

            # du_dy = (circshift(velocityX, (0,-1,0)) .- circshift(velocityX, (0,1,0))) ./ 2
            # dw_dy = (circshift(velocityZ, (0,-1,0)) .- circshift(velocityZ, (0,1,0))) ./ 2
            
            # du_dz = (circshift(velocityX, (0,0,-1)) .- circshift(velocityX, (0,0,1))) ./ 2
            # dv_dz = (circshift(velocityY, (0,0,-1)) .- circshift(velocityY, (0,0,1))) ./ 2

            # omegaX .= dw_dy .- dv_dz
            # omegaY .= du_dz .- dw_dx
            # omegaZ .= dv_dx .-du_dy
            # omegaMag .= sqrt.(omegaX.^2 .+ omegaY.^2 .+ omegaZ.^2)
            
            if ((i % 100 == 0)) || (i == simulationTime)
                Log_Simulation_Runtime(i, simulationTime)
            end
            # Update the observables
            if Plotvorticity == true   
                                    
                #Mask inlet outlet
                omegaMag[inlet] .= 0.0
                omegaMag[outlet] .= 0.0
                # Mask the cylinder region
                omegaMag[sphere] .= NaN
                
                # #3D update: 
                # omegaMag_obs[] = Float32.(omegaMag)
                # step_text_omega[] = "Time step: $i, $(floor(Int, i*delta_t))s"
                
                #2D update:
                #xy at midZ
                omega_xy_obs[] = copy(omegaMag[:,:,midZ])
                step_text_omega_xy[] = "Time step: $i, $(floor(Int, i*delta_t))s" 
                #xz at midY
                omega_xz_obs[] = copy(omegaMag[:,midY,:])
                step_text_omega_xz[] = "Time step: $i, $(floor(Int, i*delta_t))s"

            end

            if Plotvx==true
                #refresh xy slice at z=midZ
                vx_xy_obs[] = copy(velocityX[2:gridlengthX-1, 2:gridlengthY-1, midZ])
                step_text_vx_xy[] = "Time step:$i, $(floor(Int, i*delta_t))s"

                #refresh xz slice at y=midY
                vx_xz_obs[] = copy(velocityX[2:gridlengthX-1 ,midY, 2:gridlengthZ-1])
                step_text_vx_xz[] = "Time step:$i, $(floor(Int, i*delta_t))s"
            end

            if Plotvy==true 
                velocityY_obs[] = copy(velocityY) 
                step_text_vy[] = "Time step: $i, $(floor(Int, i*delta_t))s"
            end

            if Plotdebug==true
                #vx_xz_front_obs[] = copy(velocityX[:,frontY,:])
                vx_xz_front_obs[] = velocityX[:,frontY,:]

                step_text_vx_xz_front[] = "Time step:$i, $(floor(Int, i*delta_t))s"
                
                # vx_xz_nearfront_obs[] = copy(velocityX[:,nearFrontY,:])
                # step_text_vx_xz_nearfront[] = "Time step:$i, $(floor(Int, i*delta_t))s"

                #vx_xz_back_obs[] = copy(velocityX[:,backY,:])
                vx_xz_back_obs[] = velocityX[:,backY,:]

                step_text_vx_xz_back[] = "Time step:$i, $(floor(Int, i*delta_t))s"
                
                # vx_xz_nearback_obs[] = copy(velocityX[:,nearBackY,:])
                # step_text_vx_xz_nearback[] = "Time step:$i, $(floor(Int, i*delta_t))s"

                # vy_xz_front_obs[] = copy(velocityY[:,frontY,:])
                # step_text_vy_xz_front[] = "Time step:$i, $(floor(Int, i*delta_t))s"

                # vy_xz_back_obs[] = copy(velocityY[:,backY,:])
                # step_text_vy_xz_back[] = "Time step:$i, $(floor(Int, i*delta_t))s"

            end
            yield()
            #sleep(0.05)
            sleep(0.05)
        end

    end

    Log_Simulation_Tail()
end
run_JuLattice()