############################
## Main file for JuLattice #
############################
include("src/Plotter.jl")
include("src/Logger.jl")
include("src/BoundaryConditions.jl")
include("src/TurbulenceModel.jl")

using MeshGrid, GLMakie
using .Plotter, .Logger
using .BoundaryConditions
using .TurbulenceModel 

function run_JuLattice()
    ####################################  Initialize  ####################################
    ##-------- User Settings --------##
    # Cylinder Definition
    Radius   = 0.08 #0.1    # m
    D = 2 * Radius
    # Position = [length_X/4, length_Y/2, length_Z/2] # m [x,y,z]
    
    # Simulation Domain Settings
    # length_X = 21 * D              # m
    # length_Y = 11 * D              # m 
    # length_Z = 16 * D              # m

    length_X = 11 * D              # m
    length_Y = 6 * D              # m 
    length_Z = 8 * D              # m



    # Fluid Settings 
    Fluid_Density = 1000.0; #1000.0;         # kg/m^3
    # Inflow_Velocity = 1.0 #0.4;              # m/s
    Kinematic_Viscosity = 0.0004; #0.001;    # m^2/s 

    # Reynolds and Mach Number Input
    reynoldsNumber = 3000 #280 #200 #500                        # Target Reynolds number
    Mach_Number = 0.03 #0.01;                         # Target Mach number (Ma = U_lattice/c_s)
                                                # Keep Ma < 0.1 for incompressible flow!

    # Simulation Settings
    Simulation_Time = 600;  #8000               # s
    delta_x = 0.015 #0.01 ;                       # Grid spacing (physical units per lattice unit)
    CS = 0.1    # CS ↑ = eddy viscosity ↑

    # Plot Requests (Flags)
    Plotvx = false;
    Plotvy = false;
    Plotvz = false;
    Plotvorticity = false;
    Plotdebug = false;
    Plotmag = true;

    # DEBUG Stability Check
    EnableStabilityCheck = false
    CheckEvery = 200
    VelLimit = 0.2
    RhoMinLimit = 100.0
    RhoMaxLimit = 5000.0


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
    println("nodes in y: $gridlengthY")
    println("nodes in xz $gridlengthZ")


    ## Cylinder Position
    cylinder_x = Int(round((5.5 * D) / delta_x)) + 1
    cylinder_y = Int(round((length_Y/ 2 ) / delta_x)) + 1
    cylinder_z_top = length_Z * 0.75
    cylinder_z_bot = length_Z * 0.25

    # convert to lattice units
    cylinder_radius = Radius/delta_x

    # Grid-idx for is_object (nodes inside of cylinder)
    cylinder_start = 2 + Int(floor(cylinder_z_bot / delta_x))
    cylinder_end   = 2 + Int(ceil(cylinder_z_top / delta_x))

    # Fluid
    # fluiddensity = Fluid_Density
    fluiddensity = 1.0 # lattice units
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
    println("Computed lattice_inflow_velocity (lattice) = ", round(lattice_inflow_velocity, digits= 10))
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

    # Define arrays for each direction D3Q9
    # f000 = rest (0,0,0)
    f000 = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fm00, fp00 = x-axis (±1,0,0)
    fm00 = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fp00 = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # f0m0, f0p0 = y-axis (0,±1,0)
    f0m0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
    f0p0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # f00m, f00p = z-axis (0,0,±1)
    f00m = zeros(gridlengthX, gridlengthY, gridlengthZ)
    f00p = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fmm0, fmp0, fpm0, fpp0 = xy-plane edges
    fmm0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fmp0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fpm0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fpp0 = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fm0m, fm0p, fp0m, fp0p = xz-plane edges
    fm0m = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fm0p = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fp0m = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fp0p = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # f0mm, f0mp, f0pm, f0pp = yz-plane edges
    f0mm = zeros(gridlengthX, gridlengthY, gridlengthZ)
    f0mp = zeros(gridlengthX, gridlengthY, gridlengthZ)
    f0pm = zeros(gridlengthX, gridlengthY, gridlengthZ)
    f0pp = zeros(gridlengthX, gridlengthY, gridlengthZ)
    
    # # Define array for each direction after Collision (C)
    # # f000 = rest (0,0,0)
    # f000C = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # # fm00, fp00 = x-axis (±1,0,0)
    # fm00C = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fp00C = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # # f0m0, f0p0 = y-axis (0,±1,0)
    # f0m0C = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # f0p0C = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # # f00m, f00p = z-axis (0,0,±1)
    # f00mC = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # f00pC = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # # fmm0, fmp0, fpm0, fpp0 = xy-plane edges
    # fmm0C = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fmp0C = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fpm0C = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fpp0C = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # # fm0m, fm0p, fp0m, fp0p = xz-plane edges
    # fm0mC = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fm0pC = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fp0mC = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fp0pC = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # # f0mm, f0mp, f0pm, f0pp = yz-plane edges
    # f0mmC = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # f0mpC = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # f0pmC = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # f0ppC = zeros(gridlengthX, gridlengthY, gridlengthZ)

    # Define array for each direction after stream (S)
    # f000 = rest (0,0,0)
    f000S = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fm00, fp00 = x-axis (±1,0,0)
    fm00S = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fp00S = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # f0m0, f0p0 = y-axis (0,±1,0)
    f0m0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
    f0p0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # f00m, f00p = z-axis (0,0,±1)
    f00mS = zeros(gridlengthX, gridlengthY, gridlengthZ)
    f00pS = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fmm0, fmp0, fpm0, fpp0 = xy-plane edges
    fmm0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fmp0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fpm0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fpp0S = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # fm0m, fm0p, fp0m, fp0p = xz-plane edges
    fm0mS = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fm0pS = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fp0mS = zeros(gridlengthX, gridlengthY, gridlengthZ)
    fp0pS = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # f0mm, f0mp, f0pm, f0pp = yz-plane edges
    f0mmS = zeros(gridlengthX, gridlengthY, gridlengthZ)
    f0mpS = zeros(gridlengthX, gridlengthY, gridlengthZ)
    f0pmS = zeros(gridlengthX, gridlengthY, gridlengthZ)
    f0ppS = zeros(gridlengthX, gridlengthY, gridlengthZ)



    #Initialise macroscopic variables
    rho = ones(gridlengthX, gridlengthY, gridlengthZ) .* fluiddensity
    u = zeros(gridlengthX, gridlengthY, gridlengthZ)    #ux
    v = zeros(gridlengthX, gridlengthY, gridlengthZ)    #uy
    w = zeros(gridlengthX, gridlengthY, gridlengthZ)    #uz

    # create grid
    gridX, gridY, gridZ = meshgrid(1:gridlengthX, 1:gridlengthY, 1:gridlengthZ);

    #Swap of Y and X axis: (Y,X,Z) -> (X,Y,Z)
    if size(gridX) == (gridlengthY, gridlengthX, gridlengthZ) #check for format of grids
        gridX = permutedims(gridX, (2,1,3))
        gridY = permutedims(gridY, (2,1,3))
        gridZ = permutedims(gridZ, (2,1,3))
    end    

    # Initialise velocity arrays for plotting
    velocityX = zeros(gridlengthX, gridlengthY, gridlengthZ)
    velocityY = zeros(gridlengthX, gridlengthY, gridlengthZ)
    velocityZ = zeros(gridlengthX, gridlengthY, gridlengthZ)
    velocityMag = zeros(gridlengthX, gridlengthY, gridlengthZ)

    # omegaX = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # omegaY = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # omegaZ = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # omegaMag = zeros(gridlengthX, gridlengthY, gridlengthZ)

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

    ####### ####### classify nodes ####### #######
    ## create solid node mask
    is_solid = falses(gridlengthX, gridlengthY, gridlengthZ)
    is_wall = falses(gridlengthX, gridlengthY, gridlengthZ)
    is_object = falses(gridlengthX, gridlengthY, gridlengthZ)
    is_fluid = falses(gridlengthX, gridlengthY, gridlengthZ)
    
    # pre compute fluid range
    is_fluid[2:gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= true

    # Solid and object mask
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

    # fluid mask
    is_fluid .&= .!is_object
    fluid_indices = findall(is_fluid)

    println("Computing Bouzidi boundary data...")
    boundary_data = compute_object_boundary_data(
        gridlengthX, gridlengthY, gridlengthZ,
        (cylinder_x -2) * delta_x,  (cylinder_y -2) * delta_x,
        cylinder_radius * delta_x,
        cylinder_start, cylinder_end,
        is_object, delta_x
    )


    # momentum coefficients for D3Q19 weights
    inlet_add_face = (2.0 / (18.0 * lattice_speedOfSound^2)) * lattice_inflow_velocity
    inlet_add_edge = (2.0 / (36.0 * lattice_speedOfSound^2)) * lattice_inflow_velocity

    #######  #######    Initialize distribution functions FLUID NODES and SOLID NODES ####### #######
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
   
    # fS's
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

    # Force GARBAGE COLLECTION to free unused memory
    GC.gc()

    ####### #######   Plot calls    ####### ####### 
    if any((Plotvorticity, Plotvx, Plotvy, Plotvz, Plotdebug, Plotmag))
                    
        if Plotvx == true
            fig = Figure(size= (1600, 1200))
            ax1 = Axis(fig[1,1])
            ax2 = Axis(fig[2,1])
            
            vx_xy_obs, step_text_vx_xy, hm1 = Create_Plot_XY(gridlengthX-2, gridlengthY-2, 
                                                                velocityX[2:gridlengthX-1,2:gridlengthY-1,midZ]; 
                                                                title="v_x at z=$(midZ)", ax=ax1)

            vx_xz_obs, step_text_vx_xz, hm2 = Create_Plot_XZ(gridlengthX-2, gridlengthZ-2, 
                                                                velocityX[2:gridlengthX-1,midY,2:gridlengthZ-1]; 
                                                                title="v_x at y=$(midY)", ax=ax2)
            Colorbar(fig[1, 2], hm1, label = "Lattice Velocity")  
            Colorbar(fig[2, 2], hm2, label = "Lattice Velocity")
            display(fig)
        end

        if Plotmag == true
            fig_mag = Figure(size= (1000, 1000))
            ax1_mag = Axis(fig_mag[1,1])
            ax2_mag = Axis(fig_mag[2,1])
            
            mag_xy_obs, step_text_mag_xy, hm1_mag = Create_Plot_Mag_XY(gridlengthX-2, gridlengthY-2, 
                                                                velocityMag[2:gridlengthX-1,2:gridlengthY-1, midZ]; 
                                                                title="|v| at z=$(midZ)", ax=ax1_mag)

            mag_xz_obs, step_text_mag_xz, hm2_mag = Create_Plot_Mag_XZ(gridlengthX-2, gridlengthZ-2, 
                                                                velocityMag[2:gridlengthX-1,midY ,2:gridlengthZ-1]; 
                                                                title="|v| at y=$(midY)", ax=ax2_mag)
            Colorbar(fig_mag[1, 2], hm1_mag, label = "Lattice Velocity Magnitude")  
            Colorbar(fig_mag[2, 2], hm2_mag, label = "Lattice Velocity Magnitude")

            label_mag_xy = Label(fig_mag[3,1:2], text=step_text_mag_xy)
            label_mag_xz = Label(fig_mag[4,1:2], text=step_text_mag_xz)

            rowsize!(fig_mag.layout, 1, Fixed(350))
            rowsize!(fig_mag.layout, 2, Fixed(350))
            rowsize!(fig_mag.layout, 3, Fixed(30))
            rowsize!(fig_mag.layout, 4, Fixed(30))
            colsize!(fig_mag.layout, 1, Auto(0.9))
            colsize!(fig_mag.layout, 2, Auto(0.1))
            rowgap!(fig_mag.layout, 0)


            display(fig_mag)
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

    end ##### Plot calls

    ###### ######   MAIN  LOOP ###### ###### 
    println("#################################")
    println("Starting Simulation:")
    # Run Simulation Loop
    for i in 1:simulationTime

        t_debug = @elapsed begin
        
        # @inbounds Threads.@threads for idx in fluid_indices
            
        #     x, y, z = Tuple(idx)
        # DEBUG THREADING
        let f000=f000, f000S=f000S,
            fm00=fm00, fm00S=fm00S, fp00=fp00, fp00S=fp00S,
            f0m0=f0m0, f0m0S=f0m0S, f0p0=f0p0, f0p0S=f0p0S,
            f00m=f00m, f00mS=f00mS, f00p=f00p, f00pS=f00pS,
            fmm0=fmm0, fmm0S=fmm0S, fmp0=fmp0, fmp0S=fmp0S,
            fpm0=fpm0, fpm0S=fpm0S, fpp0=fpp0, fpp0S=fpp0S,
            fm0m=fm0m, fm0mS=fm0mS, fm0p=fm0p, fm0pS=fm0pS,
            fp0m=fp0m, fp0mS=fp0mS, fp0p=fp0p, fp0pS=fp0pS,
            f0mm=f0mm, f0mmS=f0mmS, f0mp=f0mp, f0mpS=f0mpS,
            f0pm=f0pm, f0pmS=f0pmS, f0pp=f0pp, f0ppS=f0ppS,
            rho=rho, u=u, v=v, w=w, is_fluid=is_fluid

        # # Iteration über alle Zellen außer die Randzellen
        @inbounds Threads.@threads for z in 2:gridlengthZ-1
            for y in 2:gridlengthY-1
                for x in 2:gridlengthX-1

                    if !is_fluid[x,y,z]
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
                    
                    # Compute equilibrium
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

                    # Compute equilibrium distribution feq
                    feq000 = rho[x,y,z] * P0_u * P0_v * P0_w / 3.0

                    # Face directions
                    feqm00 = rho[x,y,z] * Pm_u * P0_v * P0_w / 18.0
                    feqp00 = rho[x,y,z] * Pp_u * P0_v * P0_w / 18.0

                    feq0m0 = rho[x,y,z] * P0_u * Pm_v * P0_w / 18.0
                    feq0p0 = rho[x,y,z] * P0_u * Pp_v * P0_w / 18.0

                    feq00m = rho[x,y,z] * P0_u * P0_v * Pm_w / 18.0
                    feq00p = rho[x,y,z] * P0_u * P0_v * Pp_w / 18.0

                    # XY edges
                    feqmm0 = rho[x,y,z] * Pm_u * Pm_v * P0_w / 36.0
                    feqmp0 = rho[x,y,z] * Pm_u * Pp_v * P0_w / 36.0
                    feqpm0 = rho[x,y,z] * Pp_u * Pm_v * P0_w / 36.0
                    feqpp0 = rho[x,y,z] * Pp_u * Pp_v * P0_w / 36.0

                    # XZ edges
                    feqm0m = rho[x,y,z] * Pm_u * P0_v * Pm_w / 36.0
                    feqm0p = rho[x,y,z] * Pm_u * P0_v * Pp_w / 36.0
                    feqp0m = rho[x,y,z] * Pp_u * P0_v * Pm_w / 36.0
                    feqp0p = rho[x,y,z] * Pp_u * P0_v * Pp_w / 36.0

                    # YZ edges
                    feq0mm = rho[x,y,z] * P0_u * Pm_v * Pm_w / 36.0
                    feq0mp = rho[x,y,z] * P0_u * Pm_v * Pp_w / 36.0
                    feq0pm = rho[x,y,z] * P0_u * Pp_v * Pm_w / 36.0
                    feq0pp = rho[x,y,z] * P0_u * Pp_v * Pp_w / 36.0 

                    # Compute Non-equilibrium distribution parts
                    fneq000 = f000[x,y,z] - feq000
                    
                    fneqm00 = fm00[x,y,z] - feqm00
                    fneqp00 = fp00[x,y,z] - feqp00

                    fneq0m0 = f0m0[x,y,z] - feq0m0
                    fneq0p0 = f0p0[x,y,z] - feq0p0

                    fneq00m = f00m[x,y,z] - feq00m
                    fneq00p = f00p[x,y,z] - feq00p

                    fneqmm0 = fmm0[x,y,z] - feqmm0
                    fneqmp0 = fmp0[x,y,z] - feqmp0
                    fneqpm0 = fpm0[x,y,z] - feqpm0
                    fneqpp0 = fpp0[x,y,z] - feqpp0

                    fneqm0m = fm0m[x,y,z] - feqm0m
                    fneqm0p = fm0p[x,y,z] - feqm0p
                    fneqp0m = fp0m[x,y,z] - feqp0m
                    fneqp0p = fp0p[x,y,z] - feqp0p

                    fneq0mm = f0mm[x,y,z] - feq0mm
                    fneq0mp = f0mp[x,y,z] - feq0mp
                    fneq0pm = f0pm[x,y,z] - feq0pm
                    fneq0pp = f0pp[x,y,z] - feq0pp

                    # compute Smagorinsky variables 
                    pi_neq_norm = compute_pi_norm(
                        fneqm00, fneqp00, fneq0m0, fneq0p0, fneq00m, fneq00p,
                        fneqmm0, fneqmp0, fneqpm0, fneqpp0,
                        fneqm0m, fneqm0p, fneqp0m, fneqp0p,
                        fneq0mm, fneq0mp, fneq0pm, fneq0pp
                    )

                    omega_local = smagorinsky_omega(τ, CS, pi_neq_norm, rho[x,y,z])

                    # Push scheme: Stream+Collision
                    # Rest particle
                    f000S[x,y,z] = f000[x,y,z] + omega_local * (feq000 - f000[x,y,z])

                    # Face neighbors
                    fm00S[x-1,y,z] = fm00[x,y,z] + omega_local * (feqm00 - fm00[x,y,z])
                    fp00S[x+1,y,z] = fp00[x,y,z] + omega_local * (feqp00 - fp00[x,y,z])

                    f0m0S[x,y-1,z] = f0m0[x,y,z] + omega_local * (feq0m0 - f0m0[x,y,z])
                    f0p0S[x,y+1,z] = f0p0[x,y,z] + omega_local * (feq0p0 - f0p0[x,y,z])

                    f00mS[x,y,z-1] = f00m[x,y,z] + omega_local * (feq00m - f00m[x,y,z])
                    f00pS[x,y,z+1] = f00p[x,y,z] + omega_local * (feq00p - f00p[x,y,z])

                    # XY-plane edges
                    fmm0S[x-1,y-1,z] = fmm0[x,y,z] + omega_local * (feqmm0 - fmm0[x,y,z])
                    fmp0S[x-1,y+1,z] = fmp0[x,y,z] + omega_local * (feqmp0 - fmp0[x,y,z])
                    fpm0S[x+1,y-1,z] = fpm0[x,y,z] + omega_local * (feqpm0 - fpm0[x,y,z])
                    fpp0S[x+1,y+1,z] = fpp0[x,y,z] + omega_local * (feqpp0 - fpp0[x,y,z])

                    # XZ-plane edges
                    fm0mS[x-1,y,z-1] = fm0m[x,y,z] + omega_local * (feqm0m - fm0m[x,y,z])
                    fm0pS[x-1,y,z+1] = fm0p[x,y,z] + omega_local * (feqm0p - fm0p[x,y,z])
                    fp0mS[x+1,y,z-1] = fp0m[x,y,z] + omega_local * (feqp0m - fp0m[x,y,z])
                    fp0pS[x+1,y,z+1] = fp0p[x,y,z] + omega_local * (feqp0p - fp0p[x,y,z])

                    # YZ-plane edges
                    f0mmS[x,y-1,z-1] = f0mm[x,y,z] + omega_local * (feq0mm - f0mm[x,y,z])
                    f0mpS[x,y-1,z+1] = f0mp[x,y,z] + omega_local * (feq0mp - f0mp[x,y,z])
                    f0pmS[x,y+1,z-1] = f0pm[x,y,z] + omega_local * (feq0pm - f0pm[x,y,z])
                    f0ppS[x,y+1,z+1] = f0pp[x,y,z] + omega_local * (feq0pp - f0pp[x,y,z])
        # end #end fluid_nodes
                end #end x
            end #end y
        end #end z
        end #end let (debug)
        end #end elapsed
        
        if i <= 10
            println("Step $i mainloop: $(round(t_debug * 1000, digits=1))ms")
        end

        # bounce-back walls
        @inbounds for wall in wall_indices
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
        # # compute inflow populations fp00S, fpp0S, fpm0S, fp0pS, fp0mS
        # # momentum coefficients for D3Q19 weights
        # # Compute new populations
        fp00S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm00S[1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_face
        fpp0S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fmp0S[1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        fpm0S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fmm0S[1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        fp0pS[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm0pS[1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        fp0mS[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm0mS[1, 2:gridlengthY-1, 2:gridlengthZ-1] .+ inlet_add_edge
        
        # OUTLET: no-gradient bounceback 
        # # all populations that stream in -x direction from previous neighbor
        # # fm00S, fmm0S, fmp0S, fm0mS, fm0pS
        
        fm00S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm00S[gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        fmm0S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fmm0S[gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        fmp0S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fmp0S[gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        fm0mS[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm0mS[gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        fm0pS[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm0pS[gridlengthX-2, 2:gridlengthY-1, 2:gridlengthZ-1]
        
   

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

        # DEBUG STABILITY CHECK
        if EnableStabilityCheck && ((i % CheckEvery == 0) || (i == 1))
            has_nan = any(isnan, u) || any(isnan, v) || any(isnan, w) || any(isnan, rho)
            has_inf = any(isinf, u) || any(isinf, v) || any(isinf, w) || any(isinf, rho)

            umax = maximum(abs, u)
            vmax = maximum(abs, v)
            wmax = maximum(abs, w)

            rhomin = minimum(rho)
            rhomax = maximum(rho)
            rho_total = sum(rho)
            min_rho_location = argmin(rho)

            println("CHK step=$i | umax=$umax vmax=$vmax wmax=$wmax | rho_min=$rhomin rho_max=$rhomax | rho_total = $rho_total")
            println("MIN RHO AT: $min_rho_location")

            if has_nan || has_inf || (umax > VelLimit) || (vmax > VelLimit) || (wmax > VelLimit) ||
            (rhomin < RhoMinLimit) || (rhomax > RhoMaxLimit)
                println("UNSTABLE at step $i")
                println("  has_nan=$has_nan has_inf=$has_inf")
                println("  umax=$umax vmax=$vmax wmax=$wmax")
                println("  rho_min=$rhomin rho_max=$rhomax")
                break
            end
        end


        if (i % 100 == 0) || (i == simulationTime)
            Log_Simulation_Runtime(i, simulationTime)
        end

        # Plot of the field
        if any((Plotvorticity, Plotvx, Plotvy, Plotvz, Plotdebug, Plotmag)) && ((i % 200 == 0) || (i == simulationTime))

            # Log_Simulation_Runtime(i, simulationTime)   
            
            #Copy velocities for plotting
            velocityX .= u
            velocityY .= v
            velocityZ .= w

            velocityMag .= sqrt.(u.^2 .+ v.^2 .+ w.^2)

            # Set velocities inside the sphere to zero
            velocityX[is_object] .= NaN
            velocityY[is_object] .= NaN
            velocityZ[is_object] .= NaN
            velocityMag[is_object] .= NaN         
     
            # Update the observables
            if Plotvorticity == true   
                                    
            end

            if Plotmag==true
                mag_xy_obs[] = copy(velocityMag[2:gridlengthX-1, 2:gridlengthY-1, midZ])
                step_text_mag_xy[] = "Time step: $i / $simulationTime \n $(floor(Int, i*delta_t))s"

                mag_xz_obs[] = copy(velocityMag[2:gridlengthX-1 ,midY, 2:gridlengthZ-1])
                step_text_mag_xz[] = "Time step:$i, $(floor(Int, i*delta_t))s"
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
            sleep(0.01)
        end

    end

    Log_Simulation_Tail()
end
run_JuLattice()