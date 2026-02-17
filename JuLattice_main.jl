############################
## Main file for JuLattice #
############################
include("src/Plotter.jl")
include("src/Logger.jl")

using MeshGrid, GLMakie
using .Plotter, .Logger

function run_JuLattice()
    ####################################  Initialize  ####################################
    ##-------- User Settings --------##
    # Domain Settings
    length_X = 2              # m
    length_Y = 0.5            # m 
    length_Z = 1              # m

    # Sphere Definition
    Radius   = 0.1    # m
    Position = [length_X/4, length_Y/2, length_Z/2] # m [x,y,z]
    # Position = [0.5, 0.25, 0.5] # m [x,y,z]

    # Fluid Settings 
    Fluid_Density = 100.0; #1000.0;         # kg/m^3
    Inflow_Velocity = 0.4;          # m/s
    Kinematic_Viscosity = 0.002; #0.001;    # m^2/s 

    # Simulation Settings
    Simulation_Time = 5;  #8000     # s
    delta_x = 0.01;                 # Grid spacing (physical units per lattice unit)
    Mach_Number = 0.01;          # Target Mach number (Ma = U_lattice/c_s)
                                # Keep Ma < 0.1 for incompressible flow!

    # Compute Reynolds number (for reference)
    Re = (Inflow_Velocity .* 2 .* Radius)/Kinematic_Viscosity;
    Re_Log=floor(Int,Re)


    # Plot Requests (Flags)
    Plotvx = true;
    Plotvy = false;
    Plotvz = false;
    Plotvorticity = false;
    Plotdebug = false;


    ####-------- Run Simulation --------#####
    Log_Simulation_Header()

    ##-------- Compute LBM Parameters from Mach Number --------##
    # fixed lattice constant
    lattice_speedOfSound = 1 / √3; #bleibt gleich bei 3D

    # Step 1: Lattice velocity from Mach number
    lattice_inflow_velocity = Mach_Number * lattice_speedOfSound

    # Step 2: Timestep from velocity scaling
    # U_phys = (dx/dt) * U_lattice => dt = dx * U_lattice / U_phys
    delta_t = delta_x * lattice_inflow_velocity / Inflow_Velocity
    
    # Step 3: Lattice viscosity from physical viscosity
    # nu_phys = (dx²/dt) * nu_lattice => nu_lattice = nu_phys * dt/dx²
    lattice_viscosity = Kinematic_Viscosity * delta_t / (delta_x * delta_x)

    # Step 4: Relaxation time and omega from lattice viscosity
    # nu_lattice = c_s² * (tau - 0.5) => tau = nu_lattice / c_s² + 0.5
    τ = lattice_viscosity / (lattice_speedOfSound * lattice_speedOfSound) + 0.5
    omega  = 1.0 / τ

    ##-------- Convert user settings to lattice units --------##
    # Domain
    gridlengthX = ceil(Int, length_X / delta_x);
    gridlengthY = ceil(Int, length_Y / delta_x);
    gridlengthZ = ceil(Int, length_Z / delta_x);

    # # Sphere
    # sphere_radius  = Radius/delta_x;
    # sphere_position = Position ./ delta_x;

    # Cyliner
    cylinder_x = Int(round(gridlengthX / 3))
    cylinder_y = Int(round(gridlengthY / 2))
    cylinder_radius = Radius/delta_x
    cylinder_start = Int(round(gridlengthZ*0.25))
    cylinder_end = Int(round(gridlengthZ*0.75))

    # Fluid
    fluiddensity = Fluid_Density
    #fluiddensity = 100;

    # ReynoldsCheck
    lattice_Re = (lattice_inflow_velocity .* 2 .* cylinder_radius)/lattice_viscosity; #Re_lattice = U*R/v -> sollte Re entsprechen weil Größen skaliert wurden
    lattice_Re_Log=floor(Int,lattice_Re)

    # Log 
    Log_Discretization_Settings(delta_x, delta_t, lattice_Re_Log)

    # Print τ value
    println("Computed relaxation time τ  = ", round(τ, digits=4))
    println("Computed omega ω = ", round(omega, digits=4))
    println("Reynolds number check:")
    println("   Re (physical) = ", round(Re, digits=2))
    println("   Re (lattice) = ", round(lattice_Re, digits=2))

    # Simulation Settings
    simulationTime = ceil(Int, Simulation_Time / delta_t);

    Q = 19; #D3Q19
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

    #Define array for each direction after Collision+stream (S)
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

    omegaX = zeros(gridlengthX, gridlengthY, gridlengthZ)
    omegaY = zeros(gridlengthX, gridlengthY, gridlengthZ)
    omegaZ = zeros(gridlengthX, gridlengthY, gridlengthZ)
    omegaMag = zeros(gridlengthX, gridlengthY, gridlengthZ)

    #Define mid-Planes for plotting
    midY = ceil(Int, gridlengthY/2)
    midZ = ceil(Int, gridlengthZ/2)

    # more slices for debugging
    frontY = 2
    backY = gridlengthY-1
    botZ = 2
    topZ = gridlengthZ-1

    nearFrontY = 10
    nearBackY = gridlengthY-10
    nearBotZ = 10
    nearTopZ = gridlengthZ-10

    # # DEBUG: Checking for "row bug"
    # initial_u = lattice_inflow_velocity
    # u_debug = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # rho_debug = zeros(gridlengthX, gridlengthY, gridlengthZ)
    # bc_anomaly_detected = Dict{String, Bool}()
    
    # function compute_u_from_fS!(u_out, rho_out)
    #     for x in 2:gridlengthX-1
    #         for y in 2:gridlengthY-1
    #             for z in 2:gridlengthZ-1
    #                 rho_out[x,y,z] = f000S[x,y,z] + 
    #                     (fm00S[x,y,z] + fp00S[x,y,z] + f0m0S[x,y,z] + f0p0S[x,y,z] + f00mS[x,y,z] + f00pS[x,y,z]) +
    #                     (fmm0S[x,y,z] + fmp0S[x,y,z] + fpm0S[x,y,z] + fpp0S[x,y,z] + 
    #                      fm0mS[x,y,z] + fm0pS[x,y,z] + fp0mS[x,y,z] + fp0pS[x,y,z] +
    #                      f0mmS[x,y,z] + f0mpS[x,y,z] + f0pmS[x,y,z] + f0ppS[x,y,z])
                    
    #                 u_out[x,y,z] = ((-fm00S[x,y,z] + fp00S[x,y,z]) +
    #                     (-fmm0S[x,y,z] - fmp0S[x,y,z] + fpm0S[x,y,z] + fpp0S[x,y,z]) +
    #                     (-fm0mS[x,y,z] - fm0pS[x,y,z] + fp0mS[x,y,z] + fp0pS[x,y,z])) / rho_out[x,y,z]
    #             end
    #         end
    #     end
    # end

    # function check_u_rows(bc_name, timestep, y_idx, z_range, init_u, bc_detected)
    #     key = bc_name
    #     if haskey(bc_detected, key) && bc_detected[key]
    #         return
    #     end
        
    #     for z in z_range
    #         row = u_debug[2:end-1, y_idx, z]
    #         if length(unique(round.(row, digits=6))) == 1
    #             row_val = row[1]
    #             if abs(row_val) < 1e-10
    #                 continue
    #             end
    #             diff = abs(row_val - init_u)
    #             if diff > 0.1*abs(init_u)
    #                 println(" HÄNDE HOCH:")
    #                 println("   Timestep: $timestep")
    #                 println("   BC: $bc_name")
    #                 println("   Position: y=$y_idx, z=$z")
    #                 println("   u_initial: $init_u")
    #                 println("   u_after_BC: $row_val")
    #                 println("   Difference: $diff")
    #                 println("")
    #                 # bc_detected[key] = true
    #                 return
    #             end
    #         end
    #     end
    # end

    # function check_u_rows_y(bc_name, timestep, y_range, z_idx, init_u, bc_detected)
    #     key = bc_name
    #     if haskey(bc_detected, key) && bc_detected[key]
    #         return
    #     end
        
    #     for y in y_range
    #         row = u_debug[2:end-1, y, z_idx]
    #         if length(unique(round.(row, digits=6))) == 1
    #             row_val = row[1]
    #             if abs(row_val) < 1e-10
    #                 continue
    #             end
    #             diff = abs(row_val - init_u)
    #             if diff > 0.1*abs(init_u)
    #                 println(" HÄNDE HOCH")
    #                 println("   Timestep: $timestep")
    #                 println("   BC: $bc_name")
    #                 println("   Position: y=$y, z=$z_idx")
    #                 println("   u_initial: $init_u")
    #                 println("   u_after_BC: $row_val")
    #                 println("   Difference: $diff")
    #                 println("")
    #                 # bc_detected[key] = true
    #                 return
    #             end
    #         end
    #     end
    # end




    # # create object indetifier
    # sphere = (gridX .- sphere_position[1]).^2 + (gridY .- sphere_position[2]).^2 + (gridZ .- sphere_position[3]).^2 .< sphere_radius.^2
    # sphere_indices = findall(sphere)

    # # create boundary indetifiers
    # # walls = gridY .== 1 .| gridY .== gridlengthY .| gridZ .==1 .| gridZ .== gridlengthZ;
    # inlet = gridX .== 1;
    # outlet = gridX .== gridlengthX;

    # create solid node mask
    is_solid = falses(gridlengthX, gridlengthY, gridlengthZ)
    for x in 1:gridlengthX, y in 1:gridlengthY, z in 1:gridlengthZ
        # walls
        # if x==1 || x==gridlengthX || y==1 || y==gridlengthY || z==1 || z==gridlengthZ
        if y==1 || y==gridlengthY || z==1 || z==gridlengthZ

            is_solid[x, y, z] = true
            continue
        end
        
        # cylinder vertically (y-axis)

        dx = x- cylinder_x
        dy = y - cylinder_y
        if (z >= cylinder_start) && (z <= cylinder_end) && (sqrt(dx^2 + dy^2) <= cylinder_radius)
            is_solid[x, y, z] = true
        end

    end

    #Initialize distribution functions FLUID NODES and SOLID NODES
    for x in 1:gridlengthX
        for y in 1:gridlengthY
            for z in 1:gridlengthZ

                if is_solid[x, y, z]
                    ux = 0.0
                else
                    ux = lattice_inflow_velocity
                end

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
                # is_solid = (x==1 || x==gridlengthX || y==1 || y==gridlengthY || z==1 || z==gridlengthZ)

                # ux = is_solid ? 0.0 : lattice_inflow_velocity
                ux = lattice_inflow_velocity
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

    #Plot calls
    if any((Plotvorticity, Plotvx, Plotvy, Plotvz, Plotdebug))
        # #3D Plot
        # if Plotvorticity == true
        #     omegaMag_obs, step_text_omega, fig_omega = Create_Plot3D(gridlengthX, gridlengthY, gridlengthZ, omegaMag; title="|ω|")
        #     screenOmega = GLMakie.Screen()
        #     display(screenOmega, fig_omega)
        # end

        if Plotvx==true
            #xy slice at z=midZ
            vx_xy_obs, step_text_vx_xy, fig_vx_xy = Create_Plot_XY(gridlengthX, gridlengthY, velocityX[:,:,midZ]; title="v_x at z=$(midZ)")
            screen_vx_xy = GLMakie.Screen()
            display(screen_vx_xy, fig_vx_xy)

            #xz slice at y=midY
            vx_xz_obs, step_text_vx_xz, fig_vx_xz = Create_Plot_XZ(gridlengthX, gridlengthZ, velocityX[:,midY,:]; title="v_x at y=$(midY)")
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


        

    end

    println("#################################")
    println("Starting Simulation:")
    # Run Simulation Loop
    for i in 1:simulationTime
        ###### NEW STABILIZATION ######
        for x in 2:gridlengthX-1 #Iteration über alle Zellen außer die Randzellen
            for y in 2:gridlengthY-1
                for z in 2:gridlengthZ-1

                    # Check: Is node in sphere
                    # is_in_sphere = ((x- sphere_position[1])^2 +
                    #                 (y - sphere_position[2])^2 +
                    #                 (z - sphere_position[3])^2) <= sphere_radius^2

                    # if !is_in_sphere

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
                        #P0_u = -2 + 3*u2
                        P0_u = 1 - 1.5*u2
                        Pp_u = 1 + 3*u[x,y,z] + 3*u2
                        
                        Pm_v = 1 - 3*v[x,y,z] + 3*v2
                        #P0_v = -2 + 3*v2
                        P0_v = 1 - 1.5*v2
                        Pp_v = 1 + 3*v[x,y,z] + 3*v2
                        
                        Pm_w = 1 - 3*w[x,y,z] + 3*w2
                        #P0_w = -2 + 3*w2
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
                    # end

                end
            end
        end

        ###### Boundary Conditions ######
        ## post collision Loop for bounce-back
        for x in 1:gridlengthX, y in 1:gridlengthY, z in 1:gridlengthZ
            if is_solid[x, y, z]
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
        end

        ## periodic inlet/outlet
        # INLET  (left side)
        fp00S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fp00S[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1]
        fp0pS[2, 2:gridlengthY-1, 3:gridlengthZ-1] .= fp0pS[gridlengthX, 2:gridlengthY-1, 3:gridlengthZ-1]
        fp0mS[2, 2:gridlengthY-1, 2:gridlengthZ-2] .= fp0mS[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-2]
        fpm0S[2, 2:gridlengthY-2, 2:gridlengthZ-1] .= fpm0S[gridlengthX, 2:gridlengthY-2, 2:gridlengthZ-1]
        fpp0S[2, 3:gridlengthY-1, 2:gridlengthZ-1] .= fpp0S[gridlengthX, 3:gridlengthY-1, 2:gridlengthZ-1]
        
        # OUTLET (right side)
        fm00S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm00S[1, 2:gridlengthY-1, 2:gridlengthZ-1]
        fm0pS[gridlengthX-1, 2:gridlengthY-1, 3:gridlengthZ-1] .= fm0pS[1, 2:gridlengthY-1, 3:gridlengthZ-1]
        fm0mS[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-2] .= fm0mS[1, 2:gridlengthY-1, 2:gridlengthZ-2]
        fmm0S[gridlengthX-1, 2:gridlengthY-2, 2:gridlengthZ-1] .= fmm0S[1, 2:gridlengthY-2, 2:gridlengthZ-1]
        fmp0S[gridlengthX-1, 3:gridlengthY-1, 2:gridlengthZ-1] .= fmp0S[1, 3:gridlengthY-1, 2:gridlengthZ-1]
        ## post collision Loop
        


        # ## Sphere Bounce-Back 
  
        # f000S[sphere_indices] = f000[sphere_indices]
        
        # fm00S[sphere_indices] = fp00[sphere_indices]
        # fp00S[sphere_indices] = fm00[sphere_indices]
        # f0m0S[sphere_indices] = f0p0[sphere_indices]
        # f0p0S[sphere_indices] = f0m0[sphere_indices]
        # f00mS[sphere_indices] = f00p[sphere_indices]
        # f00pS[sphere_indices] = f00m[sphere_indices]
        
        # # XY-plane edges
        # fmm0S[sphere_indices] = fpp0[sphere_indices]
        # fpp0S[sphere_indices] = fmm0[sphere_indices]
        # fmp0S[sphere_indices] = fpm0[sphere_indices]
        # fpm0S[sphere_indices] = fmp0[sphere_indices]
        
        # # XZ-plane edges
        # fm0mS[sphere_indices] = fp0p[sphere_indices]
        # fp0pS[sphere_indices] = fm0m[sphere_indices]
        # fm0pS[sphere_indices] = fp0m[sphere_indices]
        # fp0mS[sphere_indices] = fm0p[sphere_indices]
        
        # # YZ-plane edges
        # f0mmS[sphere_indices] = f0pp[sphere_indices]
        # f0ppS[sphere_indices] = f0mm[sphere_indices]
        # f0mpS[sphere_indices] = f0pm[sphere_indices]
        # f0pmS[sphere_indices] = f0mp[sphere_indices]


         
        
        # #periodic inlet/outlet overwriting
        # # INLET  (left side)
        # fp00S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fp00S[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fp0pS[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fp0pS[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fp0mS[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fp0mS[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fpm0S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fpm0S[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fpp0S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fpp0S[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1]
        
        # # OUTLET (right side)
        # fm00S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm00S[1, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fm0pS[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm0pS[1, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fm0mS[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm0mS[1, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fmm0S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fmm0S[1, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fmp0S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fmp0S[1, 2:gridlengthY-1, 2:gridlengthZ-1]

        # periodic inlet/outlet
        # # INLET  (left side)
        # fp00S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fp00S[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fp0pS[2, 2:gridlengthY-1, 3:gridlengthZ-1] .= fp0pS[gridlengthX, 2:gridlengthY-1, 3:gridlengthZ-1]
        # fp0mS[2, 2:gridlengthY-1, 2:gridlengthZ-2] .= fp0mS[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-2]
        # fpm0S[2, 2:gridlengthY-2, 2:gridlengthZ-1] .= fpm0S[gridlengthX, 2:gridlengthY-2, 2:gridlengthZ-1]
        # fpp0S[2, 3:gridlengthY-1, 2:gridlengthZ-1] .= fpp0S[gridlengthX, 3:gridlengthY-1, 2:gridlengthZ-1]
        
        # # OUTLET (right side)
        # fm00S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm00S[1, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fm0pS[gridlengthX-1, 2:gridlengthY-1, 3:gridlengthZ-1] .= fm0pS[1, 2:gridlengthY-1, 3:gridlengthZ-1]
        # fm0mS[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-2] .= fm0mS[1, 2:gridlengthY-1, 2:gridlengthZ-2]
        # fmm0S[gridlengthX-1, 2:gridlengthY-2, 2:gridlengthZ-1] .= fmm0S[1, 2:gridlengthY-2, 2:gridlengthZ-1]
        # fmp0S[gridlengthX-1, 3:gridlengthY-1, 2:gridlengthZ-1] .= fmp0S[1, 3:gridlengthY-1, 2:gridlengthZ-1]
        
        # # front wall (y=1)
        # f0p0S[2:gridlengthX-1, 2, 2:gridlengthZ-1] .= f0m0S[2:gridlengthX-1, 1, 2:gridlengthZ-1]
        # fpp0S[2:gridlengthX-1, 2, 2:gridlengthZ-1] .= fmm0S[1:gridlengthX-2, 1, 2:gridlengthZ-1]
        # fmp0S[2:gridlengthX-1, 2, 2:gridlengthZ-1] .= fpm0S[3:gridlengthX, 1, 2:gridlengthZ-1]

        #     # double counted populations (special handling at top/bot)
        # f0pmS[2:gridlengthX-1, 2, 2:gridlengthZ-2] .= f0mpS[2:gridlengthX-1, 1, 3:gridlengthZ-1]
        # f0ppS[2:gridlengthX-1, 2, 3:gridlengthZ-1] .= f0mmS[2:gridlengthX-1, 1, 2:gridlengthZ-2]

        # # back wall (y=gridlengthY)
        # f0m0S[2:gridlengthX-1, gridlengthY-1, 2:gridlengthZ-1] .= f0p0S[2:gridlengthX-1, gridlengthY, 2:gridlengthZ-1]
        # fpm0S[2:gridlengthX-1, gridlengthY-1, 2:gridlengthZ-1] .= fmp0S[1:gridlengthX-2, gridlengthY, 2:gridlengthZ-1]
        # # fmm0S[2:gridlengthX-1, gridlengthY-1, 2:gridlengthZ-1] .= fpp0S[3:gridlengthX, gridlengthY, 2:gridlengthZ-1]

        #     # double counted populations (special handling at top/bot)
        # f0mmS[2:gridlengthX-1, gridlengthY-1, 2:gridlengthZ-2] .= f0ppS[2:gridlengthX-1, gridlengthY, 3:gridlengthZ-1]
        # f0mpS[2:gridlengthX-1, gridlengthY-1, 3:gridlengthZ-1] .= f0pmS[2:gridlengthX-1, gridlengthY, 2:gridlengthZ-2]
        
        # ## Bounceback walls
        # # bottom wall (z=1)
        # f00pS[2:gridlengthX-1, 2:gridlengthY-1, 2] .= f00mS[2:gridlengthX-1, 2:gridlengthY-1, 1]    # mitte
        # f0ppS[2:gridlengthX-1, 2:gridlengthY-1, 2] .= f0mmS[2:gridlengthX-1, 1:gridlengthY-2, 1]
        # f0mpS[2:gridlengthX-1, 2:gridlengthY-1, 2] .= f0pmS[2:gridlengthX-1, 3:gridlengthY, 1]
        # fp0pS[2:gridlengthX-1, 2:gridlengthY-1, 2] .= fm0mS[1:gridlengthX-2, 2:gridlengthY-1, 1]
        # fm0pS[2:gridlengthX-1, 2:gridlengthY-1, 2] .= fp0mS[3:gridlengthX, 2:gridlengthY-1, 1]

        # # # DEBUG: Check after bottom wall
        # # compute_u_from_fS!(u_debug, rho_debug)
        # # check_u_rows("bot", i, frontY, 2:5, initial_u, bc_anomaly_detected)
        # # check_u_rows("bot", i, backY, 2:5, initial_u, bc_anomaly_detected)

        # # top wall (z=gridlengthZ)
        # f00mS[2:gridlengthX-1, 2:gridlengthY-1, gridlengthZ-1] .= f00pS[2:gridlengthX-1, 2:gridlengthY-1, gridlengthZ]  # mitte
        # f0pmS[2:gridlengthX-1, 2:gridlengthY-1, gridlengthZ-1] .= f0mpS[2:gridlengthX-1, 1:gridlengthY-2, gridlengthZ]
        # f0mmS[2:gridlengthX-1, 2:gridlengthY-1, gridlengthZ-1] .= f0ppS[2:gridlengthX-1, 3:gridlengthY, gridlengthZ]
        # fp0mS[2:gridlengthX-1, 2:gridlengthY-1, gridlengthZ-1] .= fm0pS[1:gridlengthX-2, 2:gridlengthY-1, gridlengthZ]
        # fm0mS[2:gridlengthX-1, 2:gridlengthY-1, gridlengthZ-1] .= fp0pS[3:gridlengthX, 2:gridlengthY-1, gridlengthZ]

        # #  # DEBUG: Check after top wall
        # # compute_u_from_fS!(u_debug, rho_debug)
        # # check_u_rows("top", i, frontY, (gridlengthZ-5):(gridlengthZ-1), initial_u, bc_anomaly_detected)
        # # check_u_rows("top", i, backY, (gridlengthZ-5):(gridlengthZ-1), initial_u, bc_anomaly_detected)

        # # front wall (y=1)
        # f0p0S[2:gridlengthX-1, 2, 2:gridlengthZ-1] .= f0m0S[2:gridlengthX-1, 1, 2:gridlengthZ-1]
        # fpp0S[2:gridlengthX-1, 2, 2:gridlengthZ-1] .= fmm0S[1:gridlengthX-2, 1, 2:gridlengthZ-1]
        # fmp0S[2:gridlengthX-1, 2, 2:gridlengthZ-1] .= fpm0S[3:gridlengthX, 1, 2:gridlengthZ-1]

        #     # double counted populations (special handling at top/bot)
        # f0pmS[2:gridlengthX-1, 2, 2:gridlengthZ-2] .= f0mpS[2:gridlengthX-1, 1, 3:gridlengthZ-1]
        # f0ppS[2:gridlengthX-1, 2, 3:gridlengthZ-1] .= f0mmS[2:gridlengthX-1, 1, 2:gridlengthZ-2]

        # # # DEBUG: Check after front wall
        # # compute_u_from_fS!(u_debug, rho_debug)
        # # check_u_rows_y("front", i, 2:5, midZ, initial_u, bc_anomaly_detected)

        # # back wall (y=gridlengthY)
        # f0m0S[2:gridlengthX-1, gridlengthY-1, 2:gridlengthZ-1] .= f0p0S[2:gridlengthX-1, gridlengthY, 2:gridlengthZ-1]
        # fpm0S[2:gridlengthX-1, gridlengthY-1, 2:gridlengthZ-1] .= fmp0S[1:gridlengthX-2, gridlengthY, 2:gridlengthZ-1]
        # fmm0S[2:gridlengthX-1, gridlengthY-1, 2:gridlengthZ-1] .= fpp0S[3:gridlengthX, gridlengthY, 2:gridlengthZ-1]

        #     # double counted populations (special handling at top/bot)
        # f0mmS[2:gridlengthX-1, gridlengthY-1, 2:gridlengthZ-2] .= f0ppS[2:gridlengthX-1, gridlengthY, 3:gridlengthZ-1]
        # f0mpS[2:gridlengthX-1, gridlengthY-1, 3:gridlengthZ-1] .= f0pmS[2:gridlengthX-1, gridlengthY, 2:gridlengthZ-2]

        # # # DEBUG: Check after back wall
        # # compute_u_from_fS!(u_debug, rho_debug)
        # # check_u_rows_y("back", i, (gridlengthY-5):(gridlengthY-2), midZ, initial_u, bc_anomaly_detected)


        # # periodic inlet/outlet
        # # INLET  (left side)
        # fp00S[2, 2:gridlengthY-1, 2:gridlengthZ-1] .= fp00S[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fp0pS[2, 2:gridlengthY-1, 3:gridlengthZ-1] .= fp0pS[gridlengthX, 2:gridlengthY-1, 3:gridlengthZ-1]
        # fp0mS[2, 2:gridlengthY-1, 2:gridlengthZ-2] .= fp0mS[gridlengthX, 2:gridlengthY-1, 2:gridlengthZ-2]
        # fpm0S[2, 2:gridlengthY-2, 2:gridlengthZ-1] .= fpm0S[gridlengthX, 2:gridlengthY-2, 2:gridlengthZ-1]
        # fpp0S[2, 3:gridlengthY-1, 2:gridlengthZ-1] .= fpp0S[gridlengthX, 3:gridlengthY-1, 2:gridlengthZ-1]
        
        # # OUTLET (right side)
        # fm00S[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-1] .= fm00S[1, 2:gridlengthY-1, 2:gridlengthZ-1]
        # fm0pS[gridlengthX-1, 2:gridlengthY-1, 3:gridlengthZ-1] .= fm0pS[1, 2:gridlengthY-1, 3:gridlengthZ-1]
        # fm0mS[gridlengthX-1, 2:gridlengthY-1, 2:gridlengthZ-2] .= fm0mS[1, 2:gridlengthY-1, 2:gridlengthZ-2]
        # fmm0S[gridlengthX-1, 2:gridlengthY-2, 2:gridlengthZ-1] .= fmm0S[1, 2:gridlengthY-2, 2:gridlengthZ-1]
        # fmp0S[gridlengthX-1, 3:gridlengthY-1, 2:gridlengthZ-1] .= fmp0S[1, 3:gridlengthY-1, 2:gridlengthZ-1]


        # # Sphere Bounceback
        # # Cell faces - swap (vectorized)
        # fp00S[sphere_indices], fm00S[sphere_indices] = fm00S[sphere_indices], fp00S[sphere_indices]
        # f0p0S[sphere_indices], f0m0S[sphere_indices] = f0m0S[sphere_indices], f0p0S[sphere_indices]
        # f00pS[sphere_indices], f00mS[sphere_indices] = f00mS[sphere_indices], f00pS[sphere_indices]

        # # XY-plane edges (vectorized)
        # fpp0S[sphere_indices], fmm0S[sphere_indices] = fmm0S[sphere_indices], fpp0S[sphere_indices]
        # fpm0S[sphere_indices], fmp0S[sphere_indices] = fmp0S[sphere_indices], fpm0S[sphere_indices]

        # # XZ-plane edges (vectorized)
        # fp0pS[sphere_indices], fm0mS[sphere_indices] = fm0mS[sphere_indices], fp0pS[sphere_indices]
        # fp0mS[sphere_indices], fm0pS[sphere_indices] = fm0pS[sphere_indices], fp0mS[sphere_indices]

        # # YZ-plane edges (vectorized)
        # f0ppS[sphere_indices], f0mmS[sphere_indices] = f0mmS[sphere_indices], f0ppS[sphere_indices]
        # f0pmS[sphere_indices], f0mpS[sphere_indices] = f0mpS[sphere_indices], f0pmS[sphere_indices]
        

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


        # ## Reset solid cells to zeros
        # all_dir = [f000, fm00, fp00, f0m0, f0p0, f00m, f00p,
        #             fmm0, fmp0, fpm0, fpp0, fm0m, fm0p, fp0m, fp0p,
        #             f0mm, f0mp, f0pm, f0pp]

        # for dir in all_dir
        #     dir[[1, gridlengthX], 2:gridlengthY-1, 2:gridlengthZ-1] .= 0.0
        #     dir[2:gridlengthX-1, [1, gridlengthY], 2:gridlengthZ-1] .= 0.0
        #     dir[2:gridlengthX-1, 2:gridlengthY-1, [1, gridlengthZ]] .= 0.0
        # end

        ###### NEW STABILIZATION #######

        # ------------------------ vorerst ohne ------------------------ 
        # ## Apply Boundary conditions
        # #Inlet velocity bc (unknown: f_1, f_8, f_9)
        # densityGrid[inlet, :] .= (sum(distributions[inlet, [1,3,5]], dims=2).+ 2 .*sum(distributions[inlet, [2,6,7]], dims=2)) ./ (1-lattice_inflow_velocity)
        # distributions[inlet, 4] .= distributions[inlet, 2] .+ (2/3 .* densityGrid[inlet,:] .* lattice_inflow_velocity)
        # distributions[inlet, 8] .= distributions[inlet, 6] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .- (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))
        # distributions[inlet, 9] .= distributions[inlet, 7] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .+ (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))
        

        # #Outlet zero gradient bc
        # distributions[outlet, [4, 8, 9]] .= distributions[gridlengthX-1, :, [4, 8, 9]]
        # ------------------------ vorerst ohne ------------------------ 


        # Plot of the field
        # if ((i % 10 == 0)) || (i == simulationTime)
        if ((i % 20 == 0)) || (i == simulationTime)


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
                vx_xy_obs[] = copy(velocityX[:,:,midZ])
                step_text_vx_xy[] = "Time step:$i, $(floor(Int, i*delta_t))s"

                #refresh xz slice at y=midY
                vx_xz_obs[] = copy(velocityX[:,midY,:])
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