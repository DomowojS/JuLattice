############################
## Main file for JuLattice #
############################
include("src/Plotter.jl")
include("src/Logger.jl")

using MeshGrid, GLMakie
using .Plotter, .Logger

using Random

function run_JuLattice()
    ####################################  Initialize  ####################################
    ## User Settings
    # Domain Settings
    length_X = 4;              # m
    length_Y = 1;              # m 

    # Cylinder Definition
    Radius   = 0.1    # m
    Position = [1, 0.5] # m

    # Fluid Settings 
    fluiddensity = 100.0; #1000.0       # kg/m^3
    Inflow_Velocity = 0.4;      # m/s
    Kinematic_Viscosity = 0.001; # m^2/s 
    

    # Simulation Settings
    Simulation_Time = 120;     # s
    delta_x = 0.01;             # Grid spacing (physical units per lattice unit)
    Mach_Number = 0.1;          # Target Mach number (Ma = U_lattice/c_s)
                                # Keep Ma < 0.1 for incompressible flow!

    # Compute Reynolds number (for reference)
    Re = (Inflow_Velocity .* Radius)/Kinematic_Viscosity;
    Re_Log=floor(Int,Re)

    # Plot Requests
    Plotvx = true;
    Plotvy = true;
    Plotvorticity = true;

    #### Run Simulation #####
    Log_Simulation_Header()

    ##-------- Compute LBM Parameters from Mach Number --------##
    # Fixed lattice constant
    lattice_speedOfSound = 1 / √3;
    
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


    ## Convert user settings to lattice units
    # Domain
    gridlengthX  = ceil(Int, length_X / delta_x);
    gridlengthY  = ceil(Int, length_Y / delta_x);

    # Cylinder
    cylinder_radius  = Radius/delta_x;
    cylinder_position = Position ./ delta_x;

    # Verify Reynolds number consistency (lattice vs physical)
    lattice_Re = (lattice_inflow_velocity .* cylinder_radius)/lattice_viscosity;
    lattice_Re_Log=floor(Int,lattice_Re)

    #Log 
    Log_Discretization_Settings(delta_x, delta_t, lattice_Re_Log)

    # Simulation Settings
    simulationTime = ceil(Int, Simulation_Time / delta_t);

    # Output Simulation Time
    sim_minutes = floor(Int, Simulation_Time / 60)
    sim_seconds = Simulation_Time % 60
    println("physical simulation time: $(sim_minutes) min $(sim_seconds) s")
    println("Timesteps: $(simulationTime)")

    #Define arrays for each direction D2Q9
    f00 = zeros(gridlengthX, gridlengthY) #center
    fm0 = zeros(gridlengthX, gridlengthY) #left
    f0m = zeros(gridlengthX, gridlengthY) #down
    fp0 = zeros(gridlengthX, gridlengthY) #right
    f0p = zeros(gridlengthX, gridlengthY) #up
    fmm = zeros(gridlengthX, gridlengthY) #left-down
    fmp = zeros(gridlengthX, gridlengthY) #left-up
    fpp = zeros(gridlengthX, gridlengthY) #right-up
    fpm = zeros(gridlengthX, gridlengthY) #right-down

    #Define array for each direction after stream+Collision (S)
    f00S = zeros(gridlengthX, gridlengthY) #center
    fm0S = zeros(gridlengthX, gridlengthY) #left
    f0mS = zeros(gridlengthX, gridlengthY) #down
    fp0S = zeros(gridlengthX, gridlengthY) #right
    f0pS = zeros(gridlengthX, gridlengthY) #up
    fmmS = zeros(gridlengthX, gridlengthY) #left-down
    fmpS = zeros(gridlengthX, gridlengthY) #left-up
    fppS = zeros(gridlengthX, gridlengthY) #right-up
    fpmS = zeros(gridlengthX, gridlengthY) #right-down

    #Initialise macroscopic variables
    rho = ones(gridlengthX, gridlengthY) .* fluiddensity
    u = zeros(gridlengthX,gridlengthY) #vx
    v = zeros(gridlengthX, gridlengthY) #vy

    #define omega
    omega = 1.0 / τ

    #Initialise distribution functions 
    for x in 1:gridlengthX
        for y in 1:gridlengthY
            #same ux and uy for all cells so no arrays for ux, uy
            ux = lattice_inflow_velocity
            uy = 0.0
            u[x,y] = ux
            v[x,y] = uy

            rho_init = fluiddensity

            f00[x,y] = rho_init * (-2.0 + 3.0*ux*ux) * (-2.0 + 3.0*uy*uy) / 9.0

            fm0[x,y] = rho_init * (1.0 - 3.0*ux + 3.0*ux*ux) * (-2.0 + 3.0*uy*uy) / -18.0 #statt / 18
            fp0[x,y] = rho_init * (1.0 + 3.0*ux + 3.0*ux*ux) * (-2.0 + 3.0*uy*uy) / -18.0
            f0m[x,y] = rho_init * (-2.0 + 3.0*ux*ux) * (1.0 - 3.0*uy + 3.0*uy*uy) / -18.0
            f0p[x,y] = rho_init * (-2.0 + 3.0*ux*ux) * (1.0 + 3.0*uy + 3.0*uy*uy) / -18.0

            fmm[x,y] = rho_init * (1.0 - 3.0*ux + 3.0*ux*ux) * (1.0 - 3.0*uy + 3.0*uy*uy) / 36.0
            fmp[x,y] = rho_init * (1.0 - 3.0*ux + 3.0*ux*ux) * (1.0 + 3.0*uy + 3.0*uy*uy) / 36.0
            fpm[x,y] = rho_init * (1.0 + 3.0*ux + 3.0*ux*ux) * (1.0 - 3.0*uy + 3.0*uy*uy) / 36.0
            fpp[x,y] = rho_init * (1.0 + 3.0*ux + 3.0*ux*ux) * (1.0 + 3.0*uy + 3.0*uy*uy) / 36.0       
        end
    end

    f00S .= f00

    fm0S .= fm0
    fp0S .= fp0
    f0mS .= f0m
    f0pS .= f0p

    fmmS .= fmm
    fmpS .= fmp
    fpmS .= fpm
    fppS .= fpp

    #DEBUG Density check

    for x in 1:gridlengthX
        for y in 1:gridlengthY
            
            rho_check = f00[x,y] + fm0[x,y] + fp0[x,y] + f0m[x,y] + f0p[x,y] +
                        fmm[x,y] + fmp[x,y] + fpm[x,y] + fpp[x,y]
            
            rho_fluiddensity = fluiddensity

            if abs(rho_check - rho_fluiddensity) > 0.01 * rho_fluiddensity
                println("=== DENSITY CHECK ===")
                println("rho_check: $rho_check , fluiddensity: $rho_fluiddensity")
            end
        end
    end

    # create grid
    gridX, gridY = meshgrid(1:gridlengthX, 1:gridlengthY);
    gridX, gridY = gridX', gridY';

    # create object indetifier
    cylinder = (gridX.-cylinder_position[1]).^2 + (gridY.-cylinder_position[2]).^2 .< cylinder_radius.^2;
    cylinder_indices = findall(cylinder)

    # create boundary indetifiers
    walls = gridY .== 1 .| gridY .== gridlengthY;
    inlet = gridX .== 1;
    outlet = gridX .== gridlengthX;

    if any((Plotvorticity, Plotvx, Plotvy))
        if Plotvorticity==true 
            vorticity, vorticity_obs, text_obj, step_text, fig_vorticity = Create_Plot(gridlengthX, gridlengthY)
            screen1 = GLMakie.Screen()
            GLMakie.display(screen1, fig_vorticity)
        end
        if Plotvx==true 
            velocityX_obs, text_obj_vx, step_text_vx, fig_vx = Create_Plot(gridlengthX, gridlengthY, u, "X")
            screen2 = GLMakie.Screen(; position = (600, 0))
            GLMakie.display(screen2, fig_vx)
        end
        if Plotvy==true 
            velocityY_obs, text_obj_vy, step_text_vy, fig_vy = Create_Plot(gridlengthX, gridlengthY, v, "Y")
            screen3 = GLMakie.Screen()
            display(screen3, fig_vy)
        end
    end
    ####################################  Initialize  ####################################


    ####################################  Sim-Loop ####################################
    println("#################################")
    println("Starting Simulation:")
    # Run Simulation Loop
    for i in 1:simulationTime

        for x in 2:gridlengthX-1
            for y in 2:gridlengthY-1
                # Get Macroscopic values
                rho[x, y] = f00[x, y] + (((fmm[x, y] + fpp[x, y]) + (fmp[x, y] + fpm[x, y])) + 
                                        ((fm0[x, y] + fp0[x, y]) + (f0p[x, y] + f0m[x, y])))
                
                u[x, y] = (((-fmm[x, y] + fpp[x, y]) + (-fmp[x, y] + fpm[x, y])) + 
                        ((-fm0[x, y] + fp0[x, y]))) / rho[x, y]
                
                v[x, y] = (((-fmm[x, y] + fpp[x, y]) + (fmp[x, y] - fpm[x, y])) + 
                        ((f0p[x, y] - f0m[x, y]))) / rho[x, y]
                
                # Push and Collision scheme
                # Collision is computed and initialised as the corresponding cell after Push
                fmmS[x-1, y-1] = fmm[x, y] + omega * ((rho[x, y] * (1 - 3*u[x, y] + 3*(u[x, y]*u[x, y])) * 
                                                    (1 - 3*v[x, y] + 3*(v[x, y]*v[x, y]))) / 36.0 - fmm[x, y])
                
                f0mS[x, y-1] = f0m[x, y] + omega * (-0.05555555555555555 * ((-2 + 3*(u[x, y]*u[x, y])) * 
                                                    rho[x, y] * (1 + 3*(v[x, y]*v[x, y]) - 3*v[x, y])) - f0m[x, y])
                
                fpmS[x+1, y-1] = fpm[x, y] + omega * ((rho[x, y] * (1 + 3*(u[x, y]*u[x, y]) + 3*u[x, y]) * 
                                                    (1 + 3*(v[x, y]*v[x, y]) - 3*v[x, y])) / 36.0 - fpm[x, y])
                
                fm0S[x-1, y] = fm0[x, y] + omega * (-0.05555555555555555 * ((-2 + 3*(v[x, y]*v[x, y])) * 
                                                    rho[x, y] * (1 + 3*(u[x, y]*u[x, y]) - 3*u[x, y])) - fm0[x, y])
                
                f00S[x, y] = f00[x, y] + omega * (((-2 + 3*(u[x, y]*u[x, y])) * (-2 + 3*(v[x, y]*v[x, y])) * 
                                                rho[x, y]) / 9.0 - f00[x, y])
                
                fp0S[x+1, y] = fp0[x, y] + omega * (-0.05555555555555555 * ((-2 + 3*(v[x, y]*v[x, y])) * 
                                                    rho[x, y] * (1 + 3*(u[x, y]*u[x, y]) + 3*u[x, y])) - fp0[x, y])
                
                fmpS[x-1, y+1] = fmp[x, y] + omega * ((rho[x, y] * (1 + 3*(u[x, y]*u[x, y]) - 3*u[x, y]) * 
                                                    (1 + 3*(v[x, y]*v[x, y]) + 3*v[x, y])) / 36.0 - fmp[x, y])
                
                f0pS[x, y+1] = f0p[x, y] + omega * (-0.05555555555555555 * ((-2 + 3*(u[x, y]*u[x, y])) * 
                                                    rho[x, y] * (1 + 3*(v[x, y]*v[x, y]) + 3*v[x, y])) - f0p[x, y])
                
                fppS[x+1, y+1] = fpp[x, y] + omega * ((rho[x, y] * (1 + 3*(u[x, y]*u[x, y]) + 3*u[x, y]) * 
                                                    (1 + 3*(v[x, y]*v[x, y]) + 3*v[x, y])) / 36.0 - fpp[x, y])
            end
        end


        ##### Boundary Conditions #####
        ## Periodic inlet / outlet
        # inlet (left)
        fp0S[2, 2:gridlengthY-1] .= fp0S[gridlengthX, 2:gridlengthY-1]
        fppS[2, 2:gridlengthY-1] .= fppS[gridlengthX, 2:gridlengthY-1]
        fpmS[2, 2:gridlengthY-1] .= fpmS[gridlengthX, 2:gridlengthY-1]

        # outlet (right)
        fm0S[gridlengthX-1, 2:gridlengthY-1] .= fm0S[1, 2:gridlengthY-1]
        fmpS[gridlengthX-1, 2:gridlengthY-1] .= fmpS[1, 2:gridlengthY-1]
        fmmS[gridlengthX-1, 2:gridlengthY-1] .= fmmS[1, 2:gridlengthY-1]

        ## Bounceback walls
        # bottom wall (y=1) streams in (y=2)
        f0pS[2:gridlengthX-1, 2] .= f0mS[2:gridlengthX-1, 1]        # top = bottom
        # diagonals streaming back to origin node
        fppS[2:gridlengthX-1, 2] .= fmmS[1:gridlengthX-2, 1]        # right-top = left-bottom
        fmpS[2:gridlengthX-1, 2] .= fpmS[3:gridlengthX,   1]        # left-top = right-bottom                          

        # top wall (y=gridlengthY) streams in (y=gridlengthY-1)
        f0mS[2:gridlengthX-1, gridlengthY-1] .= f0pS[2:gridlengthX-1, gridlengthY]          # bottom = top
        # diagonals streaming back to origin node
        fmmS[2:gridlengthX-1, gridlengthY-1] .= fppS[3:gridlengthX, gridlengthY]            # left-bottom = right-top
        fpmS[2:gridlengthX-1, gridlengthY-1] .= fmpS[1:gridlengthX-2, gridlengthY]          # right-bottom = left-top


        ## Bounceback cylinder
        fp0S[cylinder_indices], fm0S[cylinder_indices] = fm0S[cylinder_indices], fp0S[cylinder_indices] #horizontal getauscht
        f0pS[cylinder_indices], f0mS[cylinder_indices] = f0mS[cylinder_indices], f0pS[cylinder_indices] #vertikal getauscht
        fppS[cylinder_indices], fmmS[cylinder_indices] = fmmS[cylinder_indices], fppS[cylinder_indices] #diagonal getauscht rechtsoben <-> linksunten
        fpmS[cylinder_indices], fmpS[cylinder_indices] = fmpS[cylinder_indices], fpmS[cylinder_indices] #diagonal getauscht rechtsunten <-> linksoben
        
        ##### Boundary Conditions #####
        
        #Swap: SWAP POINTERS new distributions to array
        f00, f00S = f00S, f00
        fm0, fm0S = fm0S, fm0
        f0m, f0mS = f0mS, f0m
        fp0, fp0S = fp0S, fp0
        f0p, f0pS = f0pS, f0p
        fmm, fmmS = fmmS, fmm
        fmp, fmpS = fmpS, fmp
        fpp, fppS = fppS, fpp
        fpm, fpmS = fpmS, fpm


        # old BC's
        # ## Apply Boundary conditions
        # #Inlet velocity bc (unknown: f_1, f_8, f_9)
        # densityGrid[inlet, :] .= (sum(distributions[inlet, [1,3,5]], dims=2).+ 2 .*sum(distributions[inlet, [2,6,7]], dims=2)) ./ (1-lattice_inflow_velocity)
        # distributions[inlet, 4] .= distributions[inlet, 2] .+ (2/3 .* densityGrid[inlet,:] .* lattice_inflow_velocity)
        # distributions[inlet, 8] .= distributions[inlet, 6] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .- (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))
        # distributions[inlet, 9] .= distributions[inlet, 7] .+ (1/6 .* densityGrid[inlet,:] .* lattice_inflow_velocity) .+ (1/2 .* (distributions[inlet, 3] .- distributions[inlet, 5]))


        # #Outlet zero gradient bc
        # distributions[outlet, [4, 8, 9]] .= distributions[gridlengthX-1, :, [4, 8, 9]]


            # Plot of the field
            if ((i % 10 == 0)) || (i == simulationTime)

                # # #DEBUG print min/max of u and v
                # umax = maximum(u)
                # umin = minimum(u)
                # vmax = maximum(v)
                # vmin = minimum(v)
                # # println("Step $i — u min/max: $(round(umin, sigdigits=6)) / $(round(umax, sigdigits=6)) 
                # #             | v min/max: $(round(vmin, sigdigits=6)) / $(round(vmax, sigdigits=6))")

                # #DEBUG print min/max of u and v in physical units
                # vel_factor = delta_x / delta_t
                # umin_phys = umin * vel_factor
                # umax_phys = umax * vel_factor
                # vmin_phys = vmin * vel_factor
                # vmax_phys = vmax * vel_factor
                # # println("         u min/max (phys m/s): $(round(umin_phys, sigdigits=6)) / $(round(umax_phys, sigdigits=6)) 
                # #             | v min/max (phys m/s): $(round(vmin_phys, sigdigits=6)) / $(round(vmax_phys, sigdigits=6))")




                # Set velocities inside the cylinder to zero
                u_plot = copy(u)
                v_plot = copy(v)
                u_plot[cylinder] .= NaN
                v_plot[cylinder] .= NaN

                # Compute vorticity
                fill!(vorticity, 0.0)
                dv_dx = circshift(v_plot, (-1, 0)) .- circshift(v_plot, (1, 0))
                du_dy = circshift(u_plot, (0, -1)) .- circshift(u_plot, (0, 1))
                vorticity .= dv_dx .- du_dy
                vorticity[inlet] .= 0.0
                vorticity[outlet] .= 0.0

                # Mask the cylinder region
                vorticity[cylinder] .= NaN

                if ((i % 100 == 0)) || (i == simulationTime)
                    Log_Simulation_Runtime(i, simulationTime)
                end
                # Update the observables
                if Plotvorticity==true 
                    vorticity_obs[] = copy(vorticity) 
                    step_text[] = "Time step: $i, $(floor(Int, i*delta_t))s"
                end
                if Plotvx==true 
                    velocityX_obs[] = copy(u_plot) 
                    step_text_vx[] = "Time step: $i, $(floor(Int, i*delta_t))s"
                end
                if Plotvy==true 
                    velocityY_obs[] = copy(v_plot) 
                    step_text_vy[] = "Time step: $i, $(floor(Int, i*delta_t))s"
                end

                yield()
                sleep(0.05)
            end
    end

    Log_Simulation_Tail()
end
run_JuLattice()

