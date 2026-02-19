module JuLattice
    ############################
    ## Main file for JuLattice #
    ############################
    include("src/Plotter.jl")
    include("src/Logger.jl")

    using MeshGrid, GLMakie
    using .Plotter, .Logger

    function run()
        ## User Settings
        # Domain Settings
        lengthX = 8;              # m
        lengthY = 2;              # m 

        # Rectangle Settings (we will expand on this later -> only for Re right now)
        d = 0.5;                   # m

        # Fluid Properties
        Re = 300
        machNumber = 0.1;          # Target Mach number (Ma = U / c_s)
        viscosity = 0.001;          # m²/s

        # Simulation Settings
        simulationTime = 8000;     # s
        deltaX = 0.1;              # Grid spacing (physical units per lattice unit)


        # Plot Requests
        plotU = true;
        plotV = true;
        plotVorticity = true;

        #### Run Simulation #####
        Log_Simulation_Header()

        ##-------- Compute LBM Parameters from Fluid Properties & Simulation Settings --------##
        speedOfSound = viscosity * Re / (machNumber * d);
        latticeSpeedOfSound = 1 / √3;

        deltaT = deltaX * latticeSpeedOfSound / speedOfSound
        latticeViscosity = viscosity * deltaT / (deltaX^2)
        latticeInflowVelocity = machNumber * latticeSpeedOfSound
        ## Convert user settings to lattice units
        # Domain

        # Fluid
        latticeDensity = 1;

        #Log 
        Log_Discretization_Settings(deltaX, deltaT, latticeDensity)

        # Simulation Settings
        simulationTime = ceil(Int, simulationTime / deltaT);
        Q   = 9;

        velocityVector = [     [0, -1, 0, 1, 0, -1, -1, 1, 1],
                                [0, 0, 1, 0, -1, -1, 1, 1, -1]];

        weights =   [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

        omegaBGK = 1 / (3 * latticeViscosity + 0.5);
        omegaAcoustic = 1;
        omegaGhosts = 1;

        relaxationVector = [omegaBGK, omegaAcoustic, omegaGhosts]; #should be 0 for conserved moments (00,01,10), omegaBGK for hydrodynamic moments (11,20-02), 
                                                                   #omegaAcoustic for acoustic moments (20+02) and 1 for ghost moments (12,21,22)

        # Now the main arrays. Initialissation should always be Nx x Ny (where Nx e.g. is LengthX/deltaX but always +2 (to create a "boundary node" around the actual domain for easier implementation of boundary conditions))
        # 9 separate arrays for distributions following miller indices (f00, f01, f10, f11, f20, f02, f12, f21, f22)
        # 9 + 9 + 9 separate arrays for precollision moments (m00, m01, m10, m11, m20-m02, m02+m02, m12, m21, m22), post collision (m00s, m01s, m10s, m11s, m20s-m02s, m20s+m02s, m12s, m21s, m22s) and their equilibrium values (m00_eq, m01_eq, m10_eq, m11_eq, m20_eq, m02_eq, m12_eq, m21_eq, m22_eq)
        
        # initialise transformation matrix (M) and its inverse (M_inv) to transform distribution to moments and back 

        # create boundary indetifiers
        walls = 
        inlet = 
        outlet = 

        # create fluid identifiers (only nodes which are actual fluid (not inlet, walls or solid objects))

        # Initialize macroscopic density and scale distribution
        densityGrid = sum(distributions, dims=3);
        distributions .*= fluiddensity ./ densityGrid; # or something similar

        # Initialize macroscopic velocity arrays
        velocityX   = zeros()
        velocityY   = zeros()

        # Plotting setup -> adjust such that it works for the new setup. All in one window (3 optional subplots)
        if any((Plotvorticity, Plotvx, Plotvy))
            if Plotvorticity==true 
                vorticity, vorticity_obs, text_obj, step_text, fig_vorticity = Create_Plot(gridlengthX, gridlengthY)
                screen1 = GLMakie.Screen()
                GLMakie.display(screen1, fig_vorticity)
            end
            if Plotvx==true 
                velocityX_obs, text_obj_vx, step_text_vx, fig_vx = Create_Plot(gridlengthX, gridlengthY, velocityX, "X")
                screen2 = GLMakie.Screen(; position = (600, 0))
                GLMakie.display(screen2, fig_vx)
            end
            if Plotvy==true 
                velocityY_obs, text_obj_vy, step_text_vy, fig_vy = Create_Plot(gridlengthX, gridlengthY, velocityY, "Y")
                screen3 = GLMakie.Screen()
                display(screen3, fig_vy)
            end
        end

        println("#################################")
        println("Starting Simulation:")
        # Run Simulation Loop
        for i in 1:simulationTime
            # Only for fluid nodes!!:
                # Transform from f to m

                # Recover macroscopic variables (M00 density, M10 velcotityX/density, M02 velocityY/density)

                # Relax correctly (with M02+M20 - omegaAcoustic and M02-M20 - omegaBGK)
                # Relaxing means: m_s = m - omega .* (m - m_eq) where m is the precollision moment, m_eq the equilibrium moment and m_s the post collision moment. Omega is the relaxation parameter which can be different for different moments.

                # Transform back from m to f

                # Streaming step (shift distributions according to their velocity vector)(f00 stays f0p shifts to the right f0m shifts to the left etc.)
            
            # Apply BC -> Inlet: Velocity BounceBack: Since we are "post streaming" we have to be mindful where our population has streamed to, when computing new fs.
            # Velocity bounceBack means: the post streaming westwards population is now inside the inlet node. We return it to the fluid node as eastwards population (f0m becomes f0p) same for the post streaming northwest and southwest (these are in the inlet node x+1 Y-1 and x-1 Y+1 respectively). 
            # We return them with an added Momentum according to latticeInflowVelocity with the equation: +2/(latticeSpeedOfSound*latticeSpeedOfSound) * weighti (velocityVectori * latticeInflowVelocity)

            # -> Outlet: we do zero gradient (Neumann) BC -> we just copy the post streaming values from the last fluid node (x-1) to the outlet node (x). This means we copy fp0, fpp, fpm from the last fluid node to the outlet node.
            
            # -> Walls: equilibrium BC. So we set all wall boundary nodes to their equilibrium distribution with latticeInflowVelocity.

                # Plot of the field -> Adjust to new logic described above
                if ((i % 10 == 0)) || (i == simulationTime)
                    # Set velocities inside the cylinder to zero
                    velocityX[cylinder] .= NaN
                    velocityY[cylinder] .= NaN

                    # Compute vorticity
                    fill!(vorticity, 0.0)
                    dv_dx = circshift(velocityY, (-1, 0)) .- circshift(velocityY, (1, 0))
                    du_dy = circshift(velocityX, (0, -1)) .- circshift(velocityX, (0, 1))
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
                        velocityX_obs[] = copy(velocityX) 
                        step_text_vx[] = "Time step: $i, $(floor(Int, i*delta_t))s"
                    end
                    if Plotvy==true 
                        velocityY_obs[] = copy(velocityY) 
                        step_text_vy[] = "Time step: $i, $(floor(Int, i*delta_t))s"
                    end

                    yield()
                    sleep(0.05)
                end

        end

        Log_Simulation_Tail()
    end#run

end#JuLattice