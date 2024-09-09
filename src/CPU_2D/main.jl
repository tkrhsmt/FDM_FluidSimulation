# ==================================================
# FluidSimulation CPU
# julia code
#
# DATE : 2024/8/23
# ==================================================

# using package
include("input_param.jl")
include("navier.jl")
include("poisson.jl")
include("log.jl")

using .Param
using .Navier
using .Poisson
using .Log

# ---------- setting constant ---------- #
# input filename
INPUT_FILE = pwd() * "/" * ARGS[1]

# input parameter
prm, boundary_ux, boundary_uy, boundary_pp, init!, force_ux, force_uy = inputfile_param(INPUT_FILE)

# ---------- setting vaiable ---------- #
# old velocity (x and y direction)
ux1 = zeros(prm.nx + 3, prm.ny + 2)
uy1 = zeros(prm.nx + 2, prm.ny + 3)
# intermidiate velocity (x and y direction)
ux2 = zeros(prm.nx + 3, prm.ny + 2)
uy2 = zeros(prm.nx + 2, prm.ny + 3)
# new velocity (x and y direction)
ux3 = zeros(prm.nx + 3, prm.ny + 2)
uy3 = zeros(prm.nx + 2, prm.ny + 3)
# old pressure
pp1 = zeros(prm.nx + 2, prm.ny + 2)
# new pressure
pp2 = zeros(prm.nx + 2, prm.ny + 2)
# velocity divergence
div = zeros(prm.nx + 2, prm.ny + 2)
# adams-bashforth
fx1 = zeros(prm.nx + 3, prm.ny + 2)
fy1 = zeros(prm.nx + 2, prm.ny + 3)
fx2 = zeros(prm.nx + 3, prm.ny + 2)
fy2 = zeros(prm.nx + 2, prm.ny + 3)

# ---------- setting folder ---------- #
# if there is no "data" folfer, create it
if isdir("data") == false
    mkdir("data")
end

# --------- setting function --------- #
# select poisson solver
poisson_solver = poisson_fft
if prm.poisson_scheme == 1
    poisson_solver = poisson!
end

# ---------- main program ---------- #

function MAIN_PROGRAM(prm, ux1, uy1, ux2, uy2, ux3, uy3, pp1, pp2, div, fx1, fx2, fy1, fy2)

    # initial log
    init_log(prm)

    # if starting from the middle, read the checkpoint file
    if prm.istart == 1
        # initial condition
        init!(ux1, uy1, prm)

        # first 1 step is euler method

        # 1st FS step (ux1, uy1) -> (ux2, uy2)
        first_velocity!(ux1, uy1, ux2, uy2, prm, boundary_ux, boundary_uy, force_ux, force_uy)

        # poisson equation pp1 -> pp1
        #rep = poisson!(ux2, uy2, pp1, pp2, div, prm, boundary_pp)
        rep = poisson_solver(ux2, uy2, pp1, pp2, div, prm, boundary_pp)

        # 2nd FS step (ux2, uy2) -> (ux3, uy3)
        second_velocity!(ux2, uy2, ux3, uy3, pp1, prm, boundary_ux, boundary_uy)

        # preserve adams-bashforth data
        #fx2 = (ux3 - ux1) / prm.dt
        #fy2 = (uy3 - uy1) / prm.dt

        # return next step first data
        ux1 = copy(ux3)
        uy1 = copy(uy3)

    else
        # read the checkpoint file
        ux1, uy1, fx2, fy2 = inputdata()
    end


    # adams-bashforth after 1 step
    for time in prm.istart+1:prm.iend

        # 1st FS step (ux1, uy1) -> (ux2, uy2)
        first_velocity!(ux1, uy1, ux2, uy2, prm, boundary_ux, boundary_uy, force_ux, force_uy)

        # poisson equation pp1 -> pp1
        #rep = poisson!(ux2, uy2, pp1, pp2, div, prm, boundary_pp)
        rep = poisson_solver(ux2, uy2, pp1, pp2, div, prm, boundary_pp)

        # 2nd FS step (ux2, uy2) -> (ux3, uy3)
        second_velocity!(ux2, uy2, ux3, uy3, pp1, prm, boundary_ux, boundary_uy)

        # adams-bashforth data 1
        fx1 = (ux3 - ux1) / prm.dt
        fy1 = (uy3 - uy1) / prm.dt

        # return next step first data (adams-bashforth method)
        ux1 = ux1 + prm.dt * (3 / 2 * fx1 - 1 / 2 * fx2)
        uy1 = uy1 + prm.dt * (3 / 2 * fy1 - 1 / 2 * fy2)

        # adams-bashforth data 2
        fx2 = copy(fx1)
        fy2 = copy(fy1)

        # output log
        if time % prm.ilog == 0
            calculation_log(ux1, uy1, pp1, prm, rep, time)
        end

        # output VTK file
        if time % prm.ioutput == 0
            filename = pwd() * "/data/output-$(Int(time / prm.ioutput))"
            output_vtkfile(ux1, uy1, pp1, fx2, fy2, prm, filename)
            output_log(time, prm)
        end

    end

    final_log()


end

MAIN_PROGRAM(prm, ux1, uy1, ux2, uy2, ux3, uy3, pp1, pp2, div, fx1, fx2, fy1, fy2)
