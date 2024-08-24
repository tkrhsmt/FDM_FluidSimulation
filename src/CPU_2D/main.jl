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

using .Param
using .Navier
using .Poisson

# ---------- setting constant ---------- #
# input filename
INPUT_FILE = pwd() * "/" * ARGS[1]

# input parameter
prm = inputfile_param(INPUT_FILE)

# ---------- setting vaiable ---------- #
# old velocity (x and y direction)
ux1 = zeros(prm.nx + 3, prm.ny + 2)
uy1 = zeros(prm.nx + 2, prm.ny + 3)
# new velocity (x and y direction)
ux2 = zeros(prm.nx + 3, prm.ny + 2)
uy2 = zeros(prm.nx + 2, prm.ny + 3)
# old pressure
pp1 = zeros(prm.nx + 2, prm.ny + 2)
# new pressure
pp2 = zeros(prm.nx + 2, prm.ny + 2)
# velocity divergence
div = zeros(prm.nx + 2, prm.ny + 2)

# ---------- setting folder ---------- #
# if there is no "data" folfer, create it
if isdir("data") == false
    mkdir("data")
end

# ---------- main program ---------- #

function MAIN_PROGRAM(prm)

    # initial condition
    init!(ux1, uy1, prm)

    for time in prm.istart:prm.iend

        # 1st FS step (ux1, uy1) -> (ux2, uy2)
        first_velocity!(ux1, uy1, ux2, uy2, prm)

        # poisson equation pp1 -> pp1
        poisson!(ux2, uy2, pp1, pp2, div, prm)

        # 2nd FS step (ux2, uy2) -> (ux1, uy1)
        second_velocity!(ux2, uy2, ux1, uy1, pp1, prm)

        # output VTK file
        if time % prm.ioutput == 0
            filename = pwd() * "/data/output-$(Int(time / prm.ioutput))"
            output_vtkfile(ux1, uy1, pp1, prm, filename)
        end

        println("$time output")

    end


end

MAIN_PROGRAM(prm)
