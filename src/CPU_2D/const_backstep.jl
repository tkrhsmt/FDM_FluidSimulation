# --------------------------------------------------
# FluidSimulation_CPU constant file
# julia code
# --------------------------------------------------

# number of lattices
NX = 256
NY = 64

# length of flow field
LX = 2.0
LY = 0.5

# start and end number
ISTART = 1
IEND = 10000

# time span
DT = 0.002

# output log file
ILOG = 10

# output vtk file
IOUTPUT = 100

# reynolds number
RE = 2000

# poisson solver (0 : FFT spectral method, 1 : multigrid method)
POISSON_SCHEME = 0

# pressure boundary (0 : periodic, 1 : neumann)
BX = 1
BY = 1

# boundary condition ------------------------------

function boundary_ux(ux, nx1, nx2, ny1, ny2)

    #=
    ux[1:2, :] .= 0.0
    ux[end-1:end, :] .= 0.0
    ux[:, 1] = - ux[:, 2] .+ 2.0
    ux[:, end] = - ux[:, end-1]
    =#
    ux[1, :] .= 1.0
    ux[end, :] .= ux[end - 1, :]

    ux[:, 1] = - ux[:, 2]
    ux[:, end] = - ux[:, end-1]

    nx = nx2 - nx1 + 1
    ny = ny2 - ny1 + 1
    dx = 2.0 / nx
    dy = 0.5 / ny

    for j in ny1:ny2
        y = dy * (j - ny1)
        for i in nx1:nx2
            x = dx * (i - nx1)
            if x < 0.5 && y < 0.25
                ux[i, j] = 0.0
            end
        end
    end

end

function boundary_uy(uy, nx1, nx2, ny1, ny2)

    #=
    uy[:, 1:2] .= 0.0
    uy[:, end-1:end] .= 0.0
    uy[1, :] = - uy[2, :]
    uy[end, :] = - uy[end-1, :]
    =#
    uy[:, 1:2] .= 0.0
    uy[:, end-1:end] .= 0.0

    uy[1, :] = - uy[2, :]
    uy[end, :] = uy[end-1, :]

    nx = nx2 - nx1 + 1
    ny = ny2 - ny1 + 1
    dx = 2.0 / nx
    dy = 0.5 / ny

    for j in ny1:ny2
        y = dy * (j - ny1)
        for i in nx1:nx2
            x = dx * (i - nx1)
            if x < 0.5 && y < 0.25
                uy[i, j] = 0.0
            end
        end
    end

end

function boundary_pp(pp, nx1, nx2, ny1, ny2)

    pp[1, :] = pp[2, :]
    pp[end, :] = pp[end-1, :]
    pp[:, 1] = pp[:, 2]
    #pp[:, 1] .= 0.0
    pp[:, end] = pp[:, end-1]

    pav = sum(pp)
    pp = pp .- pav

end

# init condition ------------------------------

function init!(ux1, uy1, prm)

end

# velocity force ------------------------------

function force_ux(ux, nx1, nx2, ny1, ny2, dt)

end

function force_uy(uy, nx1, nx2, ny1, ny2, dt)

end
