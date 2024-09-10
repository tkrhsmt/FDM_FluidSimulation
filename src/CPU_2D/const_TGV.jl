# --------------------------------------------------
# FluidSimulation_CPU constant file
# julia code
# --------------------------------------------------

# number of lattices
NX = 128
NY = 128

# length of flow field
LX = 2π
LY = 2π

# start and end number
ISTART = 1
IEND = 10000

# time span
DT = 0.003

# output log file
ILOG = 10

# output vtk file
IOUTPUT = 100

# reynolds number
RE = 5000

# poisson solver (0 : FFT spectral method, 1 : multigrid method)
POISSON_SCHEME = 0

# pressure boundary (0 : periodic, 1 : neumann)
BX = 0
BY = 0

# boundary condition ------------------------------

function boundary_ux(ux, nx1, nx2, ny1, ny2)

    ux[1, :] = ux[end-2, :]
    ux[2, :] = ux[end-1, :]
    ux[end, :] = ux[3, :]

    ux[:, 1] = ux[:, end-1]
    ux[:, end] = ux[:, 2]

end

function boundary_uy(uy, nx1, nx2, ny1, ny2)

    uy[:, 1] = uy[:, end-2]
    uy[:, 2] = uy[:, end-1]
    uy[:, end] = uy[:, 3]

    uy[1, :] = uy[end-1, :]
    uy[end, :] = uy[2, :]

end

function boundary_pp(pp, nx1, nx2, ny1, ny2)

    pp[1, :] = pp[end-1, :]
    pp[end, :] = pp[2, :]
    pp[:, 1] = pp[:, end-1]
    pp[:, end] = pp[:, 2]

    pav = sum(pp)
    pp = pp .- pav

end

# init condition ------------------------------

function init!(ux1, uy1, prm)

    for j in 2:prm.ny+1
        y = j / prm.ny * 2π
        for i in 2:prm.nx+1
            x = i / prm.nx * 2π
            ux1[i, j] = -sin(10*x) * cos(10*y)
            uy1[i, j] = cos(10*x) * sin(10*y)
        end
    end

end

# velocity force ------------------------------

function force_ux(ux, nx1, nx2, ny1, ny2, dt)

end

function force_uy(uy, nx1, nx2, ny1, ny2, dt)

end
