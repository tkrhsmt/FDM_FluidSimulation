# --------------------------------------------------
# FluidSimulation_CPU constant file
# julia code
# --------------------------------------------------

# number of lattices
NX = 128
NY = 64
NZ = 64

# length of flow field
LX = 2.0
LY = 0.5
LZ = 1.0

# start and end number
ISTART = 1
IEND = 60000

# time span
DT = 0.001

# output log file
ILOG = 10

# output vtk file
IOUTPUT = 200

# reynolds number
RE = 50000

# poisson solver (0 : FFT spectral method, 1 : multigrid method)
POISSON_SCHEME = 0

# pressure boundary (0 : periodic, 1 : neumann)
BX = 0
BY = 1
BZ = 0

# boundary condition ------------------------------

function boundary_ux(ux, nx1, nx2, ny1, ny2, nz1, nz2)

    ux[1, :, :] = ux[end-2, :, :]
    ux[2, :, :] = ux[end-1, :, :]
    ux[end, :, :] = ux[3, :, :]
    ux[:, :, 1] = ux[:, :, end-1]
    ux[:, :, end] = ux[:, :, 2]
    ux[:, 1, :] = - ux[:, 2, :] .- 2.0
    ux[:, end, :] = - ux[:, end-1, :] .+ 2.0

end

function boundary_uy(uy, nx1, nx2, ny1, ny2, nz1, nz2)

    uy[1, :, :] = uy[end-1, :, :]
    uy[end, :, :] = uy[2, :, :]
    uy[:, :, 1] = uy[:, :, end-1]
    uy[:, :, end] = uy[:, :, 2]
    uy[:, 1:2, :] .= 0.0
    uy[:, end-1:end, :] .= 0.0

end

function boundary_uz(uz, nx1, nx2, ny1, ny2, nz1, nz2)

    uz[1, :, :] = uz[end-1, :, :]
    uz[end, :, :] = uz[2, :, :]
    uz[:, :, 1] = uz[:, :, end-2]
    uz[:, :, 2] = uz[:, :, end-1]
    uz[:, :, end] = uz[:, :, 3]
    uz[:, 1, :] = - uz[:, 2, :]
    uz[:, end, :] = - uz[:, end-1, :]

end

# init condition ------------------------------

function init!(ux1, uy1, uz1, prm)

    dy = prm.ly / prm.ny

    for k in 1:prm.nz+1
        for j in 1:prm.ny+1
            y = dy * (j - 1)
            for i in 1:prm.nx+2
                ux1[i, j, k] = 4.0 * (y - 0.25) + rand() * 0.3
                uy1[i, j, k] = 0.0 + rand() * 0.3
                uz1[i, j, k] = 0.0 + rand() * 0.3
            end
        end
    end

end

# velocity force ------------------------------

function force_ux(ux, nx1, nx2, ny1, ny2, nz1, nz2, dt)

end

function force_uy(uy, nx1, nx2, ny1, ny2, nz1, nz2, dt)

end

function force_uz(uz, nx1, nx2, ny1, ny2, nz1, nz2, dt)

end
