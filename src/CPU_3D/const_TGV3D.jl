# --------------------------------------------------
# FluidSimulation_CPU constant file
# julia code
# --------------------------------------------------

# number of lattices
NX = 64
NY = 64
NZ = 64

# length of flow field
LX = 2π
LY = 2π
LZ = 2π

# start and end number
ISTART = 1
IEND = 10000

# time span
DT = 0.0025

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
BZ = 0

# LES or DNS
LES = true
# LES model (1: smagorinsky, 2:coherent structure model)
LES_MODEL = 2

# boundary condition ------------------------------

function boundary_ux(ux, nx1, nx2, ny1, ny2, nz1, nz2)

    ux[1, :, :] = ux[end-2, :, :]
    ux[2, :, :] = ux[end-1, :, :]
    ux[end, :, :] = ux[3, :, :]
    ux[:, 1, :] = ux[:, end-1, :]
    ux[:, end, :] = ux[:, 2, :]
    ux[:, :, 1] = ux[:, :, end-1]
    ux[:, :, end] = ux[:, :, 2]

end

function boundary_uy(uy, nx1, nx2, ny1, ny2, nz1, nz2)

    uy[1, :, :] = uy[end-1, :, :]
    uy[end, :, :] = uy[2, :, :]
    uy[:, 1, :] = uy[:, end-2, :]
    uy[:, 2, :] = uy[:, end-1, :]
    uy[:, end, :] = uy[:, 3, :]
    uy[:, :, 1] = uy[:, :, end-1]
    uy[:, :, end] = uy[:, :, 2]

end

function boundary_uz(uz, nx1, nx2, ny1, ny2, nz1, nz2)

    uz[1, :, :] = uz[end-1, :, :]
    uz[end, :, :] = uz[2, :, :]
    uz[:, 1, :] = uz[:, end-1, :]
    uz[:, end, :] = uz[:, 2, :]
    uz[:, :, 1] = uz[:, :, end-2]
    uz[:, :, 2] = uz[:, :, end-1]
    uz[:, :, end] = uz[:, :, 3]

end

# init condition ------------------------------

function init!(ux1, uy1, uz1, prm)

    for k in 1:prm.nz+1
        for j in 1:prm.ny+1
            for i in 1:prm.nx+2
                ux1[i, j, k] = rand() * 0.01
                uy1[i, j, k] = rand() * 0.01
                uz1[i, j, k] = rand() * 0.01
            end
        end
    end
end

# velocity force ------------------------------

function force_ux(ux, nx1, nx2, ny1, ny2, nz1, nz2, dt)

    dx = 2π / (nx2 - nx1 + 1)
    dy = 2π / (ny2 - ny1 + 1)

    for k in nz1:nz2
        for j in ny1:ny2
            y = j * dy
            for i in nx1:nx2
                x = i * dx
                fx = -sin(x)*cos(y)
                ux[i, j, k] -= fx * dt
            end
        end
    end

end

function force_uy(uy, nx1, nx2, ny1, ny2, nz1, nz2, dt)

    dx = 2π / (nx2 - nx1 + 1)
    dy = 2π / (ny2 - ny1 + 1)

    for k in nz1:nz2
        for j in ny1:ny2
            y = j * dy
            for i in nx1:nx2
                x = i * dx
                fy = cos(x)*sin(y)
                uy[i, j, k] -= fy * dt
            end
        end
    end

end

function force_uz(uz, nx1, nx2, ny1, ny2, nz1, nz2, dt)

end
