# --------------------------------------------------
# FluidSimulation_CPU constant file
# julia code
# --------------------------------------------------

# number of lattices
NX = 64
NY = 64
NZ = 64

# length of flow field
LX = 1.0
LY = 1.0
LZ = 1.0

# start and end number
ISTART = 1
IEND = 10000

# time span
DT = 0.001

# output log file
ILOG = 10

# output vtk file
IOUTPUT = 100

# reynolds number
RE = 1000

# poisson solver (0 : FFT spectral method, 1 : multigrid method)
POISSON_SCHEME = 0

# pressure boundary (0 : periodic, 1 : neumann)
BX = 1
BY = 1
BZ = 1

# boundary condition ------------------------------

function boundary_ux(ux, nx1, nx2, ny1, ny2, nz1, nz2)

    ux[1:2, :, :] .= 0.0
    ux[end-1:end, :, :] .= 0.0
    ux[:, 1, :] = - ux[:, 2, :] .+ 2.0
    ux[:, end, :] = - ux[:, end-1, :]
    ux[:, :, 1] = - ux[:, :, 2]
    ux[:, :, end] = - ux[:, :, end-1]

end

function boundary_uy(uy, nx1, nx2, ny1, ny2, nz1, nz2)

    uy[:, 1:2, :] .= 0.0
    uy[:, end-1:end, :] .= 0.0
    uy[1, :, :] = - uy[2, :, :]
    uy[end, :, :] = - uy[end-1, :, :]
    uy[:, :, 1] = - uy[:, :, 2]
    uy[:, :, end] = - uy[:, :, end-1]

end

function boundary_uz(uz, nx1, nx2, ny1, ny2, nz1, nz2)

    uz[1, :, :] = - uz[2, :, :]
    uz[end, :, :] = - uz[end-1, :, :]
    uz[:, 1, :] = - uz[:, 2, :]
    uz[:, end, :] = - uz[:, end-1, :]
    uz[:, :, 1:2] .= 0.0
    uz[:, :, end-1:end] .= 0.0

end

# init condition ------------------------------

function init!(ux1, uy1, uz1, prm)

    #ux1 .= 1.0
end

# velocity force ------------------------------

function force_ux(ux, nx1, nx2, ny1, ny2, nz1, nz2, dt)

end

function force_uy(uy, nx1, nx2, ny1, ny2, nz1, nz2, dt)

end

function force_uz(uz, nx1, nx2, ny1, ny2, nz1, nz2, dt)

end
