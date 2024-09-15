# --------------------------------------------------
# MODULE    Param
# julia code
#
# DATE : 2024/8/23
# --------------------------------------------------

module Param

using WriteVTK
using HDF5

export inputfile_param, output_vtkfile, inputdata

struct parameter

    nx :: Int64
    ny :: Int64
    nz :: Int64

    lx :: Float64
    ly :: Float64
    lz :: Float64

    istart :: Int64
    iend :: Int64

    dt :: Float64

    ilog :: Int64
    ioutput :: Int64

    re :: Float64

    poisson_scheme :: Int64

    bx :: Int64
    by :: Int64
    bz :: Int64

    les :: Bool
    les_model :: Int64

end

function output_vtkfile(ux, uy, uz, pp, fx, fy, fz, prm, filename)

    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny
    dz = prm.lz / prm.nz

    x = dx:dx:dx*prm.nx
    y = dy:dy:dy*prm.ny
    z = dz:dz:dz*prm.nz

    # output vtk file
    vtk_grid(filename, x, y, z) do vtk
        vtk["ux"] = 0.5*(ux[2:prm.nx+1, 2:prm.ny+1, 2:prm.nz+1] + ux[3:prm.nx+2, 2:prm.ny+1, 2:prm.nz+1])
        vtk["uy"] = 0.5*(uy[2:prm.nx+1, 2:prm.ny+1, 2:prm.nz+1] + uy[2:prm.nx+1, 3:prm.ny+2, 2:prm.nz+1])
        vtk["uz"] = 0.5*(uz[2:prm.nx+1, 2:prm.ny+1, 2:prm.nz+1] + uz[2:prm.nx+1, 2:prm.ny+1, 3:prm.nz+2])
        vtk["pp"] = pp[2:prm.nx+1, 2:prm.ny+1, 2:prm.nz+1]
    end

    # output checkout file
    @time h5open("checkpoint", "w") do file
        write(file,"ux",ux)
        write(file,"uy",uy)
        write(file,"uz",uz)
        write(file,"fx",fx)
        write(file,"fy",fy)
        write(file,"fz",fz)
    end

end

function inputfile_param(INPUT_FILE)

    # include parameter file
    include(INPUT_FILE)
    # setting parameter struct
    prm = parameter(NX, NY, NZ, LX, LY, LZ, ISTART, IEND, DT, ILOG, IOUTPUT, RE, POISSON_SCHEME, BX, BY, BZ, LES, LES_MODEL)

    return prm, boundary_ux, boundary_uy, boundary_uz, init!, force_ux, force_uy, force_uz
end

function inputdata()

    file = h5open("checkpoint", "r")
    ux = read(file,"ux")
    uy = read(file,"uy")
    uz = read(file,"uz")
    fx = read(file,"fx")
    fy = read(file,"fy")
    fz = read(file,"fz")

    close(file)
    return ux, uy, uz, fx, fy, fz
end

end
