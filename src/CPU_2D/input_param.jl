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

    lx :: Float64
    ly :: Float64

    istart :: Int64
    iend :: Int64

    dt :: Float64

    ilog :: Int64
    ioutput :: Int64

    re :: Float64

    poisson_scheme :: Int64

    bx :: Int64
    by :: Int64

end

function output_vtkfile(ux, uy, pp, fx, fy, prm, filename)

    dx = prm.lx / prm.nx
    dy = prm.ly / prm.ny

    x = dx:dx:dx*prm.nx
    y = dy:dy:dy*prm.ny

    # output vtk file
    vtk_grid(filename, x, y) do vtk
        vtk["ux"] = 0.5*(ux[2:prm.nx+1, 2:prm.ny+1] + ux[3:prm.nx+2, 2:prm.ny+1])
        vtk["uy"] = 0.5*(uy[2:prm.nx+1, 2:prm.ny+1] + uy[2:prm.nx+1, 3:prm.ny+2])
        vtk["pp"] = pp[2:prm.nx+1, 2:prm.ny+1]
    end

    # output checkout file
    @time h5open("checkpoint", "w") do file
        write(file,"ux",ux)
        write(file,"uy",uy)
        write(file,"fx",fx)
        write(file,"fy",fy)
    end

end

function inputfile_param(INPUT_FILE)

    # include parameter file
    include(INPUT_FILE)
    # setting parameter struct
    prm = parameter(NX, NY, LX, LY, ISTART, IEND, DT, ILOG, IOUTPUT, RE, POISSON_SCHEME, BX, BY)

    return prm, boundary_ux, boundary_uy, boundary_pp, init!, force_ux, force_uy
end

function inputdata()

    file = h5open("checkpoint", "r")
    ux = read(file,"ux")
    uy = read(file,"uy")
    fx = read(file,"fx")
    fy = read(file,"fy")

    close(file)
    return ux, uy, fx, fy
end

end
