using HDF5 
using DelimitedFiles

export GetOrbUpDn

export WriteVecMat
export WriteVec
export WriteMat
export ReadMat

export ReadHDF5
export WriteHDF5

export FlattenArr

export WriteStatus
export WriteStatusItem

function GetOrbUpDn( nspinorb )
    return ( [ i for i in 1:2:nspinorb ] , [ i for i in 2:2:nspinorb ] )
end

function WriteVecMat( w_grid, M_grid, fname ) 
    open(fname, "w" ) do io
        writedlm( io, [ real(w_grid) real(M_grid) imag(M_grid) ] )
    end
end

function ReadMat( fname  ; dtype=ComplexF64 )
    M   = readedlm( fname )
    return M[:,1:2:end] + M[:,2:2:end]*im
end

function WriteMat( M, fname ) 
    open(fname, "w" ) do io
        writedlm( io, [ real(M) imag(M) ] )
    end
end

function WriteVec( vec, fname ) 
    open(fname, "w" ) do io
        writedlm( io, vec )
    end
end

function ReadHDF5( fname , obj_fullpath )
    fid = h5open( fname, "r" ) 
    val = read( fid, obj_fullpath )
    close(fid)
    return val
end

function ReadHDF5( fname , obj_fullpath, T::DataType )
    val = ReadHDF5(fname, obj_fullpath)
    if (T==Vector{Matrix{ComplexF64}}) | (T==Vector{Matrix{Float64}})
        return ConvertArr3DToVecMat(val)
    elseif (T==Vector{Vector{ComplexF64}}) | (T==Vector{Vector{Float64}})
        return ConvertArr2DToVecVec(val)
    elseif (T==Vector{Vector{Vector{ComplexF64}}}) | (T==Vector{Vector{Vector{Float64}}})
        return ConvertArr3DToVecVecVec(val)
    else
        return val
    end
end

function ConvertVecMatToArr3D( VM::Vector{Matrix{T}} ) where T
    dimV, (dim1,dim2) = (size(VM)[1], size(VM[1]))
    res = zero( Array{T,3}(undef,dimV,dim1,dim2) )
    for iv=1:dimV, i1=1:dim1, i2=1:dim2
        res[iv,i1,i2]   = VM[iv][i1,i2]
    end
    return res
end

function ConvertVecVecVecToArr3D( VM::Vector{Vector{Vector{T}}} ) where T
    dimV, dim1, dim2 = ( size(VM)[1], size(VM[1])[1], size(VM[1][1])[1] )
    res = zero( Array{T,3}(undef,dimV,dim1,dim2) )
    for iv=1:dimV, i1=1:dim1, i2=1:dim2
        res[iv,i1,i2]   = VM[iv][i1][i2]
    end
    return res
end

function ConvertVecVecToArr2D( VM::Vector{Vector{T}} ) where T <: Number
    dimV, dim1= (size(VM)[1], size(VM[1])[1])
    @show (dimV, dim1)
    @show T
    res = zero( Array{T,2}(undef,dimV,dim1) )
    for iv=1:dimV, i1=1:dim1
        res[iv,i1]   = VM[iv][i1]
    end
    return res
end

function ConvertArr3DToVecMat( A3::Array{T,3} ) where T
    dimV, dim1,dim2 = size(A3)
    res = [ A3[iv,:,:] for iv=1:dimV ]
    return res
end

function ConvertArr3DToVecVecVec( A3::Array{T,3} ) where T
    dimV, dim1,dim2 = size(A3)
    res = [ [ A3[iv,i1,:] for i1=1:dim1 ] for iv=1:dimV ]
    return res
end

function ConvertArr2DToVecVec( A2::Array{T,2} ) where T
    dimV, dim1 = size(A2)
    res = [ A2[iv,:] for iv=1:dimV ]
    return res
end

function WriteHDF5( fname , obj_fullpath, val ; mode="cw" )
    fid = h5open( fname, mode ) 
    if haskey( fid, obj_fullpath )
        delete_object( fid, obj_fullpath )
    end
    fid[obj_fullpath] = val
    close(fid)
end

function WriteHDF5( fname , obj_fullpath, val::Vector{Matrix{T}} ; mode="cw" ) where T
    val_formatted   = ConvertVecMatToArr3D(val)
    fid = h5open( fname, mode ) 
    # @show typeof(val_formatted)
    if haskey( fid, obj_fullpath )
        delete_object( fid, obj_fullpath )
    end
    fid[obj_fullpath] = val_formatted
    close(fid)
end

function WriteHDF5( fname , obj_fullpath, val::Vector{Vector{T}} ; mode="cw" ) where T <: Number
    val_formatted   = ConvertVecVecToArr2D(val)
    fid = h5open( fname, mode ) 
    # @show typeof(val_formatted)
    if haskey( fid, obj_fullpath )
        delete_object( fid, obj_fullpath )
    end
    fid[obj_fullpath] = val_formatted
    close(fid)
end

function WriteHDF5( fname , obj_fullpath, val::Vector{Vector{Vector{T}}} ; mode="cw" ) where T <: Number
    val_formatted   = ConvertVecVecVecToArr3D(val)
    fid = h5open( fname, mode ) 
    # @show typeof(val_formatted)
    if haskey( fid, obj_fullpath )
        delete_object( fid, obj_fullpath )
    end
    fid[obj_fullpath] = val_formatted
    close(fid)
end

function WriteHDF5_CFData( fname , obj_fullpath, CFD ; mode="cw" )
    fid = h5open( fname, mode ) 
    if haskey( fid, obj_fullpath )
        delete_object( fid, obj_fullpath )
    end
    fid[obj_fullpath*"/AbsAVec"]     = CFD[1]
    fid[obj_fullpath*"/TriMatDiag"]  = CFD[2]
    fid[obj_fullpath*"/TriMatUpper"] = CFD[3]
    close(fid)
end

function ReadHDF5_CFData( fname , obj_fullpath, CFD )
    fid = h5open( fname, "r" ) 
    AAV = read( fid, obj_fullpath*"/AbsAVec" )
    TMD = read( fid, obj_fullpath*"/TriMatDiag" )
    TMU = read( fid, obj_fullpath*"/TriMatUpper" )
    close(fid)
    return AAV, TMD, TMU
end

function FlattenArr( arr ) 
    return collect( Iterators.flatten( arr ) ) 
end

function WriteStatus( fname, data_arr )
    open( fname, "a" ) do io
        writedlm( io, data_arr )
    end
end

function WriteStatusItem( fname, data )
    open( fname, "a" ) do io
        println( io, data )
    end
end

