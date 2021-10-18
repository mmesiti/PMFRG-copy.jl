"""
This file contains functions used for reading and writing FRG data to disc.
Stored data can be used to re-launch an incomplete FRG calculation.
"""

"""Saves Vertices to a compressed HDF5 file in a Group "Lam"."""
function saveState(Filename::String,State::ArrayPartition,Lam,mode = "cw")
    Vertices = State.x
    Names = "fint","gamma","Va","Vb","Vc"
    # Filename = string(DirName,"/$(string(round(Lam,digits =3))).h5")
    try
        h5open(Filename,mode) do f
            for (Name,V) in zip(Names,Vertices)
                f["$Name",blosc = 9] = V
            end
        end
        h5write(Filename,"Lam",Lam)
    catch e
        @warn "Saving state was unsuccessfull with exception:\n $e "
    end
end

"""Reads Vertices from file"""
function readState(Filename::String)
    Names = "fint","gamma","Va","Vb","Vc"
    State = ArrayPartition(
        (h5read(Filename,"$N") for N in Names)...
    )
    return State
end

function readLam(Filename::String)
    Lam = h5read(Filename,"Lam")
    return Lam
end

EssentialParamFields() = (
    :T,
    :N,
    :Ngamma,
    :accuracy,
    :Lam_min,
    :Lam_max,
    :ex_freq,
    :lenIntw,
    :lenIntw_acc
)

"""Saves important information about computation parameters so that they can be reconstructed"""
function saveParams(Filename,Par::Params)
    Fields = EssentialParamFields()
    for F in Fields
        h5write(Filename,"Params/$F",getfield(Par,F))
    end
    h5write(Filename,"Geometry/Name",Par.System.Name)
    h5write(Filename,"Geometry/NLen",Par.System.NLen)
    h5write(Filename,"Geometry/couplings",Par.System.couplings)
    h5write(Filename,"Geometry/Npairs",Par.System.Npairs)
end

function readParams(Filename,Geometry;modifyParams...)
    Fields = EssentialParamFields()
    Kwargs = Dict((F => h5read(Filename,"Params/$F") for F in Fields)...)
    Par = Params(;System = Geometry,Kwargs...,modifyParams...)
    return Par
end

function setupDirectory(DirPath,Par)
    DirPath = generateUniqueName(DirPath,Par)
    println("Checkpoints saved at $DirPath")
    # CheckPath = joinpath(DirPath,"Checkpoints")
    mkpath(DirPath)
    return DirPath 
end

function saveCurrentState(DirPath::String,State::AbstractArray,saved_Values::DiffEqCallbacks.SavedValues,Lam::Real,Par::Params)
    saveState(joinpath(DirPath,"CurrentState.h5"),State,Lam,"w")
    saveParams(joinpath(DirPath,"CurrentState.h5"),Par)
    saveObs(joinpath(DirPath,"CurrentState.h5"),saved_Values,"Observables")
end

function UniqueDirName(FullPath)
    newpath = FullPath
    versionPath(index) = string(FullPath,"(v_",index,")")
    while ispath(newpath)
        it = findfirst("(v_",newpath)
        if it !== nothing
            currindex = parse(Int,newpath[it[end]+1:end-1])
            newpath =versionPath(currindex+1)
        else
            newpath =versionPath(1)
        end
    end
    return newpath
end


function generateUniqueName(Directory::String,Par::Params)
    @unpack Name,T,N = Par
    Name = "$(Name)_N=$(N)_T=$T"
    Name = UniqueDirName(joinpath(Directory,Name))
    return Name
end

function setupFromCheckpoint(Filename::String,Geometry)
    State = readState(Filename)
    Lam = readLam(Filename)
    Par = readParams(Filename,Geometry;Lam_max = Lam)
    @unpack N,Ngamma,Npairs,VDims,couplings,T,NUnique = Par
    println("Reading Checkpoint from $Filename")
    # println("starting with ",generateName(Par))
    X = CreateX(3,VDims)
    XTilde = CreateX(4,VDims)
    return State,(X,XTilde,Par)
end

function SolveFRG_Checkpoint(Filename::String,Geometry;kwargs...)
    State,setup = setupFromCheckpoint(Filename,Geometry) #Package parameter and pre-allocate arrays
    Par = setup[3]
    FilePath = dirname(Filename)
    launchPMFRG!(State,setup,getDeriv!,Par;CheckpointDirectory = FilePath,kwargs...)
end

"""Saves Observables"""
function saveObs(Filename,Obs::DiffEqCallbacks.SavedValues,Group = "")
    Fields = fieldnames(eltype(Obs.saveval))
    ObsArr = StructArray(Obs.saveval)
    for F in Fields
        arr = convertToArray(getproperty(ObsArr,F))
        h5write(Filename,joinpath(Group,string(F)),arr)
    end
    h5write(Filename,joinpath(Group,"Lambda"),Obs.t)
end

function convertToArray(VecOfArray::AbstractVector{VT}) where {T,N,VT <: AbstractArray{T,N}}
    cat(VecOfArray...,dims = N+1)
end