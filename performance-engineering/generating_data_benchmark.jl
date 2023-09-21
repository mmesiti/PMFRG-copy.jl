using SpinFRGLattices, PMFRG
using SpinFRGLattices.SquareLattice
using StaticArrays


# Number of nearest neighbor bonds
# up to which correlations are treated in the lattice.
# For NLen = 5, all correlations C_{ij} are zero
#if sites i and j are separated by more than 5 nearest neighbor bonds.
NLenToy = 5
NLen = 14
J1 = 1
J2 = 0.1
# Construct a vector of couplings:
# nearest neighbor coupling is J1 (J2)
# and further couplings to zero.
# For finite further couplings simply provide a longer array,
# i.e [J1,J2,J3,...]
couplings = [J1, J2]

# create a structure that contains all information about the geometry of the problem.

SystemToy = getSquareLattice(NLenToy, couplings)
@time "generate" compute_intensive_toy_expr = PMFRG.generate_compute_intensive(SystemToy)
@time "eval" compute_intensive_toy = eval(compute_intensive_toy_expr)
function trigger_compile_compute_intensive_toy()
    Npairs = SystemToy.Npairs
    NUnique = SystemToy.NUnique
    a = zeros(Float64, (Npairs, 1, 1, 1))# ::Array{T,4},
    b = zeros(Float64, (Npairs, 1, 1, 1))# ::Array{T,4},
    c = zeros(Float64, (Npairs, 1, 1, 1))# ::Array{T,4},
    Va12 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vb12 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vc12 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Va34 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vb34 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vc34 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vc21 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vc43 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Props = SMatrix{NUnique,NUnique,Float64}(zeros(Float64, (NUnique, NUnique)))
    is = 1 #::Integer,
    it = 1 #::Integer,
    iu = 1 #::Integer
    compute_intensive_toy(a,
        b,
        c,
        Va12,
        Vb12,
        Vc12,
        Va34,
        Vb34,
        Vc34,
        Vc21,
        Vc43,
        Props,
        is,
        it,
        iu)
end

@time "compilation" trigger_compile_compute_intensive_toy()

println("Warm up")

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    SystemToy, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    T=0.5, # Temperature for the simulation.
    N=10, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy=1e-3, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    MinimalOutput=true,
)

tempdir = "temp"
println("Removing data from previous runs ($tempdir)")
rm(tempdir, recursive=true, force=true)
mainFile = "$tempdir/" * PMFRG.generateFileName(Par, "_testFile") # specify a file name for main Output
flowpath = "$tempdir/flows/" # specify path for vertex checkpoints

Solution, saved_values = SolveFRG(
    compute_intensive_toy,
    Par,
    MainFile=mainFile,
    CheckpointDirectory=flowpath,
    method=DP5(),
    VertexCheckpoints=[],
    CheckPointSteps=3,
);



println("Warmup done, timing real problem now.")
System = getSquareLattice(NLen, couplings)
@time "generate" compute_intensive_expr = PMFRG.generate_compute_intensive(System)
@time "eval" compute_intensive = eval(compute_intensive_expr)
function trigger_compile_compute_intensive()
    Npairs = System.Npairs
    NUnique = System.NUnique
    a = zeros(Float64, (Npairs, 1, 1, 1))# ::Array{T,4},
    b = zeros(Float64, (Npairs, 1, 1, 1))# ::Array{T,4},
    c = zeros(Float64, (Npairs, 1, 1, 1))# ::Array{T,4},
    Va12 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vb12 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vc12 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Va34 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vb34 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vc34 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vc21 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Vc43 = zeros(Float64, (Npairs,)) # ::Vector{T},
    Props = SMatrix{NUnique,NUnique,Float64}(zeros(Float64, (NUnique, NUnique)))
    is = 1 #::Integer,
    it = 1 #::Integer,
    iu = 1 #::Integer
    compute_intensive(a,
        b,
        c,
        Va12,
        Vb12,
        Vc12,
        Va34,
        Vb34,
        Vc34,
        Vc21,
        Vc43,
        Props,
        is,
        it,
        iu)
end

@time "compilation" trigger_compile_compute_intensive()




Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    T=0.5, # Temperature for the simulation.
    N=25, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy=1e-3, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    MinimalOutput=true,
)

tempdir = "temp"
println("Removing data from previous runs ($tempdir)")
rm(tempdir, recursive=true, force=true)
mainFile = "$tempdir/" * PMFRG.generateFileName(Par, "_testFile") # specify a file name for main Output
flowpath = "$tempdir/flows/" # specify path for vertex checkpoints


@time Solution, saved_values = SolveFRG(
    compute_intensive,
    Par,
    MainFile=mainFile,
    CheckpointDirectory=flowpath,
    method=DP5(),
    VertexCheckpoints=[],
    CheckPointSteps=3,
);

println("Done")
