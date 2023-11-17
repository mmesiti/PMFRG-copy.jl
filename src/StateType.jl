

"""Struct containing information about the (physical) ODE State, i.e. vertices"""
struct StateType{T}
    s::ArrayPartition{T}
    #f_int::Array{T,1}
    #γ::Array{T,2}
    #Γ::VertexType{T}
end

#############################
# Methods needed by solve() #
#############################
function Base.length(state::StateType)
    length(state.s)
end

function Base.iterate(state::StateType, args...)
    iterate(state.s, args...)
end

function Base.zero(state::StateType, args...)
    zero(state.s,args...) # TODO: parallelize
end

function Base.similar(state::StateType, args...)
    StateType(similar(state.s,args...))
end

function RecursiveArrayTools.recursivefill!(state::StateType, args...)
    println("args to RecursiveFill!: $(typeof(state.s)), $(args)")
    for arr in state.s.x
        fill!(arr,args...)
    end
    state
end

##################################
# end of methods needed by solve #
##################################

function StateType(NUnique::Int, Ngamma::Int, VDims::Tuple, type = Float64::Type)
    f = zeros(type, NUnique) # fint
    γ = zeros(type, NUnique, Ngamma) # gamma
    Γ = VertexType(VDims,type)

    s = ArrayPartition(f, γ, Γ.a, Γ.b, Γ.c)
    return StateType(s)
end

StateType(Par::PMFRGParams) = StateType(
    Par.System.NUnique,
    Par.NumericalParams.Ngamma,
    getVDims(Par),
    _getFloatType(Par),
)

function StateType(f_int::Array{T,1},
          γ::Array{T,2},
          Γa::Array{T,4},
          Γb::Array{T,4},
          Γc::Array{T,4}) where T
    StateType(ArrayPartition(f_int, γ, Γa, Γb, Γc))
end

RecursiveArrayTools.ArrayPartition(x::StateType) = x.s
#StateType(Arr::ArrayPartition) = StateType(Arr.x...)
