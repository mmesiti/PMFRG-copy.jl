
function generate_compute_intensive_addX!(system)
    generate_compute_intensive_addX!(system.Npairs, system.Nsum, system.siteSum, system.NUnique)
end

function generate_compute_intensive_special_addX!(system)
    generate_compute_intensive_special_addX!(system.Npairs, system.Nsum, system.siteSum)
end

"Generates programmatically the compute-intensive code from a lattice specification"
function generate_compute_intensive_addX!(Npairs, Nsum, siteSum, NUnique)
    S_ki = siteSum.ki
    S_kj = siteSum.kj
    S_xk = siteSum.xk
    S_m = siteSum.m

    exp = quote
        function compute_intensive_addX!(a::Array{T,4},
            b::Array{T,4},
            c::Array{T,4},
            Va12::Vector{T},
            Vb12::Vector{T},
            Vc12::Vector{T},
            Va34::Vector{T},
            Vb34::Vector{T},
            Vc34::Vector{T},
            Vc21::Vector{T},
            Vc43::Vector{T},
            Props::SMatrix{$NUnique,$NUnique,T},
            is::Integer,
            it::Integer,
            iu::Integer) where {T}
        end
    end
    # Some ugly expr-fu t
    # o get to the list of expressions
    # inside "compute_intensive"
    func = exp.args[2]
    fbody = func.args[2].args

    for Rij = 1:Npairs
        #loop over all left hand side inequivalent pairs Rij

        push!(fbody, :(Xa_sum = 0.0))
        push!(fbody, :(Xb_sum = 0.0))
        push!(fbody, :(Xc_sum = 0.0))
        for k_spl = 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort!
            ki = S_ki[k_spl, Rij]
            kj = S_kj[k_spl, Rij]
            m = S_m[k_spl, Rij]
            xk = S_xk[k_spl, Rij]

            push!(fbody, :(@inbounds Ptm = Props[$xk, $xk] * $m))

            push!(fbody, :(@inbounds Xa_sum += (+Va12[$ki] * Va34[$kj] + Vb12[$ki] * Vb34[$kj] * 2) * Ptm))

            push!(fbody, :(@inbounds Xb_sum += (+Va12[$ki] * Vb34[$kj] + Vb12[$ki] * Va34[$kj] + Vb12[$ki] * Vb34[$kj]) * Ptm))

            push!(fbody, :(@inbounds Xc_sum += (+Vc12[$ki] * Vc34[$kj] + Vc21[$ki] * Vc43[$kj]) * Ptm))
        end
        push!(fbody, :(@inbounds a[$Rij, is, it, iu] += Xa_sum))
        push!(fbody, :(@inbounds b[$Rij, is, it, iu] += Xb_sum))
        push!(fbody, :(@inbounds c[$Rij, is, it, iu] += Xc_sum))
    end
    return exp
end

function generate_compute_intensive_special_addX!(Npairs, Nsum, siteSum)
    S_ki = siteSum.ki
    S_kj = siteSum.kj

    S_m = siteSum.m

    exp = quote
        function compute_intensive_special_addX!(a::Array{Float64,4},
            b::Array{Float64,4},
            c::Array{Float64,4},
            Va12::Vector{Float64},
            Vb12::Vector{Float64},
            Vc12::Vector{Float64},
            Va34::Vector{Float64},
            Vb34::Vector{Float64},
            Vc34::Vector{Float64},
            Vc21::Vector{Float64},
            Vc43::Vector{Float64},
            Props::SMatrix{1,1,Float64,1},
            is::Int64,
            it::Int64,
            iu::Int64)
        end
    end
    # Some ugly expr-fu t
    # o get to the list of expressions
    # inside "compute_intensive"
    func = exp.args[2]
    fbody = func.args[2].args

    for Rij = 1:Npairs
        #loop over all left hand side inequivalent pairs Rij

        push!(fbody, :(Xa_sum = 0.0))
        push!(fbody, :(Xb_sum = 0.0))
        push!(fbody, :(Xc_sum = 0.0))
        push!(fbody, :(Prop = only(Props)))
        for k_spl = 1:Nsum[Rij]
            #loop over all Nsum summation elements defined in geometry. This inner loop is responsible for most of the computational effort!
            ki = S_ki[k_spl, Rij]
            kj = S_kj[k_spl, Rij]
            m = S_m[k_spl, Rij]

            push!(fbody, :(@inbounds Xa_sum += (+Va12[$ki] * Va34[$kj] + Vb12[$ki] * Vb34[$kj] * 2) * $m))
            push!(fbody, :(@inbounds Xb_sum += (+Va12[$ki] * Vb34[$kj] + Vb12[$ki] * Va34[$kj] + Vb12[$ki] * Vb34[$kj]) * $m))
            push!(fbody, :(@inbounds Xc_sum += (+Vc12[$ki] * Vc34[$kj] + Vc21[$ki] * Vc43[$kj]) * $m))
        end
        push!(fbody, :(@inbounds a[$Rij, is, it, iu] += Xa_sum * Prop))
        push!(fbody, :(@inbounds b[$Rij, is, it, iu] += Xb_sum * Prop))
        push!(fbody, :(@inbounds c[$Rij, is, it, iu] += Xc_sum * Prop))
    end
    return exp
end
