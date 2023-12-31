function test_SDE(
    Par::PMFRG.ParquetParams = Params(
        getPolymer(2),
        Parquet(),
        T = 1.0,
        N = 30,
        Ngamma = 30,
        accuracy = 1e-5,
        Lam_min = 0.0,
        Lam_max = 100.0,
        usesymmetry = false,
        lenIntw = 800,
        lenIntw_acc = 1000,
        maxIterSDE = 300,
        maxIterBSE = 80,
        MinimalOutput = true,
    );
    tol = 1e-5,
)
    Workspace = PMFRG.SetupParquet(Par)
    Lam = 0.0

    (; State, Γ0, X, B0, BX, Par, Buffer) = Workspace
    @inline Prop(x, nw) = 1 / 6 * PMFRG.iG_(State.γ, x, Lam, nw, Par.NumericalParams.T)

    test_FirstBSEIteration!(Workspace, Lam)
    gamma1 = copy(State.γ)
    gamma2 = copy(State.γ)


    PMFRG.compute1PartBubble!(gamma1, B0, Prop, Par)
    PMFRG.compute1PartBubble_BS!(gamma2, State.Γ, Γ0, Prop, Par)
    @testset "SDE gamma" begin
        for i in eachindex(gamma1, gamma2)
            @test gamma1[i] ≈ gamma2[i] atol = tol
        end
    end
end

function test_FirstBSEIteration!(Workspace, Lam)
    (; OldState, State, Γ0, X, B0, BX, Par, Buffer) = Workspace

    getProp! = PMFRG.constructPropagatorFunction(Workspace, Lam)

    PMFRG.PMFRG.computeLeft2PartBubble!(B0, Γ0, Γ0, State.Γ, getProp!, Par, Buffer)
end

function test_SDE_FP(Par::PMFRG.ParquetParams)
    Workspace1 = PMFRG.SetupParquet(Par)
    Workspace2 = PMFRG.SetupParquet(Par)
    Lam = 0.0
    test_FirstBSEIteration!(Workspace1, Lam)
    test_FirstBSEIteration!(Workspace2, Lam)

    PMFRG.iterateSDE!(Workspace1, Lam)
    # PMFRG.iterateSDE_FP!(Workspace2,Lam)
    @testset "SDE FixedPoint" begin
        State1 = Workspace1.State
        State2 = Workspace2.State
        for i in eachindex(State1.γ, State2.γ)
            @test State1.γ[i] ≈ State2.γ[i]
        end
        for field in (:a, :b, :c)
            V1 = getproperty(State1.Γ, field)
            V2 = getproperty(State2.Γ, field)
            for i in eachindex(V1, V2)
                @test V1[i] ≈ V2[i]
            end
        end

    end
end


function test_iterate!(Workspace::PMFRG.ParquetWorkspace, Lam::Real; SDECheatfactor = 1)
    (; OldState, State, I, Γ0, X, B0, BX, Par, Buffer) = Workspace

    maxIterBSE = Par.Options.maxIterBSE


    Tol_Vertex = 1E16 #Initial tolerance level

    getProp! = constructPropagatorFunction(Workspace, Lam)
    iter = 0


    PMFRG.writeTo!(OldState.Γ, State.Γ)

    # PMFRG.iterateSDE!(Workspace,Lam)

    PMFRG.computeLeft2PartBubble!(B0, Γ0, Γ0, State.Γ, getProp!, Par, Buffer)
    PMFRG.computeLeft2PartBubble!(BX, X, X, State.Γ, getProp!, Par, Buffer)

    PMFRG.getXFromBubbles!(X, B0, BX) #TODO when applicable, this needs to be generalized for beyond-Parquet approximations 
    PMFRG.getVertexFromChannels!(State.Γ, I, X)
    PMFRG.symmetrizeVertex!(State.Γ, Par)

    # PMFRG.iterateSDE!(Workspace,Lam,SDECheatfactor = SDECheatfactor)
    # PMFRG.iterateSDE_FP!(Workspace,Lam)
    PMFRG.iterateSDE_FP!(State.γ, State.Γ, Γ0, Lam, Par)
    Tol_Vertex = reldist(OldState.Γ, State.Γ)

    isnan(Tol_Vertex) && @warn ": BSE: Solution diverged after $iter iterations\n"
    println(maximum.((PMFRG.ArrayPartition(OldState), PMFRG.ArrayPartition(State))))
    return Workspace
end

function test_iterate_FP!(Workspace::PMFRG.ParquetWorkspace, Lam::Real)
    (; OldState, State, I, Γ0, X, B0, BX, Par, Buffer) = Workspace

    maxIterBSE = Par.Options.maxIterBSE

    OldStateArr, StateArr = PMFRG.ArrayPartition.((OldState, State))

    function FixedPointFunction!(State_Arr, OldState_Arr)
        anyisnan(OldState_Arr) && return State_Arr
        gamma = OldState_Arr.x[2]

        State = StateType(State_Arr.x...)
        OldState = StateType(OldState_Arr.x...)
        getProp! = constructPropagatorFunction(gamma, Lam, Par)

        PMFRG.computeLeft2PartBubble!(B0, Γ0, Γ0, OldState.Γ, getProp!, Par, Buffer)
        PMFRG.computeLeft2PartBubble!(BX, X, X, OldState.Γ, getProp!, Par, Buffer)

        PMFRG.getXFromBubbles!(X, B0, BX) #TODO when applicable, this needs to be generalized for beyond-Parquet approximations 
        PMFRG.getVertexFromChannels!(State.Γ, I, X)
        PMFRG.symmetrizeVertex!(State.Γ, Par)


        WS_State = PMFRG.ArrayPartition(Workspace.State)
        WS_State .= State_Arr
        # PMFRG.iterateSDE_FP!(Workspace,Lam)
        PMFRG.iterateSDE_FP!(State.γ, State.Γ, Γ0, Lam, Par)

        println(maximum.((OldState_Arr, State_Arr)))
        return State_Arr
    end

    # FixedPointFunction!(StateArr,OldStateArr)
    s = afps!(FixedPointFunction!, OldStateArr, iters = 1, vel = 0.0, ep = 1.0, tol = 0.0)
    # StateArr .= s.x
    # OldStateArr .= StateArr
    # isnan(Tol_Vertex) && @warn ": BSE: Solution diverged after $iter iterations\n"
    return Workspace
end

function test_BSE(Par::PMFRG.ParquetParams)
    Workspace1 = PMFRG.SetupParquet(Par)
    Workspace2 = PMFRG.SetupParquet(Par)
    Lam = 0.0

    for _ = 1:4
        test_iterate!(Workspace1, Lam)
        test_iterate_FP!(Workspace2, Lam)
        println("FP!:")
    end
    Gamma1 = Workspace1.State.Γ
    Gamma2 = Workspace2.State.Γ
    # @testset "BSE test" begin
    #     for k in fieldnames(VertexType)
    #         G1,G2 = getfield.((Gamma1,Gamma2),Ref(k))
    #         for i in eachindex(G1,G2)
    #             @test G1[i] ≈ G2[i]
    #         end
    #     end
    # end
    println(maximum(Gamma1.c))
    println(maximum(Gamma2.c))
    println(maximum(Workspace1.State.γ))
    println(maximum(Workspace2.State.γ))
end
