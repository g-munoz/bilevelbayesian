using FileIO, JLD2

include("BayesianVertex.jl")

function checkBounded(A,B,b,bigM)

    m = size(A,1);
    nx = size(A,2);
    ny = size(B,2);

    xbox = Inf*ones(nx,2)

    model = direct_model(Gurobi.Optimizer())
    set_silent(model)
    @variable(model, x[1:nx]);
    @variable(model, y[1:ny]);
    
    #next to lines are to force vectors to be interpreted matrices
    A = reshape(A,m,nx)
    B = reshape(B,m,ny)

    @constraint(model, A*x + B*y .<=b);
    
    for i = 1:nx
        #Upper bound
        @objective(model, Max, x[i]);
        optimize!(model)
        status = termination_status(model)
        if status != OPTIMAL
            newrowA = zeros(1,nx)
            newrowA[i] = 1
            newrowB = zeros(1,ny)

            A = vcat(A,newrowA)
            B = vcat(B,newrowB)
            b = vcat(b,bigM)

            xbox[i,2] = bigM
        else
            xbox[i,2] = objective_value(model)
        end
        
        #Lower bound
        @objective(model, Min, x[i]);
        optimize!(model)
        status = termination_status(model)
        if status != OPTIMAL
            newrowA = zeros(1,nx)
            newrowA[i] = -1
            newrowB = zeros(1,ny)

            A = vcat(A,newrowA)
            B = vcat(B,newrowB)
            b = vcat(b,bigM)

            xbox[i,1] = -bigM
        else
            xbox[i,1] = objective_value(model)
        end
    end

    for i = 1:ny
        #Upper bound
        @objective(model, Max, y[i]);
        optimize!(model)
        status = termination_status(model)
        if status != OPTIMAL
            newrowA = zeros(1,nx)
            newrowB = zeros(1,ny)
            newrowB[i] = 1

            A = vcat(A,newrowA)
            B = vcat(B,newrowB)
            b = vcat(b,bigM)
        end
        #Lower bound
        @objective(model, Min, y[i]);
        optimize!(model)
        status = termination_status(model)
        if status != OPTIMAL
            newrowA = zeros(1,nx)
            newrowB = zeros(1,ny)
            newrowB[i] = -1

            A = vcat(A,newrowA)
            B = vcat(B,newrowB)
            b = vcat(b,bigM)
        end
    end

    return A,B,b,xbox
end


function createFollowerModel(x,A,B,b)
    m = size(A,1);
    nx = size(A,2);
    ny = size(B,2);

    model = direct_model(Gurobi.Optimizer())
    set_silent(model)
    set_optimizer_attribute(model, "DualReductions", 0)

    @variable(model, y[1:ny]);
    
    #next to lines are to force vectors to be interpreted matrices
    A = reshape(A,m,nx)
    B = reshape(B,m,ny)

    @constraint(model, B*y .<= b - A*x);

    return model, y
end


function samplevector_unitsphere(d)
    v = zeros(d,1) 

    while norm(v) < .0001 #avoid numerical issues
        for i=1:d
            v[i] = randn()
        end
    end
    v = v / norm(v)  # normalize to unit norm
    return v
end

function samplevector_box(nx, xbox)
    v = zeros(nx,1) 

    v = xbox[:,1] + (xbox[:,2] - xbox[:,1]).*rand(nx)
    return v
end

function evalLeader(x, A, B, b, Fx, Fy, N, csample)

    m = size(A,1);
    nx = size(A,2);
    ny = size(B,2);

    Y = zeros(ny,1)

    model,yvar = createFollowerModel(x,A,B,b)

    for i=1:N 
        c = csample[i]
        ysol = solveFollower(c,model,yvar)
        Y = hcat(Y,ysol);
    end
    Y = Y[:,begin+1:end]
    if ny == 1
        avgY = (1/N)*sum(Y)
    else
        avgY = (1/N)*sum(Y,dims=2)
    end
    leaderobj = dot(Fx,x) + dot(Fy,avgY)
    return leaderobj
end

function solveFollower(c,model,yvar)
    ny = length(yvar)
    @objective(model, Min, sum(c.*yvar));
    optimize!(model)
    status = termination_status(model)
    if status != OPTIMAL
        println(status)
    end
    ysol = JuMP.value.(yvar)
    return ysol
end



bolibinstances = [
"LiuHart1994",
"GlackinEtal2009",
"BialasKarwan1984a",
"ClarkWesterberg1988",
"LanWenShihLee2007",
"Bard1984b",
"WangJiaoLi2005", #Unbounded
"MershaDempe2006Ex2",#Unbounded
"TuyEtal1994", #Unbounded #issue with t
"BialasKarwan1984b",
"HaurieSavardWhite1990",
"HuHuangZhang2009",
"Bard1991Ex2",
"VisweswaranEtal1996",
"BardFalk1982Ex2", #issue with t
"BenAyedBlair1990b",
"AnandalinghamWhite1990",
"BenAyedBlair1990a",
"ClarkWesterberg1990b",
"MershaDempe2006Ex1",
# #"TuyEtal2007Ex3", BIG
"Bard1984a",
"TuyEtal1993", #issue with t
"CandlerTownsley1982"] #issue with t

bolibdir = "BOLIBver2/JuliaExamples/"


coralinstances = [##"knapsack",
"linderoth",
##"milp_10_20_50_2310",
##"milp_4_20_10_0110",
"moore90_2",#
"moore90"
]

coraldir = "CoralLib/notInterdiction/"

timesexact = Dict()
timesstoch = Dict()
valsexact = Dict()
valsstoch = Dict()

dims = Dict()
ncons = Dict()

exactsupportsize = Dict()
stochsupportsize = Dict()
totalfaces = Dict()
stochtotalfaces = Dict()

errstoch = Dict()

timesfacelattice = Dict()

vertex = Dict()

instancelist = bolibinstances
dir = bolibdir

couplingtolower = true
if ARGS[1] == "knapsack"
    instancelist = ["ContinuousKnapsack_instance.jl"]
    couplingtolower = false
else
    instancelist = [ARGS[1]]
end

for name in instancelist
    println("\n========\n Running ", name, "\n=========\n")
    #print(@sprintf("Would read BOLIBver2/JuliaExamples/%s.jl\n",name))
    #include(@sprintf("%s%s.jl",dir,name))

    
    # knapflag = false
    # if name == "knapsack"
    #     name = "ContinuousKnapsack_instance.jl"
    #     knapflag = true
    # end
    include(@sprintf("%s",name))
    
    #if !knapflag
        (ncoupling,) = size(Gx)
        if couplingtolower && ncoupling >= 1
            local A = vcat(Gx,gx)
            local B = vcat(Gy,gy)
            local b = vcat(-bG,-bg)
        else
            local A = gx
            local B = gy
            local b = -bg
        end
    #else
    #     n = parse(Int, ARGS[3]);

    #     a = rand(1:4,n);
    #     #a = zeros(n);
    #     #for i=1:n
    #     #    a[i] = i%4 + 1
    #     #end

    #     A = zeros(2*n+3,1);
    #     A[1] = 1;
    #     A[2] = -1;
    #     A[3] = -1;

    #     B = vcat(zeros(3,n),Matrix(1.0I, n, n),Matrix(-1.0I, n, n));

    #     B[3,1:n] = a;

    #     ub = sum(a);
    #     lb = 0;

    #     b = vcat( [ub; -lb; 0], ones(n,1), zeros(n,1));

    #     #objective is Max sum(y) - b, then we make it min
    #     Fx = [1]
    #     Fy = -ones(n,1)

    #     Gx = []
    #     Gy = []
    # end
    
    A,B,b,xbox = checkBounded(A,B,b,1000) #We add artificial bounds if we need to
    nx = size(A,2);
    ny = size(B,2);
    m = size(A,1);

    dims[name] = nx+ny;
    ncons[name] = m;
    

    #### Common computations ####
    print("Running big M computation... ")
    bigMarray = computeBigMs(A,B,b)
    print("Done\n")

    print("\nRunning faces computation... ")
    start = time()
    (K,s, K_MC, s_MC) = getFaces_pmk(A,B,b) #Warning: this function is implicitly assuming that all inequalities are facets    
    timesfacelattice[name] = time() - start
    print("Done\n")
    #continue
    ######## Sampling of follower cost ###########

    N = 100
    csample = zeros(N, ny)
    for i=1:N 
        csample[i,:] = samplevector_unitsphere(ny)
    end 

    runtimeexact = parse(Int,ARGS[2])

    #### Exact method ####
    runExact = true
    minval = Inf
    if runExact
        start = time()
        local (X,used_faces) = AllVertex(A, B, b,true, K, s, bigMarray, runtimeexact)
        vertex[name] = X

	exactsupportsize[name] = used_faces
	totalfaces[name] = s
    stochtotalfaces[name] = s_MC

        print("\nVertices are:\n")
        @show vertex
        print("\nBox:\n")
        @show xbox
        nvertex = size(X,2);
        
        print("\nComputing leader values in points\n")
        
        for i =1:nvertex
            x = X[:,i]
            leadervalue = evalLeader(x, A, B, b, Fx, Fy, N,csample)
            #println(leadervalue)
            if leadervalue < minval
                minval = leadervalue
            end
        end

        timesexact[name] = time() - start
        valsexact[name] = minval
    end
    #############################

    #### Stochastic method ####
    runStoch = true
    Nsamplesx = 100
    if runStoch
        start = time()
        bestx = zeros(nx,1)
        bestxval = Inf

        seenlabels = Set()
        
        z_all = zeros(s_MC) #To count the faces used

        for i=1:Nsamplesx
            x_fix = samplevector_box(nx, xbox)
            #x_fix = [0.5; 2.0]
            found, zvals = findLabels(A, B, b, K_MC, s_MC, bigMarray, x_fix)

            #print(size(z_all), " ", size(zvals))
            

            if found

                z_all += zvals
                Ind = findall(a->a>=0.5, vec(zvals));
                push!(seenlabels,Set(Ind))

                d = zeros(nx,1)
                success, t = findSteps(A, B, b, K_MC, s_MC, x_fix, zvals)
                if !success #something strange happened in steps. Happening in CandlerTownsley1982
                    continue
                end
                # @show t
                # @show x_fix
                # @show found
                # @show zvals

                ej = zeros(nx,1)
                val2 = evalLeader(x_fix, A, B, b, Fx, Fy, N,csample)

                ## new attempt ####
                linA = zeros(nx+1,nx+1)
                linrhs = zeros(nx+1)
                ###################

                for j=1:nx
                    ej[j] = 1
                    val1 = evalLeader(x_fix+t[j]*ej, A, B, b, Fx, Fy, N,csample) 
                    d[j] = (val1-val2)/t[j]
                    
                    #############
                    # linA[j, :] = transpose(vcat(x_fix + t[j]*ej,1))
                    # linrhs[j] = val1
                    #############

                    ej[j] = 0
                end
                

                #############
                # linA[nx+1, :] = transpose(vcat(x_fix,1))
                # linrhs[nx+1] = val2
                # newd = linA\linrhs

                # xoptnew, optvalnew = solveLeader(A, B, b, K, s, newd[1:nx], zvals)
                # optvalnew += newd[nx+1]
                #############

                # @show d
                # @show newd[1:nx]

                xopt, optval = solveLeader(A, B, b, K_MC, s_MC, d, zvals)

                ## sanity check ##
                # @show x_fix
                # println("D problem gave value ", optval)
                # @show xopt
                # println("Eval leader on opt gave ", evalLeader(xopt, A, B, b, Fx, Fy, N))
                
                # println("new D problem gave value ", optvalnew)
                # @show xoptnew
                # println("Eval leader on opt gave ", evalLeader(xoptnew, A, B, b, Fx, Fy, N))
                ##################

                trueval = evalLeader(xopt, A, B, b, Fx, Fy, N,csample)
                if trueval < bestxval
                    bestxval = trueval
                    bestx = xopt
                end
            end
        end

        used_faces = count(a->a>=0.5, vec(z_all))

        timesstoch[name] = time() - start
        valsstoch[name] = bestxval

	    stochsupportsize[name] = used_faces

        ## Error computation
        seencount = 0
        for i=1:Nsamplesx
            x_fix = samplevector_box(nx, xbox)
            found, zvals = findLabels(A, B, b, K_MC, s_MC, bigMarray, x_fix)

            if found 
                Ind = findall(a->a>=0.5, vec(zvals));
                if Set(Ind) in seenlabels
                    println("Seen!")
                    seencount += 1
                end
            end
        end
        errstoch[name] = 1- seencount/Nsamplesx

    end

    # println("========================================")
    # println("========================================")
    # println("Optimal exact value for leader ", minval)
    # println("Optimal stoch value for leader ", bestxval)
    # println("========================================")
    # println("========================================")

    #break
end

println("\n========\n Summary \n=========\n")
for name in instancelist
    println("SUMMARY: ", name, " Vals ", valsexact[name], " ", valsstoch[name], " Gap ", (valsexact[name] - valsstoch[name])/valsexact[name], " Err ", errstoch[name], " Times ", timesexact[name], " ", timesstoch[name], " FL time ", timesfacelattice[name], " Dim ", dims[name], " Cons ", ncons[name]," Exact Supp ", exactsupportsize[name], " Stoch Supp ", stochsupportsize[name], " NFaces ", totalfaces[name], " NFacesStoch ", stochtotalfaces[name] )
end
