using FileIO, JLD2

include("BayesianVertex.jl")

function checkBounded(A,B,b,bigM) #Check if polyhedron is bounded, otherwise add bigM as bound

    m = size(A,1);
    nx = size(A,2);
    ny = size(B,2);

    xbox = Inf*ones(nx,2)

    model = direct_model(Gurobi.Optimizer())
    set_silent(model)
    @variable(model, x[1:nx]);
    @variable(model, y[1:ny]);
    
    #the next to lines are to force vectors to be interpreted matrices
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

function createFollowerModel(x,A,B,b) # Create follower problem for a given leader vector x
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

function evalLeader(x, A, B, b, Fx, Fy, N, csample) # evaluate the leader objective using the given samples of the follower

    m = size(A,1);
    nx = size(A,2);
    ny = size(B,2);

    Y = zeros(ny,1)

    model,yvar = createFollowerModel(x,A,B,b)
    set_silent(model)
    
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

function solveFollower(c,model,yvar) # Solve the follower problem
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

# Variables mostly used for reports
timeexact = 0
timestoch = 0
valexact = 0
valstoch = 0

dim = 0
ncons = 0

exactsupportsize = 0
stochsupportsize = 0
totalfaces = 0
stochtotalfaces = 0

errstoch = 0

timefacelattice = 0

vertex = 0

couplingtolower = true
if ARGS[1] == "knapsack"
    name = "ContinuousKnapsack_instance.jl"
    couplingtolower = false
else
    name = ARGS[1]
end

println("\n========\n Running ", name, "\n=========\n")

include(@sprintf("%s",name))

(ncoupling,) = size(Gx)
if couplingtolower && ncoupling >= 1
    A = vcat(Gx,gx)
    B = vcat(Gy,gy)
    b = vcat(-bG,-bg)
else
    A = gx
    B = gy
    b = -bg
end

A,B,b,xbox = checkBounded(A,B,b,1000) #We add artificial bounds if we need to
nx = size(A,2);
ny = size(B,2);
m = size(A,1);

dim = nx+ny;
ncons = m;

#### Common computations ####
print("Running big M computation... ")
bigMarray = computeBigMs(A,B,b)
print("Done\n")

print("\nRunning faces computation... ")
start = time()
(K,s, K_MC, s_MC) = getFaces_pmk(A,B,b) #Warning: this function is implicitly assuming that all inequalities are facets    
timefacelattice = time() - start
print("Done\n")
#continue
######## Sampling of follower cost ###########

N = 100
csample = zeros(N, ny)
for i=1:N 
    csample[i,:] = samplevector_unitsphere(ny)
end 

runtimeexact = parse(Int,ARGS[2])

totalfaces = s
stochtotalfaces = s_MC

#### Exact method ####
runExact = true
minval = Inf
if runExact
    start = time()
    local (X,used_faces) = AllVertex(A, B, b,true, K, s, bigMarray, runtimeexact)
    vertex = X

    exactsupportsize = used_faces

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
            global minval = leadervalue
        end
    end

    timeexact = time() - start
    valexact = minval
end
#############################

#### Stochastic method ####
runStoch = true
Nsamplesx = 100
bestx = zeros(nx,1)
bestxval = Inf
z_all = zeros(s_MC) #To count the faces used
if runStoch
    start = time()
    
    seenlabels = Set()
    
    for i=1:Nsamplesx
        x_fix = samplevector_box(nx, xbox)

        found, zvals = findLabels(A, B, b, K_MC, s_MC, bigMarray, x_fix)

        if found
            global z_all += zvals
            Ind = findall(a->a>=0.5, vec(zvals));
            push!(seenlabels,Set(Ind))

            d = zeros(nx,1)
            success, t = findSteps(A, B, b, K_MC, s_MC, x_fix, zvals)
            if !success 
                continue
            end

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
                global bestxval = trueval
                global bestx = xopt
            end
        end
    end

    used_faces = count(a->a>=0.5, vec(z_all))

    timestoch = time() - start
    valstoch = bestxval

    stochsupportsize = used_faces

    ## Error computation
    seencount = 0
    for i=1:Nsamplesx
        x_fix = samplevector_box(nx, xbox)
        found, zvals = findLabels(A, B, b, K_MC, s_MC, bigMarray, x_fix)

        if found 
            Ind = findall(a->a>=0.5, vec(zvals));
            if Set(Ind) in seenlabels
                #println("Seen!")
                global seencount += 1
            end
        end
    end
    errstoch = 1- seencount/Nsamplesx
end


println("\n========\n Summary \n=========\n")

println("SUMMARY: ", name, " ValExact ", valexact, " ValStoch ", valstoch, " Gap ", (valexact - valstoch)/valexact, " ErrStoch ", errstoch, " TimeExact ", timeexact, " TimeStoch ", timestoch, " FaceLat-Time ", timefacelattice, " Dim ", dim, " NCons ", ncons," ExactSupp ", exactsupportsize, " StochSupp ", stochsupportsize, " NFaces ", totalfaces, " NFacesStoch ", stochtotalfaces)
