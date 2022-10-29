using FileIO, JLD2

include("BayesianVertex.jl")

function createFollowerModel(x,A,B,b)
    m = size(A,1);
    nx = size(A,2);
    ny = size(B,2);

    model = direct_model(Gurobi.Optimizer())
    set_silent(model)
    @variable(model, y[1:ny]);
    
    #next to lines are to force vectors to be interpreted matrices
    A = reshape(A,m,nx)
    B = reshape(B,m,ny)

    @constraint(model, B*y .<= b - A*x);

    return model, y
end


function samplevector(d)
    v = zeros(d,1) 

    while norm(v) < .0001 #avoid numerical issues
        for i=1:d
            v[i] = randn()
        end
    end
    v = v / norm(v)  # normalize to unit norm
end

function evalLeader(x, A, B, b, Fx, Fy, N)

    m = size(A,1);
    nx = size(A,2);
    ny = size(B,2);

    Y = zeros(ny,1)

    model,yvar = createFollowerModel(x,A,B,b)

    for i=1:N 
        c = samplevector(ny)
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
    ysol = JuMP.value.(yvar)
    return ysol
end



bolibinstances = [#"LiuHart1994",
#"GlackinEtal2009",
#"BialasKarwan1984a",
#"ClarkWesterberg1988",
#"LanWenShihLee2007",
#"Bard1984b",
"WangJiaoLi2005",
"MershaDempe2006Ex2",
"TuyEtal1994",
"BialasKarwan1984b",
"HaurieSavardWhite1990",
"HuHuangZhang2009",
"Bard1991Ex2",
"VisweswaranEtal1996",
"BardFalk1982Ex2",
"BenAyedBlair1990b",
"AnandalinghamWhite1990",
"BenAyedBlair1990a",
"ClarkWesterberg1990b",
"MershaDempe2006Ex1",
##"TuyEtal2007Ex3",
"Bard1984a",
"TuyEtal1993",
"CandlerTownsley1982"]

bolibdir = "BOLIBver2/JuliaExamples/"

#instancelist = ["ClarkWesterberg1988"]

coralinstances = [##"knapsack",
"linderoth",
##"milp_10_20_50_2310",
##"milp_4_20_10_0110",
"moore90_2"#
#"moore90"
]

coraldir = "CoralLib/notInterdiction/"


times = Dict()
vertex = Dict()

instancelist = bolibinstances
dir = bolibdir

for name in instancelist
    #print(@sprintf("Would read BOLIBver2/JuliaExamples/%s.jl\n",name))
    include(@sprintf("%s%s.jl",dir,name))
    couplingtolower = true
    (ncoupling,) = size(Gx)
    if couplingtolower && ncoupling >= 1
        local A = vcat(Gx,gx)
        local B = vcat(Gy,gy)
        local d = vcat(-bG,-bg)
    else
        local A = gx
        local B = gy
        local d = -bg
    end
    println("\nRunning ", name,"\n")
    #@show size(A)
    times[name] = @elapsed local X = AllVertex(A,B,d,true,true)
    #FileIO.save(@sprintf("%s_chamvtx.jld2",name),"X",X)
    vertex[name] = X

    print("\nVertices are:\n")
    @show vertex
    nvertex = size(X,2);
    N = 2
    print("\nComputing leader values in points\n")
    for i =1:nvertex
        x = X[:,i]
        leadervalue = evalLeader(x, A, B, d, Fx, Fy, N)
        println(leadervalue)
    end

end

#@show times
