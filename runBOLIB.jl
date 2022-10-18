using FileIO, JLD2

include("BayesianVertex.jl")

instancelist = ["LiuHart1994",
"GlackinEtal2009",
"BialasKarwan1984a",
"ClarkWesterberg1988",
"LanWenShihLee2007",
"Bard1984b",
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
#"TuyEtal2007Ex3",
"Bard1984a",
"TuyEtal1993",
"CandlerTownsley1982"]

#instancelist = ["ClarkWesterberg1988"]

times = Dict()
vertex = Dict()

for name in instancelist
    #print(@sprintf("Would read BOLIBver2/JuliaExamples/%s.jl\n",name))
    include(@sprintf("BOLIBver2/JuliaExamples/%s.jl",name))
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
end

@show times
@show vertex