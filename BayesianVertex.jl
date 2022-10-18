using LinearAlgebra, JuMP, Gurobi, Polymake, Printf


function computeBigMs(A,B,b)

    #(m,nx) = size(A);
    m = size(A,1);
    nx = size(A,2);
    ny = size(B,2);

    #(ny,) = size(B');

    bigMs = -1e6*ones(m)

    model = direct_model(Gurobi.Optimizer())
    set_silent(model)
    @variable(model, x[1:nx]);
    @variable(model, y[1:ny]);
    
    #next to lines are to force vectors to be interpreted matrices
    A = reshape(A,m,nx)
    B = reshape(B,m,ny)

    @constraint(model, A*x + B*y .<=b);

    #@constraint(model, [i = 1:m], A[i]*x + B[i]*y <=b[i])
    
    for i = 1:m
        @objective(model, Min, sum(A[i,j]*x[j] for j=1:nx) + sum( B[i,j]*y[j] for j=1:ny));
        optimize!(model)
        status = termination_status(model)
        #println(sum(A[i,j]*x[j] for j=1:nx) + sum( B[i,j]*y[j] for j=1:ny), " ", objective_value(model))
        if status == OPTIMAL
            bigMs[i] += (objective_value(model) - 1)
        else
            bigMs[i] = -1e6
        end
    end
    return bigMs
end

function verifysol(ver_model, zverify, z_val, outverify, out_val, xverify, yverify)

    s = length(zverify)
    for i = 1 : s
        fix(zverify[i], z_val[i]; force = true)
    end
    fix(outverify, out_val; force = true)

    optimize!(ver_model)
    status = termination_status(ver_model)

    if status == INFEASIBLE
        return false, nothing, nothing
    end

    xval = JuMP.value.(xverify)
    yval = JuMP.value.(yverify)

    #println("Opt status: ",status," Val: ", objective_value(ver_model))
    return true, xval, yval
end

function AllVertex(A,B,b,cbFlag,bigMflag)

    #(m,nx) = size(A,1);
    m = size(A,1);
    nx = size(A,2);
    ny = size(B,2);

    ###############################################
    #We recover all the faces of D = {Ax + By <= b}
    #of dimension nx. These are given by the nonempty
    #intersection of ny affine spaces with linearly
    #independent normal vectors.
    ###############################################

    bigMarray = -1e6*ones(m)
    if bigMflag
        bigMarray = computeBigMs(A,B,b)
    end

    (K,s) = getFaces_pmk(A,B,b,m,nx,ny) #Warning: this function is implicitly assuming that all inequalities are facets
    
    (m,zvar,xvar,outvar,yvar) = MIPVertexSearchBase(A,B,b,K,s,nx,ny,bigMarray);

    ### this copy is for a verifier model
    ### for technical reasons couldn't use copy method
    (ver_model,zvar_ver,xvar_ver,outvar_ver,yvar_ver) = MIPVertexSearchBase(A,B,b,K,s,nx,ny,bigMarray);
    set_silent(ver_model)

    #println("Tolerances:")
    #println(get_optimizer_attribute(m, "FeasibilityTol"))
    #println(get_optimizer_attribute(m, "IntFeasTol"))
    #set_optimizer_attribute(m, "FeasibilityTol", 1e-4)
    #set_optimizer_attribute(m, "IntFeasTol", 1e-4)

    function greedy_callback_GRB(cb_data,cb_where::Cint)

        if cb_where != GRB_CB_MIPSOL
            return
        end
        
        Gurobi.load_callback_variable_primal(cb_data, cb_where)
        s = length(zvar)
        z_val = callback_value.(Ref(cb_data), zvar)
        out_val = callback_value(cb_data, outvar)
        
        if out_val >= 1
            return
        end

        feasible, x_val, y_val = verifysol(ver_model, zvar_ver, z_val, outvar_ver, out_val, xvar_ver, yvar_ver)
        if !feasible
            print("Error: Gurobi provided an invalid solution")
            return
        end

        improved = false
        for i = 1 : s
            if z_val[i] >= 0.5
                continue
            end
            z_new = copy(round.(z_val))
            z_new[i] = 1
            #@show zvar
            #@show z_val
            #@show z_new

            #println("Submitting sol")
            status = MOI.submit(m, MOI.HeuristicSolution(cb_data), zvar, z_new)
            #println("Submitted a heuristic solution, and the status was: ", status)
            #sleep(5)
            if status == MOI.HEURISTIC_SOLUTION_ACCEPTED
                z_val = z_new
                improved = true
                println("Went in!")
                return
                #sleep(5)
            else
                feasible, x_val, y_val = verifysol(ver_model, zvar_ver, z_new, outvar_ver, out_val, xvar_ver, yvar_ver)
                if feasible
                    println("Warning: Gurobi rejected a valid better solution, not stopping")
                    z_val = z_new
                    improved = true
                    return
                    
                    println("Test if passing more info fixes it")
                    # @show size(zvar)
                    # @show size(z_val)
                    # @show size(xvar)
                    # @show size(x_val)
                    # @show size(yvar)
                    # @show size(y_val)
                    status = MOI.submit(m, MOI.HeuristicSolution(cb_data), vcat(zvar,xvar,vec(yvar),outvar), vcat(z_new,x_val,vec(y_val),out_val))
                    println("Submitted a heuristic solution, and the status was: ", status,"  out ", out_val)

                    jumpvars = vcat(zvar,xvar,vec(yvar),outvar)
                    jumpvals = vcat(z_new,x_val,vec(y_val),out_val)
                    dictionary = Dict(jumpvars .=> jumpvals)
                    println(primal_feasibility_report(m, dictionary))

                    #println(vcat(zvar,xvar,vec(yvar),outvar), vcat(z_new,x_val,vec(y_val),out_val))
                    #sleep(5)
                #else
                #    println("Well rejected solution")
                #    sleep(5)
                end
            end
        end
        if !improved
            println("Solution was not improved. Should be maximal")
            #sleep(5)
            GRBterminate(backend(m))
        else
            println("A solution was improved") #should terminate too
            #sleep(5)
        end

        return
    end


    if cbFlag
        MOI.set(m, Gurobi.CallbackFunction(), greedy_callback_GRB)
    end
    X = zeros(nx,1)

    iter = 1
    while true
        status = optimize!(m)
        z = JuMP.value.(zvar)
        x = JuMP.value.(xvar)
        out = JuMP.value.(outvar);
        value = objective_value(m)
        if out==1
            break
        end
        println("\nCollected sol ", iter,"\n")
        iter += 1
        X = hcat(X,x);

        #We force at least one constraint z[k] to be active in the complement
        #of z. If this is not possibe, out is activated
        Ind = findall(a->a==0, vec(z));
        @constraint(m,sum(zvar[Ind]) + outvar>=1);

        @constraint(ver_model,sum(zvar_ver[Ind]) + outvar_ver>=1);
        #@constraint(m, sum(zvar) <= value); #this may be a bad idea, since it is adding a ineq parallel to obj
        #write_to_file(m, "modeliter.lp")
    end
    X = X[:,begin+1:end]
  
    return X
end

######################################################################################

function getFaces_pmk(A,B,b,m,nx,ny)
    P = polytope.Polytope(INEQUALITIES=hcat(b,-A,-B)) ##Careful here on how Polymake takes the input
    
    HD = P.HASSE_DIAGRAM
    allfaces_str = @sprintf("%s", HD.FACES) ## String with all faces
    allfaces_vec = split(allfaces_str, "\n")
    allfaces_vec = allfaces_vec[begin+1:end-1] # The first and last lines are not faces
    
    faces_dim_nx = @convert_to Array{Int} graph.nodes_of_rank(HD,nx+1) #rank is one more than dimension    
    #@show faces_dim_nx
    #@show size(faces_dim_nx)

    for i = 1:nx
        faces_rank_i = @convert_to Array{Int} graph.nodes_of_rank(HD,i)
        faces_dim_nx = vcat(faces_dim_nx, faces_rank_i)
    end

    #@show faces_dim_nx
    #@show size(faces_dim_nx)
    #@show size(allfaces_vec)

    s = length(faces_dim_nx)
    
    K = zeros(s,m)
    
    for i = 1 : s
    	#@show allfaces_vec[i][begin+1:end-1]
    	face_string = allfaces_vec[faces_dim_nx[i]+1][begin+1:end-1] #index of the faces of dim x start at 0. Begin and end are to remove brackets and then be able to cast
    	
    	vertices_in_face = @convert_to Array{Int} face_string
    	vertices_in_face = vertices_in_face .+ 1 #again, the vertex indices start at 0 so we correct them
    	
    	ineqs = P.VERTICES_IN_INEQUALITIES[:,vertices_in_face] #boolean sparse matrix, row = inequality. Returns all inequalities tight at each vertex
    	
    	ind = zeros(1,m) #this will indicate the inequalities that are tight for ALL vertices
    	for row = 1 : m
    		#this piece of code should be done in one line, but the structure of polymake SparseBoolVec is strange
    		alltrue = true
    		for val in ineqs[row,:]
    			alltrue = alltrue && val 
    		end	
    		if alltrue
    			ind[row] = 1
    		end
    	end
    	K[i,:] = ind	

    end
    return K, s
end

######################################################################################

function MIPVertexSearchBase(A,B,b,K,s,nx,ny,bigMarray)

    ###########################################
    # nx -> X-dimension
    # ny -> Y-dimension
    # s  -> Number of nx-faces of D = {(x,y) : Ax+By <= b}
    ###########################################
    
        m = direct_model(Gurobi.Optimizer())
    
        @variable(m, x[1:nx]);
        @variable(m, y[1:ny,1:s]);
        @variable(m, 0<=z[1:s]<=1, Int);
        @variable(m, 0<=out<=1, Int);

        nrows = size(A,1);
        A = reshape(A,nrows,nx)
        B = reshape(B,nrows,ny)
    
        @constraint(m, [i = 1:s], A*x + B*y[1:ny,i] .<=b);
        #println(bigMarray)
        
        for i = 1:s
            Ind = findall(a->a==1, vec(K[i,:]));
            #if z[i] is activated, then sum(A[Ind,:]*x + B[Ind,:]*y[1:ny,i] - b[Ind]) == 0.
            #Otherwise, if z[i] = 0, then the constraint is trivially verified.
            #The constant -1e6 is a Big-M constraint.
            #println(bigMarray[Ind]*(1-z[i]))

            # Constraints are now decoupled. It leads to better performance
            for j in Ind
                @constraint(m, sum( A[j,k]*x[k] for k=1:nx) + sum(B[j,k]*y[k,i] for k=1:ny)- b[j] - bigMarray[j]*(1-z[i])  >= 0);
            end
            #@constraint(m, sum( A[Ind,:]*x + B[Ind,:]*y[1:ny,i] - b[Ind] - bigMarray[Ind]*(1-z[i]) ) .>= 0);
        end
    
        @constraint(m, sum(z) <= 2*s*(1-out));
    
        @objective(m,Max, sum(z));
        return m, z, x, out, y;
end
