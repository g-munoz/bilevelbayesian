using LinearAlgebra, JuMP, Gurobi, Polymake, Printf


function verifysol(ver_model, zverify, z_val, outverify, out_val, xverify, yverify)

    s = length(zverify)
    for i = 1 : s
        fix(zverify[i], z_val[i]; force = true)
    end
    fix(outverify, out_val; force = true)

    #println("\nGoing to optimize\n")
    optimize!(ver_model)
    #println("\nDone\n")
    status = termination_status(ver_model)

    if status == INFEASIBLE
        return false, nothing, nothing
    end

    xval = JuMP.value.(xverify)
    yval = JuMP.value.(yverify)

    #println("Opt status: ",status," Val: ", objective_value(ver_model))
    return true, xval, yval
end

function AllVertex(A,B,b,cbFlag)

    (m,nx) = size(A);
    (ny,) = size(B');

    ###############################################
    #We recover all the faces of D = {Ax + By <= b}
    #of dimension nx. These are given by the nonempty
    #intersection of ny affine spaces with linearly
    #independent normal vectors.
    ###############################################

    (K,s) = getFaces_pmk(A,B,b,m,nx,ny) #Warning: this function is implicitly assuming that all inequalities are 
    
    (m,zvar,xvar,outvar,yvar) = MIPVertexSearchBase(A,B,b,K,s,nx,ny);

    ### this copy is for a verifier model
    ### for technical reasons couldn't use copy method
    (ver_model,zvar_ver,xvar_ver,outvar_ver,yvar_ver) = MIPVertexSearchBase(A,B,b,K,s,nx,ny);
    set_silent(ver_model)

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
            if z_val[i] >= 1
                continue
            end
            z_new = copy(round.(z_val))
            z_new[i] = 1
            #@show zvar
            #@show z_val
            #@show z_new
            status = MOI.submit(m, MOI.HeuristicSolution(cb_data), zvar, z_new)
            #println("I submitted a heuristic solution, and the status was: ", status)
            #sleep(5)
            if status == MOI.HEURISTIC_SOLUTION_ACCEPTED
                z_val = z_new
                improved = true
                println("Went in!")
                return
                sleep(5)
            else
                feasible, x_val, y_val = verifysol(ver_model, zvar_ver, z_new, outvar_ver, out_val, xvar_ver, yvar_ver)
                if feasible
                    println("Error: Gurobi rejected a valid better solution, not stopping")
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
                    println("I submitted a heuristic solution, and the status was: ", status,"  out ", out_val)

                    jumpvars = vcat(zvar,xvar,vec(yvar),outvar)
                    jumpvals = vcat(z_new,x_val,vec(y_val),out_val)
                    dictionary = Dict(jumpvars .=> jumpvals)
                    println(primal_feasibility_report(m, dictionary))
                    #sleep(10)
                #else
                #    println("Well rejected solution")
                    #sleep(5)
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
        #if rand() < 0.1
            # You can terminate the callback as follows:
        #    GRBterminate(backend(model))
        #end
        return
    end


    if cbFlag
        #MOI.set(m, MOI.HeuristicCallback(), greedy_callback)
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

function MIPVertexSearchBase(A,B,b,K,s,nx,ny)

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
    
        @constraint(m, [i = 1:s], A*x + B*y[1:ny,i] .<=b);
    
        for i = 1:s
            Ind = findall(a->a==1, vec(K[i,:]));
            #if z[i] is activated, then sum(A[Ind,:]*x + B[Ind,:]*y[1:ny,i] - b[Ind]) == 0.
            #Otherwise, if z[i] = 0, then the constraint is trivially verified.
            #The constant -1e6 is a Big-M constraint.
            @constraint(m, sum( A[Ind,:]*x + B[Ind,:]*y[1:ny,i] - b[Ind] )>= -1e6*(1-z[i])  );
        end
    
        @constraint(m, sum(z) <= 2*s*(1-out));
    
        @objective(m,Max, sum(z));
    
        return m, z, x, out, y;
end
