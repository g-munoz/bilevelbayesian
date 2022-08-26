using LinearAlgebra, JuMP, Gurobi, Polymake, Printf

function AllVertex(A,B,b)

    (m,nx) = size(A);
    (ny,) = size(B');

    ###############################################
    #We recover all the faces of D = {Ax + By <= b}
    #of dimension nx. These are given by the nonempty
    #intersection of ny affine spaces with linearly
    #independent normal vectors.
    ###############################################

    (K,s) = getFaces_pmk(A,B,b,m,nx,ny) #Warning: this function is implicitly assuming that all inequalities are 
    
    (m,zvar,xvar,outvar) = MIPVertexSearchBase(A,B,b,K,s,nx,ny);
    X = zeros(nx,1)

    while true
        status = optimize!(m)
        z = JuMP.value.(zvar)
        x = JuMP.value.(xvar)
        out = JuMP.value.(outvar);
        value = objective_value(m)
        if out==1
            break
        end

        X = hcat(X,x);

        #We force at least one constraint z[k] to be active in the complement
        #of z. If this is not possibe, out is activated
        Ind = findall(a->a==0, vec(z));
        @constraint(m,sum(zvar[Ind]) + outvar>=1);

        @constraint(m, sum(zvar) <= value); #this may be a bad idea, since it is adding a ineq parallel to obj
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
    
        m = Model(Gurobi.Optimizer)
    
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
    
        return m, z, x, out;
end
