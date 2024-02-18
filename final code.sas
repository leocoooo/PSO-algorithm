/********** PSO Algorithm **********/

proc iml;

/*** defining some functions used in this paper ***/

start Rastrigin(X);
	y=2*10 + (x[,1]**2-10*cos(2*constant("pi")*x[,1])) +
	(x[,2]**2-10*cos(2*constant("pi")*x[,2]));
	return y;
finish Rastrigin; /* solution ={0,0} */

start sphere(X); 
    y=(X[, 1] - 1) ** 2 + (X[, 2] - 1) ** 2+ (X[, 3] - 1) ** 2 - 4;
	return y;
finish sphere;  /* solution ={1,1,1} */

/*** Initialization phase ***/

dim=2; /* Indicating the dimension of the function */
S=100; /* The size S of the swarm */
ub=5; /* Defining the upper band */
min_position = j(1, dim, -ub);
max_position = j(1, dim, ub); /* Defining the limits of the initial position */
X = j(S, dim, .); /* Initializing an array to store the generated values */

/* Loop to generate particles automatically */
do j=1 to S; 
	do i = 1 to dim;
    	X[j, i] = min_position[i] + randfun(1, "Uniform") * 2 * max_position[i];
	end;
end; 

Pb = X;/* Assigning each particle its personal best position, which is its current position in the first step  */
V = j(dim, S, 0)`;/* Initializing the particle velocity vector to zero */

/* defining the global best Gb for the first stage */
Gb_image = rastrigin(X[1,]);
do i = 2 to S;
	if rastrigin(X[i,]) < Gb_image then do;
	Gb_image = rastrigin(X[i,]);
	Gb = X[i,];
	end;
end; 

/* Functions */

/* Function that updates the Pb : Personnal Best (from t to t+1)*/ 
start update_Pb(X,Pb,S);
	do i = 1 to S;
   		if rastrigin(X[i,]) < rastrigin(Pb[i,]) then do;
      		Pb[i,] = X[i,];
   		end;
	end;
	return(Pb);
finish update_Pb;

/* Function that updates the Gb : Global Best (from t to t+1)*/ 
start update_Gb(Gb,Pb,S);
	do i = 1 to S;
   		if rastrigin(Pb[i,]) < rastrigin(Gb) then do;
      		Gb = Pb[i,];
   		end;
	end;
	return(Gb);
finish update_Gb;

/* Function that updates the velocity of the particles (from t to t+1) */
start update_V(w,c1,c2,X,V,Pb,Gb);
	r1=RAND('UNIForm');
	r2=RAND('UNIForm');
	V_updated = w*V+ c1*r1*(Pb-X) + c2*r2*(Gb-X); /* update of the velocity  */
	return(V_updated);
finish update_V;


/* Function that updates the position of X (the particles)(from t to t+1) */
start update_X(X,V_updated);
   X_updated = X + V_updated; /* update of the position */
   return X_updated;
finish update_X;

/* Penalty function that resets the coordinates of particles that exceed the search space */
/* This function is much more complicated in the case of problems with constraints */
start Penalty(X,dim,S,ub);
    do j = 1 to dim;
        do i = 1 to S;
            if abs(X[i,j]) > ub then X[i,j] = (-ub) + randfun(1, "Uniform") * 2 * ub;
        end;
    end;
    return(X);
finish;

/*** Final loop with a stopping criterion based on a maximum number of iterations***/

nbr_iter_max =500;/* max iteration criterion*/
nbr_iter=1;

do while(nbr_iter < nbr_iter_max);
	/* Initialising dynamic parameters */
	w = 0.9 - nbr_iter*(0.9-0.4)/nbr_iter_max;
	c1 = 2 - nbr_iter*(2-0.1)/nbr_iter_max;
	c2 = 0.1 + nbr_iter*(2-0.1)/nbr_iter_max;
	/* Reinitialization of V, X, Pb, Gb */
	V = update_V(w,c1,c2,X,V,Pb,Gb); /* new V*/
	X = update_X(X,V); /* new X */
	X = Penalty(X,dim,S,ub); /* X penalized */
	Pb = update_Pb(X,Pb,S); /* new Pb */
	Gb = update_Gb(Gb,Pb,S); /* new Gb */

	nbr_iter = nbr_iter + 1; /* count of iteration */
end;

y_Gb = rastrigin(Gb); /* Computing the image of the global best*/
call NLPFDD(f,gradient,h,"rastrigin",Gb); /* Computing the gradient of the global best */
print Gb y_Gb gradient ; 

quit;
