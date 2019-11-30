# By a polymake format, we mean a Python-list of Python-ordered-tuples, where those Python-ordered-tuples are to be interpreted as either polyhedral computation points, or polyhedral computation inequaliities

# 28-vertex prismatoid, or 28-facet spindle
PolymakeFormat = [
(1,18,0,0,0,1),
(1,-18,0,0,0,1),
(1,0,0,30,0,1),
(1,0,0,-30,0,1),
(1,0,0,0,30,1),
(1,0,0,0,-30,1),
(1,0,5,0,25,1),
(1,0,5,0,-25,1),
(1,0,-5,0,25,1),
(1,0,-5,0,-25,1),
(1,0,0,18,18,1),
(1,0,0,18,-18,1),
(1,0,0,-18,18,1),
(1,0,0,-18,-18,1),
(1,0,0,18,0,-1),
(1,0,0,-18,0,-1),
(1,0,30,0,0,-1),
(1,0,-30,0,0,-1),
(1,30,0,0,0,-1),
(1,-30,0,0,0,-1),
(1,25,0,0,5,-1),
(1,25,0,0,-5,-1),
(1,-25,0,0,5,-1),
(1,-25,0,0,-5,-1),
(1,18,18,0,0,-1),
(1,18,-18,0,0,-1),
(1,-18,18,0,0,-1),
(1,-18,-18,0,0,-1)]

def zerolistmaker(n):
    listofzeroes = [0]*n
    return listofzeroes

# dStepTheorem(A, v, a): applies Santos' $d$-step Theorem
# A is in polymake format, a Python-list of Python-ordered-tuples to be interpreted as either polyhedral computation points or polyhedral computation inequaliities.
# v is the vertex number (vertex numbers start at 0) to suspend: this vertex number should be in the "first half"
# a is the vertex number (vertex numbers start at 0) to push: this vertex number should be in the "last half"
# returns in polymake format, a Python-list of Python-ordered-tuples to be interpreted as either polyhedral computation points or polyhedral computation inequaliities.
def dStepTheorem(A, v, a):
    kappa = 43
    
    P = Polyhedron(ieqs=A)
    if P.is_compact() == false:
        print "ERROR: Polyhedron is not bounded."
        return
    
    #First convert each row to a list so they are easier to work with
    for i in range(len(A)):
        A[i] = list(A[i])
        
    #Add 0 to the end of each row, effectively embedding the resulting figure in space one dimension higher
    for i in range(len(A)):
        A[i] += [0]
    
    #Perform the one point suspension on the user-specified vector v
    A[v][len(A[v])-1] = kappa
    u=[]
    for i in range(len(A[v])):
        u.append(A[v][i])
    u[len(u)-1] = -kappa
    A.append(u)
    
    #Push the vertex a in the direction of uw
    A[a][len(A[a])-1] = 1
    
    #Finally, convert each row back to a vector
    for i in range(len(A)):
        A[i] = vector(A[i])
    
    return A        


# C is in polymake format, a Python-list of Python-ordered-tuples, where those Python-ordered-tuples are to be interpreted as either polyhedral computation points, or polyhedral computation inequaliities.
def InterfaceClean(C):
    
    P = Polyhedron(ieqs = C)
    
    if P.is_compact() == false:
        print "ERROR: Polyhedron not bounded."
        return
    
    
    # Translate Polytope
    pushvector = min(P.bounding_box())
    pushentry = floor(min(pushvector))
    pushvert = []
    for i in range (len(C[0])-1):
        pushvert += [-pushentry]
    pushvert
    pushvert = vector(pushvert)
    Q = Polyhedron(vertices=[pushvert])
    R = P+Q
    Rineqs = R.inequalities_list()
    
    A = []
    b = []

    for i in range (len(Rineqs)):
        A.append(Rineqs[i][1:])
        b.append(-Rineqs[i][0])

    d = len(A[1])
    for i in range (len(A)):
        A[i] += zerolistmaker(len(A))
        for j in range (d):
            A[i][i+d] = -1

    for i in range(len(A)):
        A[i] = vector(A[i])
    
    b = vector(b)
    A = matrix(A) 
    return [A,b,R]

# Give a matrix LHS integer matrix A and a vector b (optional) which defines P = {y >=0 : Ay = b}
def preprocess(A, b=[]): 
    K = []

    for i in range (A.ncols()):
        k=[]
        for j in range (A.nrows()):
            if A[j,i]==0:
                k.append(0)
            if A[j,i]!=0:
                q = floor(log(abs(A[j,i]),2))
                k.append(q)
        K.append(max(k))

    Ksum = 0

    for i in range (len(K)):
        Ksum += K[i]

    m = Ksum+A.nrows()
    n = Ksum+A.ncols()

    rowsofC = []
    D = []
    for i in range (A.ncols()):
        L = len(rowsofC)
        for j in range (L, L+K[i]):
            row = zerolistmaker(n)
            row[i+j] = 2
            row[i+j+1] = -1
            D += [0]
            rowsofC.append(row)
            if j == L+K[i]-1:
                break

    D += b
    for i in range (A.nrows()):
        row = []
        for j in range (A.ncols()):
            row2 = floor(A[i,j]).digits(2)
            while len(row2) < K[j]+1:
                row2 += [0]
            row += row2
        rowsofC.append(row)

    C = matrix(rowsofC)
    
    return [C,D]

############################################

def PlaneSumEntryForbidden(A, b):
    m = A.nrows()
    n = A.ncols()

    # Check that b is a vector with m entries
    if len(b) != m:
        print "ERROR: Vector incorrect length. Length of b=", len(b)
        return

    # Check that A only has integer entries
    for i in range(m):
        for j in range(n):
            if A[i,j] not in ZZ:
                print "ERROR: Matrix not integral. The", i,j,"th entry is ", A[i,j]
                return

    # Use the polyhedron computation features to simultaneously determine (1) if P is empty (2) if P is bounded and (3) the value of U
    equationsforpolyhedron = []
    inequalitiesforpolyhedron = []
    for i in range(m): # For each row of A, which is an equation
        # now, prepend the RHS value
        rhs = b[i]
        firstentry = -rhs
        this_equation = [firstentry]
        currentrow = A[i,:]
        for j in range(n):
            this_equation.append(A[i,j])
        equationsforpolyhedron.append(this_equation)
    for j in range(n):
        this_inequality = zerolistmaker(n+1)
        this_inequality[j+1] = 1
        inequalitiesforpolyhedron.append(this_inequality)

    P = Polyhedron(ieqs=inequalitiesforpolyhedron, eqns=equationsforpolyhedron)
    
    if P.dim() == -1:
        print "ERROR: the way the LHS matrix and RHS are defined creates the empty set for P."
        return 

    if P.is_compact() == False:
        print "ERROR: the way the LHS matrix and RHS are defined creates an unbounded polyhedron"
        return
    
    Pboundingbox = P.bounding_box() 
    upperbounds = Pboundingbox[1]
    
    U = max(upperbounds)

    # Construct rj. The list rj will store the values of r_j. (The index for j will start at 0.)
    rj = []
    for j in range(n):
        positivesum = 0 # Reset sum of positive coefficients from handling previous variable
        negativesum = 0 # Reset sum of positive coefficients from handling previous variable
        for i in range(m): # in column j, it's now time to examine each entry in the column (there are m entries)
            entry = A[i,j]
            if entry > 0:
                positivesum += entry
            if entry < 0:
                negativesum += abs(entry)
        new_rj_value = max(positivesum, negativesum)
        rj.append(new_rj_value)

    # r is the sum of the rj's
    r = sum(rj)

    # R = {1,...,r} in the paper, but we'll create a list starting at 0 and ending at r-1. 
    R = range(r)

    h = m + 1

    # H = {1,...,h} in the paper, but we'll create a list starting at 0 and ending at h-1
    H = range(h)

    # Start the collected of enabled entries. 
    E = [] #enabled
    #F = [] #forbidden

    # u has r entries, all of which are U
    u = []
    for i in range(r):
        u.append(U)

    # v has r entries, all of which are U
    v = []
    for j in range(r):
        v.append(U)

    # In this code, the indexing of k is off by one (since we will start with k=0).
    w = []
    for k in range(m):
        # First, find the value of \sum_{j \in J-}|a_{k,j}|, which will need to be multiplied by U
        Umultiplier = 0 # Clear the value from the previous run, to start setting up the sum.
        for j in range(n):
            entry = A[k,j]
            if entry < 0:
                Umultiplier += abs(entry)
        w.append(b[k] + U*Umultiplier)
    # At this point, the first m values of w were set correctly. Now set the m+1 value.
    sum_of_first_m_ws = 0
    for k in range(m):
        sum_of_first_m_ws += w[k]
    w.append(r*U - sum_of_first_m_ws)

    # Determine the set of enabled and/or forbidden entries (E and F)
    for j1 in range(n):
        for j2 in range(n):
            

            # Determine how much we'll need to shift by
            shift1 = 0
            for j in range(j1):
                shift1 += rj[j]
            shift2 = 0
            for j in range(j2):
                shift2 += rj[j]

            if j1 == j2:
                # We are in an R_j \times R_j \times H block. Some (but not all) of the entries will be enabled
                j = j1
                if rj[j] == 1:
                    # This is the rare case when r_j = 1
                    # If j=1, then R_j=\{1\}; the box R_1 \times R_1 \times H is a single line \{1\} \times \{1\} \times H
                    # If j=1, this single line will have exactly two enabled entries (1,1,k^+) and (1,1,k^-)
                    used = 0
                    for k in range(h):
                        if k < h-1:
                            # We are in a "usual" slice, which encodes an equation.
                            entry = A[k,j]
                            if entry == 1:
                                E.append((shift1,shift2,k))
                                used += 1
                            elif entry == -1:
                                E.append((shift1,shift2,k))
                                used += 1
                        else:
                            # We are in the "top" slice, where we deal with slacks
                            # At this point, used is either 1 or 2.
                            if used == 1:
                                E.append((shift1,shift2,k))
                else:

                    # Now, we handle the k+'s, k-'s, and additional forbidden entries

                    # k+'s first
                    used = 0
                    for i in range(m): # We are about to go through all entries in column j.
                        entry = A[i,j]
                        if entry > 0: # positive entry. The *value* of this integer entry will determine the number of used
                            for basic_counter in range(entry):
                                E.append((shift1+used, shift1+used, i))
                                # all other heights (beside i) need to be forbidden
                                used += 1

                    # k-'s next: modular arithmetic implemented using percent sign
                    used = 0
                    for i in range(m): # We are about to go through all entries in column j.
                        entry = A[i,j]
                        if entry < 0: # negative entry. The *value* of this integer entry (in absolute value) will determine the number of used
                            entry_value = abs(entry)
                            for basic_counter in range(entry_value):
                                E.append((shift1+used, shift1+1+(used%rj[j]), i))

    return [u, v, w, E, rj]

############################################

def PSToLSInterface(u,v,w,E,F):
    a=u
    b=v
    c=w
    r=len(u)
    h=len(w)
    e={}
    U= min(max(a),max(b))
    for i in range(r):
        for j in range(r):
            for k in range(h):
                if (i,j,k) in E: 
                    e[(i,j,k)] = U
                else:
                    e[(i,j,k)] = 0
    
    return [a,b,c,e]

#############################################

def LineSumEntryFree(A,B,C,e):

    # A,B,C should be lists. e should be a dictionary.

    l = len(A)
    m = len(B)
    n = len(C)
    r = l*m
    c = l+m+n
    U = min(max(A),max(B))
    u = matrix(QQ,r,c)
    v = matrix(QQ,r,3)
    w = matrix(QQ,c,3)

    for i in range(l):
        for j in range(m):
            for k in range(n):
                u[i*m+j,k]=e[(i,j,k)]

    for h in range(l):
        for i in range(m):
            u[i+h*m,n+h] = U

    for h in range(l):
        for i in range(m):
            for j in range(m):
                if i == j:
                    u[i+h*m,j+l+n] = U

    v[:,0] = U

    for i in range(l):
        for j in range(m):
            sum = 0
            for triple in e:
                if triple[0] == i and triple[1] == j:
                    sum += e[triple]
            v[i*m+j,1] = sum

    v[:,2] = U

    for i in range (n):
        w[i,0] = C[i]
        w[i,2] = 0


    for i in range (n):
        sum = 0
        for triple in e:
            if triple[2] == i:
                sum += e[triple]
        w[i,1] = sum - C[i]

    for i in range (n,n+l):
        w[i,0] = U*m-A[i-n]
        w[i,1] = 0
        w[i,2] = A[i-n]

    for i in range (n+l,c):
        w[i,0] = 0
        w[i,1] = B[i-n-l]
        w[i,2] = l*U-B[i-n-l]

    return[u,v,w]

#################################################################################################
#################################################################################################
#################################################################################################

Aprime, dprime, Rprime = InterfaceClean(PolymakeFormat)
#Aprime, dprime = preprocess(Aprime, dprime)
u,v,w,E,R = PlaneSumEntryForbidden(Aprime, dprime)
#a,b,c,e = PSToLSInterface(u,v,w,E,F)
#U,V,W = LineSumEntryFree(a,b,c,e)

#################################################################################################
#################################################################################################
#################################################################################################

PolymakeFormatDim6 = dStepTheorem(PolymakeFormat, 0, 23)
Abarprime, dbarprime, Rbarprime = InterfaceClean(PolymakeFormatDim6)
#Abarprime, dbarprime = preprocess(Abarprime, dbarprime)
u2,v2,w2,E2,R2 = PlaneSumEntryForbidden(Abarprime, dbarprime)
#a2,b2,c2,e2 = PSToLSInterface(u,v,w,E,F)
#U,V,W = LineSumEntryFree(a,b,c,e)
