import numpy as np

def Lanczos_repo(mat, *args):
    # Knobs
    symTol = 1e-8
    # Check input matrix size
    m, n = mat.shape
    if (m != n):
        raise ValueError('Input matrix must be square')
    # Make sure input matrix is symmetric
    if np.max(np.max(np.abs(mat - mat.T))) > symTol:
        raise ValueError('Input matrix is not symmetric to working precision')
    # Parse user inputs
    if (len(args) == 1):
        if (len(args[0]) == n):
            k = n
            v = args[0]
        else:
            k = args[0]
            v = np.random.randn(n)
    elif (len(args) == 2):
            k = args[0]
            v = args[1]
    else:
        k = n
        v = np.random.randn(n)
    # Initialize variables
    Q = np.zeros((n,k))
    q = v / np.linalg.norm(v)
    Q[:,0] = q
    d = np.zeros(k)
    od = np.zeros(k-1)
    # Perform Lanczos iterations
    for i in range(k):
        z = mat @ q
        d[i] = q.T @ z

        #----------------------------------------------
        # Uncomment only one of the following 3 options
        #----------------------------------------------
        # Option 1: Full re-orthogonalization
        #z = z - Q[:,0:i+1] @ (Q[:,0:i+1].T @ z)
        
        # Option 2: Full re-orthogonalization (x2)
        z = z - Q[:,0:i+1] @ (Q[:,0:i+1].T @ z)
        z = z - Q[:,0:i+1] @ (Q[:,0:i+1].T @ z)
        
        # Option 3: No re-orthogonalization
        #z = z - d[i] * q
        #if (i > 0):
        #    z = z - od[i-1] * Q[:,i-1]
        #----------------------------------------------

        if (i != k-1):
            od[i] = np.linalg.norm(z)
            q = z / od[i]
            Q[:,i + 1] = q
    # Construct T
    T = np.diag(d) + np.diag(od,-1) + np.diag(od,1)
    # Return user-requested information
    if (len(args) == 0 or len(args) == 1):
        return T
    elif (len(args) == 2):
        return T, Q
