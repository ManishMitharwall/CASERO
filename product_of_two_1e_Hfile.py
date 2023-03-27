import numpy as np



def product_Of_two_1e_Hfile(O1,O2):
    n1 = O1.shape[0]; n2 = O2.shape[0]; O3 = []
    for l1 in range(n1):
        w1 = O1[l1,0]; c1= O1[l1,1]; d1= O1[l1,4]
        for l2 in range(n2):
            w2 = O2[l2,0]; c2= O2[l2,1]; d2= O2[l2,4]
            if d1 != c2 and d1 != d2:
                O3_new =  [w1*w2 , c1 , c2 , d2 , d1]
                O3.append(O3_new)
            if d1 == c2:
                O3_new =  [w1*w2 ,c1, 0 , 0, d2]
                O3.append(O3_new)
                if d1 != d2:
                    O3_new = [-1*(w1*w2), c1, c2, d1 , d2]
                    O3.append(O3_new)
    return O3



