"""
A matrix will be an N sized list of 4 element lists.
Each individual list will represent an [x, y, z, 1] point.
For multiplication purposes, consider the lists like so:
x0  x1      xn
y0  y1      yn
z0  z1  ... zn
1  1        1
"""
import math

#print the matrix such that it looks like
#the template in the top comment
def print_matrix( matrix ):
    bucket = []
    for r in matrix[0]:
        bucket.append("")
    for c in matrix:
        for r in range(len(c)):
            #print(str(c[r]))
            bucket[r] += str(c[r]) + " "
    for i in range(len(bucket)):
        print(bucket[i])

#turn the paramter matrix into an identity matrix
#you may assume matrix is square
def ident( matrix ):
    n = len(matrix)
    for c in range(n):
        for r in range(n):
            if  c == r:
                matrix[c][r] = 1
            else:
                matrix[c][r] = 0

#multiply m1 by m2, modifying m2 to be the product
#m1 * m2 -> m2
#ASSUMING SQUARE MATRIX
def matrix_mult( m1, m2 ):
    cols = len(m2) #num cols of m2
    rows = len(m1) #num rows of m2 same as num cols of m1
    #list to store values of current column as temp holder
    curr_col = []
    for r in range(rows):
        curr_col.append(0)
    for c in range(cols): #for each col
        for r in range(rows):
            #store values of col in temp holder
            curr_col[r] = m2[c][r]
        for r in range(rows):
            m2[c][r] = 0 #set new value to zero
            for i in range(rows):
                #add product of corresponding values
                m2[c][r]+= curr_col[i] * m1[i][r]

def new_matrix(rows = 4, cols = 4, n = 0):
    m = []
    for c in range( cols ):
        m.append( [] )
        for r in range( rows ):
            m[c].append( n )
    return m

def duplicate(m): #duplicate matrix
    n = [] #to be returned
    for v in m:
        z = [] #one row
        for i in v:
            z.append(i)
        n.append(z)
    return n

def replace(m1,m2): #change contents of m1 to be m2
    for i in range(len(m1)):
        for j in range(len(m1)):
            m1[i][j] = m2[i][j]

def mag(v): #magnitude of vector
    s = 0
    for i in v:
        s += i * i
    return math.sqrt(s)

def norm(v): #normal vector
    s = mag(v)
    for i in range(len(v)):
        if v[i]!=0:
            v[i] = v[i] / s
    return v

def vscale(vector,k): #multiply vector by constant k, returns new vector
    return [ vector[i]*k for i in range(len(vector)) ] #list comprehension

def vadd(v1,v2): #adds elements from v2 to v1 so (v1+v2), returns new vector
    return [ v1[i]+v2[i] for i in range(len(v1)) ] #list comprehension

def vsubtract(v1,v2): #subtract elements of v2 from v1 so (v1-v2), returns new vector
    return [ v1[i]-v2[i] for i in range(len(v1)) ] #list comprehension

def dot(v1,v2):
    s = 0
    for i in range(len(v1)):
        s += v1[i]*v2[i]
    return s

def cross(v1,v2): #cross product of 3d vectors
    x = v1[1]*v2[2] - v1[2]*v2[1]
    y = v1[2]*v2[0] - v1[0]*v2[2]
    z = v1[0]*v2[1] - v1[1]*v2[0]
    return [x,y,z]

def surf(polygon,index):
    a = polygon[index]
    b = polygon[index+1]
    c = polygon[index+2]
    d = [0,0,0]
    e = [0,0,0]
    for i in range(3):
        d[i] = b[i]-a[i]
        e[i] = c[i]-a[i]
    return norm(cross(d,e))
