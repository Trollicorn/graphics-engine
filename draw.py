from display import *
from matrix import *
from lighting import *
from math import cos,sin,pi

def draw_lines( matrix, screen, zbuffer,color ):
    for i in range(0,len(matrix)-1,2):
        draw_line(matrix[i][0],matrix[i][1],matrix[i][2],matrix[i+1][0],matrix[i+1][1],matrix[i+1][2],screen,zbuffer,color)

def add_edge( matrix, args): #[x0, y0, z0, x1, y1, z1]
    add_point(matrix,args[0],args[1],args[2])
    add_point(matrix,args[3],args[4],args[5])

def add_point( matrix, x, y, z=0 ):
    matrix.append([x,y,z,1])

def add_poly(polygon,x0,y0,z0,x1,y1,z1,x2,y2,z2):
    add_point(polygon,x0,y0,z0)
    add_point(polygon,x1,y1,z1)
    add_point(polygon,x2,y2,z2)

def draw_polygons(polygons, screen, zbuffer,color, view, ambient, light, areflect, dreflect, sreflect):
    for i in range(0,len(polygons)-1,3):
        norm = surf(polygons,i)
        if norm[2] > 0:
            color = get_lighting(norm, view, ambient, light, areflect, dreflect, sreflect)
            scanline(polygons[i],polygons[i+1],polygons[i+2],screen,zbuffer,color)

def scanline(c0,c1,c2,screen,zbuffer,color):
    corners = [c0,c1,c2]
    corners.sort(key=lambda x:x[1])
    bot = corners[0]
    mid = corners[1]
    top = corners[2]
    Bx = x0 = x1 = bot[0]
    Bz = z0 = z1 = bot[2]
    By = int(bot[1])
    Tx = top[0]
    Tz = top[2]
    Ty = int(top[1])
    Mx = mid[0]
    Mz = mid[2]
    My = int(mid[1])
    BtoT = Ty - By * 1.0 + 1
    BtoM = My - By * 1.0 + 1
    MtoT = Ty - My * 1.0 + 1
    switch = False
    dx0 = (Tx-Bx)/BtoT if BtoT != 0 else 0
    dz0 = (Tz-Bz)/BtoT if BtoT != 0 else 0
    dx1 = (Mx-Bx)/BtoM if BtoM != 0 else 0
    dz1 = (Mz-Bz)/BtoM if BtoM != 0 else 0
    for y in range (By,Ty+1):
        if not switch and y >= My:
            dx1 = (Tx-Mx)/MtoT if MtoT != 0 else 0
            dz1 = (Tz-Mz)/MtoT if MtoT != 0 else 0
            x1 = Mx
            z1 = Mz
            switch = True
        draw_line(int(x0),y,z0,int(x1),y,z1,screen,zbuffer,color)
        x0 += dx0
        z0 += dz0
        x1 += dx1
        z1 += dz1

def draw_line( x0, y0, z0, x1, y1, z1, screen, zbuffer, color ):
    #swap points if going right -> left
    if x0 > x1:
        xt = x0
        yt = y0
        zt = z0
        x0 = x1
        y0 = y1
        z0 = z1
        x1 = xt
        y1 = yt
        z1 = zt
    x = x0
    y = y0
    z = z0
    A = 2 * (y1 - y0)
    B = -2 * (x1 - x0)
    wide = False
    tall = False
    if ( abs(x1-x0) >= abs(y1 - y0) ): #octants 1/8
        wide = True
        loop_start = x
        loop_end = x1
        dx_east = dx_northeast = 1
        dy_east = 0
        d_east = A
        distance = x1 - x + 1
        if ( A > 0 ): #octant 1
            d = A + B/2
            dy_northeast = 1
            d_northeast = A + B
        else: #octant 8
            d = A - B/2
            dy_northeast = -1
            d_northeast = A - B
    else: #octants 2/7
        tall = True
        dx_east = 0
        dx_northeast = 1
        distance = abs(y1 - y) + 1
        if ( A > 0 ): #octant 2
            d = A/2 + B
            dy_east = dy_northeast = 1
            d_northeast = A + B
            d_east = B
            loop_start = y
            loop_end = y1
        else: #octant 7
            d = A/2 - B
            dy_east = dy_northeast = -1
            d_northeast = A - B
            d_east = -1 * B
            loop_start = y1
            loop_end = y
    dz = (z1 - z0) / distance if distance != 0 else 0
    while ( loop_start < loop_end ):
        plot( screen, zbuffer, color, x, y, z )
        if ( (wide and ((A > 0 and d > 0) or (A < 0 and d < 0))) or
             (tall and ((A > 0 and d < 0) or (A < 0 and d > 0 )))):
            x+= dx_northeast
            y+= dy_northeast
            d+= d_northeast
        else:
            x+= dx_east
            y+= dy_east
            d+= d_east
        z+= dz
        loop_start+= 1
    plot( screen, zbuffer, color, x, y, z )

#------------SOLIDS--------------
def box(polygon, args): #[x,y,z,w,h,d]
    x = args[0]
    y = args[1]
    z = args[2]
    w = args[3] #width
    h = args[4] #height
    d = args[5] #depth
    ex = x + w
    ey = y - h
    ez = z - d
    #parallel xy plane
    add_poly(polygon,ex, y, z, x, y, z, x,ey, z)
    add_poly(polygon, x,ey, z,ex,ey, z,ex, y, z)
    add_poly(polygon, x, y,ez,ex, y,ez,ex,ey,ez)
    add_poly(polygon,ex,ey,ez, x,ey,ez, x, y,ez)
    #parallel yz plane
    add_poly(polygon, x, y,ez, x,ey,ez, x,ey, z)
    add_poly(polygon, x,ey, z, x, y, z, x, y,ez)
    add_poly(polygon,ex, y, z,ex,ey, z,ex,ey,ez)
    add_poly(polygon,ex,ey,ez,ex, y,ez,ex, y, z)
    #parallel xz plane
    add_poly(polygon, x, y,ez, x, y, z,ex, y, z)
    add_poly(polygon,ex, y, z,ex, y,ez, x, y,ez)
    add_poly(polygon, x,ey, z, x,ey,ez,ex,ey,ez)
    add_poly(polygon,ex,ey,ez,ex,ey, z, x,ey, z)

def sphere(polygon,args): #[x,y,z,r]
    x = args[0]
    y = args[1]
    z = args[2]
    r = args[3]
    n = 20
    p = p_sphere(x,y,z,r,n)
    for i in range(len(p)-1):
        a = i + 1
        b = (i+n)%(n*n)
        c = (i+n+1)%(n*n)
        if i % n == 0:
            add_poly(polygon,p[i][0],p[i][1],p[i][2],
                             p[a][0],p[a][1],p[a][2],
                             p[c][0],p[c][1],p[c][2])
        elif i % n == n-2:
            add_poly(polygon,p[i][0],p[i][1],p[i][2],
                             p[a][0],p[a][1],p[a][2],
                             p[b][0],p[b][1],p[b][2])
        else:
            add_poly(polygon,p[i][0],p[i][1],p[i][2],
                             p[a][0],p[a][1],p[a][2],
                             p[b][0],p[b][1],p[b][2])
            add_poly(polygon,p[a][0],p[a][1],p[a][2],
                             p[c][0],p[c][1],p[c][2],
                             p[b][0],p[b][1],p[b][2])

def torus(polygon,args): #[x,y,z,r1,r2]
    x = args[0]
    y = args[1]
    z = args[2]
    r1 = args[3] #small circles
    r2 = args[4] #big circle
    n = 20
    p = p_torus(x,y,z,r1,r2,n)
    for i in range(len(p)):
        a = (i+1)%n+i//n*n
        b = (a + n) % (n*n) #((i+1)%n+i//n*n+n) %(n*n)
        c = (i+n)%(n*n)
        add_poly(polygon,p[i][0],p[i][1],p[i][2],
                         p[a][0],p[a][1],p[a][2],
                         p[c][0],p[c][1],p[c][2])
        add_poly(polygon,p[a][0],p[a][1],p[a][2],
                         p[b][0],p[b][1],p[b][2],
                         p[c][0],p[c][1],p[c][2])

def p_sphere(x,y,z,r,n):
    points = []
    num = float(n)
    for i in range(int(num)):
        phi = 2*pi*i/num
        cosphi = cos(phi)
        sinphi = sin(phi)
        for j in range(int(num)):
            theta = pi*j/(num-1)
            sintheta = sin(theta)
            costheta = cos(theta)
            points.append([int(r*costheta)+x, int(r*sintheta*cosphi)+y, int(r*sintheta*sinphi)+z])
    return points

def p_torus(x,y,z,r1,r2,n):
    points = []
    num = float(n)
    for i in range(int(num)):
        phi = 2*pi*i/num
        cosphi = cos(phi)
        sinphi = sin(phi)
        for j in range(int(num)):
            theta = 2*pi*j/num
            sintheta = sin(theta)
            costheta = cos(theta)
            points.append([int(cosphi*(r1*costheta+r2))+x, int(r1*sintheta)+y, int(sinphi*(r1*costheta+r2))+z])
    return points

#-------------CURVES---------------
def circle(edge, args): #[cx,cy,cz,r]
    cx = args[0]
    cy = args[1]
    cz = args[2]
    r = args[3]
    num = 100.0
    for t in range(int(num)):
        tr = t/num
        x0 = int(r*cos(2*pi*tr)+cx)
        y0 = int(r*sin(2*pi*tr)+cy)
        x1 = int(r*cos(2*pi*(tr+1.0/num))+cx)
        y1 = int(r*sin(2*pi*(tr+1.0/num))+cy)
        stuff = [x0,y0,cz,x1,y1,cz]
        add_edge(edge,stuff)

def bezier(edge, args): #[x0,y0,x1,y1,x2,y2,x3,y3]
    x0 = args[0]
    y0 = args[1]
    x1 = args[2]
    y1 = args[3]
    x2 = args[4]
    y2 = args[5]
    x3 = args[6]
    y3 = args[7]
    b = [[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]]
    #|-1  3 -3  1|
    #| 3 -6  3  0|
    #|-3  3  0  0|
    #| 1  0  0  0|
    g = [[x0,x1,x2,x3],[y0,y1,y2,y3]]
    matrix_mult(b,g)
    # now g = [[ax,bx,cx,dx],[ay,by,cy,dy]]
    funx = lambda t: int(g[0][0]*t*t*t + g[0][1]*t*t + g[0][2]*t + g[0][3])
    funy = lambda t: int(g[1][0]*t*t*t + g[1][1]*t*t + g[1][2]*t + g[1][3])
    num = 100.0
    for t in range(int(num)):
        tr = t/num
        lx0 = funx(tr)
        ly0 = funy(tr)
        lx1 = funx(tr+1.0/num)
        ly1 = funy(tr+1.0/num)
        stuff = [lx0,ly0,0,lx1,ly1,0]
        add_edge(edge,stuff)

def hermite(edge, args): #[x0,y0,x1,y1,rx0,ry0,rx1,ry1]
    x0 = args[0]
    y0 = args[1]
    x1 = args[2]
    y1 = args[3]
    rx0 = args[4]
    ry0 = args[5]
    rx1 = args[6]
    ry1 = args[7]
    h = [[2,-3,0,1],[-2,3,0,0],[1,-2,1,0],[1,-1,0,0]]
    #| 2 -2  1  1|
    #|-3  3 -2 -1|
    #| 0  0  1  0|
    #| 1  0  0  0|
    g = [[x0,x1,rx0,rx1],[y0,y1,ry0,ry1]]
    matrix_mult(h,g)
    # now g = [[ax,bx,cx,dx],[ay,by,cy,dy]]
    funx = lambda t: int(g[0][0]*t*t*t + g[0][1]*t*t + g[0][2]*t + g[0][3])
    funy = lambda t: int(g[1][0]*t*t*t + g[1][1]*t*t + g[1][2]*t + g[1][3])
    num = 100.0
    for t in range(int(num)):
        tr = t/num
        lx0 = funx(tr)
        ly0 = funy(tr)
        lx1 = funx(tr+1.0/num)
        ly1 = funy(tr+1.0/num)
        stuff = [lx0,ly0,0,lx1,ly1,0]
        add_edge(edge,stuff)
