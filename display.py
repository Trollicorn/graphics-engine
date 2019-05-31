from subprocess import Popen, PIPE
from matrix import new_matrix
from os import remove

#constants
XRES = 500
YRES = 500
MAX_COLOR = 255
RED = 0
GREEN = 1
BLUE = 2
windows = False

DEFAULT_COLOR = [0,0,0]

def new_screen( width = XRES, height = YRES ):
    return [[DEFAULT_COLOR[:] for x in range(width)] for y in range(height)]

def new_zbuffer( width = XRES, height = YRES ):
    return new_matrix(width,height,float("-inf"))

def plot( screen, zbuffer, color, x, y, z ):
    newy = YRES - 1 - int(y)
    z = int(z*1000)/1000.0
    if ( x >= 0 and x < XRES and newy >= 0 and newy < YRES ):
        if z > zbuffer[newy][int(x)]:
            zbuffer[newy][int(x)] = z
            screen[newy][int(x)] = color[:]

def clear_screen( screen ):
    for y in range( len(screen) ):
        for x in range( len(screen[y]) ):
            screen[y][x] = DEFAULT_COLOR[:]

def clear_zbuffer(zbuffer):
    for y in range( len(zbuffer) ):
        for x in range( len(zbuffer[y]) ):
            zbuffer[y][x] = float("-inf")

def save_ppm( screen, fname ):
    f = open( fname, 'w' )
    ppm = 'P3\n' + str(len(screen[0])) +' '+ str(len(screen)) +' '+ str(MAX_COLOR) +'\n'
    for y in range( len(screen) ):
        row = ''
        for x in range( len(screen[y]) ):
            pixel = screen[y][x]
            row+= str( pixel[ RED ] ) + ' '
            row+= str( pixel[ GREEN ] ) + ' '
            row+= str( pixel[ BLUE ] ) + ' '
        ppm+= row + '\n'
    f.write( ppm )
    f.close()

def save_extension( screen, fname ):
    ppm_name = fname[:fname.find('.')] + '.ppm'
    save_ppm( screen, ppm_name )
    p = Popen( ['convert', ppm_name, fname ], stdin=PIPE, stdout = PIPE )
    p.communicate()
    remove(ppm_name)

def make_anim(bname):
    p = Popen( ['convert', '-delay', '1.7', 'anim/' + bname+'*',bname+'.gif'], stdin=PIPE, stdout = PIPE )
    p.communicate()
    p = Popen( ['animate', '-delay', '1.7', 'anim/' + bname+'*'], stdin=PIPE, stdout = PIPE )
    p.communicate()

def display( screen ):
    ppm_name = 'pic.ppm'
    save_ppm( screen, ppm_name )
    if windows:
        p = Popen( ['imdisplay.exe', ppm_name], stdin=PIPE, stdout = PIPE )
    else:
        p = Popen( ['display', ppm_name], stdin=PIPE, stdout = PIPE )
    p.communicate()
    remove(ppm_name)
