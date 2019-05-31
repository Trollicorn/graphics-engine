from compiler import mdl
from display import *
from matrix import *
from draw import *
from lighting import *
from transform import *

def run(filename):
    """
    This function runs an mdl script
    """
    p = mdl.parseFile(filename)
    if p:
        (commands, symbols) = p
    else:
        print("Parsing failed.")
        return
    view = [0,
            0,
            1];
    ambient = [50,
               50,
               50]
    light = [[0.5,
              0.75,
              1],
             [255,
              255,
              255]]
    color = [50, 50, 50]
    areflect = [0.1,
                0.1,
                0.1]
    dreflect = [0.5,
                0.5,
                0.5]
    sreflect = [0.5,
                0.5,
                0.5]

    step_3d = 100
    consts = ''
    coords = []
    coords1 = []
    symbols['.white'] = ['constants',
                         {'red': [0.2, 0.5, 0.5],
                          'green': [0.2, 0.5, 0.5],
                          'blue': [0.2, 0.5, 0.5]}]
    reflect = '.white'
    transform = {
        "move": translate,
        "rotate": rotate,
        "scale": dilate
    }
    shape = {
        "line": add_edge,
        "circle": circle,
        "hermite": hermite,
        "bezier": bezier
    }
    solid = {
        "box": box,
        "sphere": sphere,
        "torus": torus
    }
    frames = 1
    basename = 'anim/frame'
    bname = 'frame'
    knobs = []
#    print(symbols)
#    print csystems
    for command in commands:
        op = command['op']
        if op == 'frames':
            frames = command['args'][0]
        elif op == 'basename':
            bname = command['args'][0]
            basename = 'anim/'+bname
###########---------------------CHECK VARY COUNT, IMPLEMENT OPITIMIZATION LATER
    if frames != 1:
        knobs = [{} for i in range(int(frames))]
        for command in commands:
            op = command['op']
            if op == 'vary':
                args = command['args']
                rang = args[1]-args[0]
                knob = command['knob']
                start = args[2]
                end = args[3]
                step = (end-start)/rang
                for d in range(int(args[0]),int(args[1])+1):
                        knobs[d][knob] = start
                        start += step
#        print(knobs)
#    print(commands)
    for f in range(int(frames)):
        print('frame',f,'of',int(frames-1))
        tmp = new_matrix()
        ident( tmp )
        csystems = [ [x[:] for x in tmp] ]
        screen = new_screen()
        zbuffer = new_zbuffer()
        tmp = []
        for command in commands:
            op = command['op']
        #    print(op)
            if op == 'constants' or op == 'vary' or op == 'basename' or op == 'frames':
                pass
            elif op == 'push':
                csystems.append(duplicate(csystems[-1]))
            elif op == 'pop':
                del csystems[-1]
            elif op in transform:
                if command['knob'] is None:
                    transform[op](csystems[-1],command['args'],1)
                else:
                    transform[op](csystems[-1],command['args'],knobs[f][command['knob']])
            elif op in solid:
                solid[op](tmp,command['args'])
                matrix_mult(csystems[-1],tmp)
                if command['constants'] is None:
                    draw_polygons(tmp,screen,zbuffer,color,view, ambient, light, areflect, dreflect, sreflect)
                else:
                    const = symbols[command['constants']][1]
                    reflects = [[const['red'][i],const['green'][i],const['blue'][i]] for i in range(3)]
                    draw_polygons(tmp,screen,zbuffer,color,view, ambient, light, reflects[0], reflects[1], reflects[2])
                tmp = []
            elif op in shape:
                shape[op](tmp,command['args'])
                matrix_mult(csystems[-1],tmp)
                draw_lines(tmp,screen,zbuffer,color)
                tmp = []
            elif op == 'display':
                display(screen)
            elif op == 'save' and frames == 1:
                save_extension(screen,command['args'][0]+'.png')
        #    print("done")
        if frames != 1:
            basenum = str(f)
            if f < 100:
                basenum = '0' + basenum
            if f < 10:
                basenum = '0' + basenum
            saved = basename+basenum+'.png'
            save_extension(screen,saved)
            print("saved as",saved)
    if frames != 1:
        make_anim(bname)
        print("saved as",bname+'.gif')
            #else:
        #    print(command['op'])
        #    print(command)
        #print(knobs)
