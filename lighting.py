from matrix import *
  # IMPORANT NOTE

  # Ambient light is represented by a color value

  # Point light sources are 2D arrays of doubles.
  #      - The fist index (LOCATION) represents the vector to the light.
  #      - The second index (COLOR) represents the color.

  # Reflection constants (ka, kd, ks) are represened as arrays of
  # doubles (red, green, blue)

AMBIENT = 0
DIFFUSE = 1
SPECULAR = 2
LOCATION = 0
COLOR = 1
SPECULAR_EXP = 4

#lighting functions
def get_lighting(normal, view, ambient, light, areflect, dreflect, sreflect ):
    A = calculate_ambient(ambient,areflect)
    D = calculate_diffuse(light,dreflect,normal)
    S = calculate_specular(light,sreflect,view,normal)
    color = [A[i]+D[i]+S[i] for i in range(3)]
    return limit_color(color)

def calculate_ambient(alight, areflect): # A = I*K
    return [alight[i]*areflect[i] for i in range(3)]

def calculate_diffuse(light, dreflect, normal): # D = I * K * (N dot L)
    langle = dot(norm(normal),norm(light[LOCATION]))
    colors = light[COLOR]
    return [colors[i]*dreflect[i]*langle for i in range(3)]

def calculate_specular(light, sreflect, view, normal): # S = I * K * ( (2*N*(N dot L) - L) dot V )^e
    normie = norm(normal)
    lighty = norm(light[LOCATION])
    langle = dot(normie,lighty)
    colors = light[COLOR]
    viewy = norm(view)
    return [(colors[i]*sreflect[i]*dot(vsubtract(vscale(normie,dot(normie,lighty)*2),lighty),viewy)**SPECULAR_EXP) for i in range(3)] #list comprenshion

def limit_color(color):
    for c in range(3):
        color[c] = 255 if color[c] > 255 else 0 if color[c] < 0 else int(color[c])
    return color
