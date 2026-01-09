from vpython import *
#Web VPython 2.9

# Tyler Knight
# 12/10/2024
# Final Project
# Simulating two galaxies colliding and measuring change in number of orbiting stars

# Setting the scene
scene.align = 'left'
scene.userspin = False
scene.autoscale = False
scene.fov = .01
scene.background = vec(0,-6,0)
scene.height = 550
scene.width = scene.height
scene.range = 36                #Initial zoom
print_options(width=545, height=150)

#Stars in Orbit Graph
orbitGraph = graph(title="Number of Stars Over Time", xtitle="Time", ytitle="Number of Stars", xmin=0, ymin=0, ymax=3000, align="right", width=900, height=300)
capStar = gcurve(color=color.green, label='Stars Still in Orbit')
sInHalo = gcurve(color=color.blue, label='Stars Within Either Dark Matter Halo')
tStars = gcurve(color=color.red, label='Total Stars')

#Initializing Constants
nStar = 1000
nGC = 2         # of Galactic centers
mGC = 5         #mass
cDM = 0.09      #constant to make dark matter 90% of total mass at r=10
mDM = cDM*mGC   #mass of DM without r^2
mStar = (mGC + 100*mDM)/nStar       #mass of a star for kinetic energy
vel1 = -10      #Initial Vel of GC1
vel2 = 0        #Initial Vel of GC2
x = 20          #for GC1
y = 3           #for GC2
dt = 0.001
halfdt = dt/2
halfdt2 = halfdt*halfdt
pi2 = pi*pi
g = 4*pi2
initialR = .6
gap = 0.25
orbit_num = 2*pi*initialR/gap
theta = gap/initialR
running = False

#Lists
GC = []
dMat = []
stars = []
curve1List = []
curve2List = []

## Loop to Initialize Positions ##
for i in range(nGC):
    GC.append(sphere(pos=vec(x*cos((i)*pi/2),y*sin(-(i)*pi/2),0), radius=(.5), shininess=0, emissive=False, vel=vec(0,0,0), acc=vec(0,0,0)))
    if i <= 0:
        GC[i].color = color.yellow      #GC1
    else:
        GC[i].color = color.white       #GC2
        
    #dark matter halo
    dMat.append(sphere(pos=GC[i].pos, radius=(10), color=color.blue, shininess=0, emissive=True, opacity=0.3, vel=vec(0,0,0), acc=vec(0,0,0)))
    
    #generate nested list
    stars[i] = []
    
    #variables for positioning stars
    initialR = .6
    orbit_num = 2*pi*initialR/gap
    theta = gap/initialR
    
    # star positions and colors
    for j in range(nStar):
        randStar = int(random()*10)
        if randStar == 0 or randStar == 1:
            starColor = vec(1,0,0)
            starSize = 0.12
        elif randStar == 2 or randStar == 3:
            starColor = vec(0.4,0.2,0.6)
            starSize = 0.1
        elif randStar == 4 or randStar == 5 or randStar == 6:
            starColor = vec(0,0,1)
            starSize = 0.08
        elif randStar == 7 or randStar == 8 or randStar == 9:
            starColor = vec(1,1,0)
            starSize = 0.08
        stars[i].append(sphere(pos=vec(0,0,0), radius=starSize, color=starColor, shininess=0, emissive=True, vel=vec(0,0,0), acc=vec(0,0,0)))
        if j > orbit_num:
            initialR += gap       #new radius to generate stars
            orbit_num = 2*pi*initialR/gap    #new 
            theta = gap/initialR    #new change in theta
            orbit_num += j      #makes if statement false
        if j <= orbit_num:
            stars[i][j].pos = GC[i].pos + vec(initialR*cos(j*theta),initialR*sin(j*theta),0)
            
#Trails
curve1 = curve(pos=curve1List, color=color.yellow, radius=0.1)
curve2 = curve(pos=curve2List, color=color.white, radius=0.1)


## Starting Velocities ##
GC[0].vel = vec(vel1,0,0)
GC[1].vel = vec(vel2,0,0)

for i in range(nGC):
    for j in range(nStar):
        r = (stars[i][j].pos-GC[i].pos).mag
        r2 = r*r
        stars[i][j].vel = GC[i].vel + sqrt(g*((mGC/r)+mDM*r))*vec(-(stars[i][j].pos.y-GC[i].pos.y),(stars[i][j].pos.x-GC[i].pos.x),0)/r
        
        
## Acceleration Calculations ##
def calculate_acceleration():
    for i in range(nGC):
        for j in range(i):
            if i != j:
                
                rGCVec = GC[i].pos - GC[j].pos
                rGCMag = (GC[i].pos - GC[j].pos).mag
                r2GC = rGCMag*rGCMag
                
                # Acceleration due to changing mass inside dark matter halo #
                if rGCMag <= 10:
                    GC[i].acc = -g*(rGCVec/rGCMag)*((mGC/r2GC) + mDM)
                    GC[j].acc = g*(rGCVec/rGCMag)*((mGC/r2GC) + mDM)
                
                # Acceleration outside dark matter halo with constant mass #
                elif rGCMag > 10:
                    GC[i].acc = -g*(rGCVec/rGCMag)*((mGC/r2GC) + (mDM*100/r2GC))
                    GC[j].acc = g*(rGCVec/rGCMag)*((mGC/r2GC) + (mDM*100/r2GC))

        #acc of stars
        for j in range(nStar):
            stars[i][j].acc = vec(0,0,0)
            rVec1 = stars[i][j].pos - GC[0].pos
            rVec2 = stars[i][j].pos - GC[1].pos
            radius1 = rVec1.mag
            radius2 = rVec2.mag
            r12 = radius1*radius1
            r22 = radius2*radius2
            
            #within GC1 DM halo
            if radius1 <= 10:
                stars[i][j].acc += -g*(rVec1/radius1)*(mGC/r12 + mDM)
            #outside GC1 DM halo
            elif radius1 > 10:
                stars[i][j].acc += -g*(rVec1/radius1)*(mGC/r12 + mDM*100/r12)
            #within GC2 DM halo
            if radius2 <= 10:
                stars[i][j].acc += -g*(rVec2/radius2)*(mGC/r22 + mDM)
            #outside GC2 DM halo
            elif radius2 > 10:
                stars[i][j].acc += -g*(rVec2/radius2)*(mGC/r22 + mDM*100/r22)


# Verlet algorithm single step
def singleStep():
    for i in range(nGC):
        GC[i].pos += GC[i].vel*dt + GC[i].acc*halfdt2
        GC[i].vel += GC[i].acc*halfdt
        for j in range(nStar):
            stars[i][j].pos += stars[i][j].vel*dt + stars[i][j].acc*halfdt2
            stars[i][j].vel += stars[i][j].acc*halfdt
            
    calculate_acceleration()
    
    #other half of vel calc
    for i in range(nGC):
        GC[i].vel += GC[i].acc*halfdt
        for j in range(nStar):
            stars[i][j].vel += stars[i][j].acc*halfdt
        dMat[i].pos = GC[i].pos
    
    
## Energy Functions ##
#Potential Energy of GCs
def potentialE():
    pE = 0
    r = (GC[0].pos - GC[1].pos).mag #+ vec(0.001,0.001,0)
    r2 = r*r
    r3 = r*r*r
    r4 = r2*r2
    if r <= 10:
        pE = (g*(-(mGC*mGC/r)+2*mGC*mDM*r+(1/3)*mDM*mDM*r3) + g*((mGC*mGC/5)-(40*mGC*mDM)-((2/3)*mDM*mDM*1000)))
    elif r > 10:
        pE = (g*((mGC*mGC/r)-(200*mGC*mDM/r)-((10000/3)*mDM*mDM/r)))
    return pE

#Kinetic Energy of GCs
def kineticE():
    kE = 0
    r = (GC[0].pos - GC[1].pos).mag
    r3 = r*r*r
    v1 = (GC[0].vel).mag
    v2 = (GC[1].vel).mag
    kE = (1/2)*(mGC+(mDM*100))*(v1*v1) + (1/2)*(mGC+(mDM*100))*(v2*v2)
    return kE
    
#Total Energy of Galactic Centers
def totalE():
    tE = kineticE() + potentialE()
    return tE
    

#Energys and number of stars orbiting
def starEnergy():
    global negTE, posTE, inHalo, outHalo
    negTE = 0
    posTE = 0
    inHalo = 0
    outHalo = 0
    for i in range(nGC):
        for j in range(nStar):
            sPE = 0
            r1 = (stars[i][j].pos - GC[0].pos).mag
            r2 = (stars[i][j].pos - GC[1].pos).mag
            if r1 <= r2:
                if r1 <= 10:
                    sPE += g*mStar*(-(mGC/r1)+mDM*r1+(mGC/5)-20*mDM)
                elif r1 > 10:
                    sPE += g*mStar*((mGC/r1)-(100/r1))
            elif r1 > r2:
                if r2 <= 10:
                    sPE += g*mStar*(-(mGC/r2)+mDM*r2+(mGC/5)-20*mDM)
                elif r2 > 10:
                    sPE += g*mStar*((mGC/r2)-(100/r2))
            v = (stars[i][j].vel - GC[i].vel).mag
            v2 = v*v
            sKE = (1/2)*mStar*v2
            sTE = sPE + sKE
            if sTE <= 0:
                negTE += 1
            elif sTE > 0:
                posTE += 1
            if r1 <= 10 or r2 <= 10:
                inHalo += 1
            else:
                outHalo += 1
            
    print(f"Stars in Orbit = {negTE}\nStars in Halo = {inHalo}")
    
    
## Functions for Buttons and Sliders ##
#Start/Stop Button
def startStop():
    global running
    running = not running
    if running:
        ssButton.text = "Stop"
        ssButton.background = color.red
        ssButton.color = color.white
    else:
        ssButton.text = "Start"
        ssButton.background = color.green
        ssButton.color = color.black

#Reset and change to slider values
def reset():
    x = xSlider.value
    y = ySlider.value
    for i in range(nGC):
        initialR = .6
        orbit_num = 2*pi*initialR/gap
        theta = gap/initialR

        GC[i].pos = vec(x*cos((i)*pi/2),y*sin(-(i)*pi/2),0)
        GC[i].clear_trail()
        
        dMat[i].pos = GC[i].pos
        for j in range(nStar):
            if j > orbit_num:
                initialR += gap                     #new radius to generate stars
                orbit_num = 2*pi*initialR/gap       #new 
                theta = gap/initialR                #new change in theta
                orbit_num += j                      #makes if statement false
            if j <= orbit_num:
                stars[i][j].pos = GC[i].pos + vec(initialR*cos(j*theta),initialR*sin(j*theta),0)
                
    vel1 = velocitySlider.value
    vel2 = velocitySlider2.value

    GC[0].vel = vec(vel1,0,0)
    GC[1].vel = vec(vel2,0,0)
            
    for i in range(nGC):
        for j in range(nStar):
            r = (stars[i][j].pos-GC[i].pos).mag
            stars[i][j].vel = GC[i].vel + sqrt(g*((mGC/r)+mDM*r))*vec(-(stars[i][j].pos.y-GC[i].pos.y),(stars[i][j].pos.x-GC[i].pos.x),0)/r

    curve1.clear()
    curve2.clear()


## wtext Functions ##
#Zoom
def zoom():
    zoomSliderReadout.text = zoomSlider.value

#Vel1 
def adjustVelocity():
    velocitySliderReadout.text = velocitySlider.value
    
#Vel2
def adjustVelocity2():
    velocitySlider2Readout.text = velocitySlider2.value
    
# x-pos GC1
def adjustX():
    xSliderReadout.text = -xSlider.value

# y-pos GC2
def adjustY():
    ySliderReadout.text = -ySlider.value

#Speed up/Slow down
def changeRate():
    rateSliderReadout.text = rateSlider.value
            
        
## Button and Slider Objects ##
#start/stop button
ssButton = button(text="Start", bind=startStop, background = color.green)

#reset button
resetButton = button(text="Reset", bind=reset)
scene.append_to_caption("\n\n")

#Zoom Slider
scene.append_to_caption("   Zoom")
scene.append_to_caption("\n\n")
zoomSlider = slider(left=0, min=2, max=300, step=2, value=36, bind=zoom, length=300)
scene.append_to_caption(" Range = ")
zoomSliderReadout = wtext(text="36")

#Potential Energy wtext
pE_text = wtext(text=f"                   Potential Energy = ")
scene.append_to_caption("\n\n")

#Rate Slider
scene.append_to_caption("   Rate")
scene.append_to_caption("\n\n")
rateSlider = slider(left=0, min=2, max=300, step=2, value=60, bind=changeRate, length=300)
scene.append_to_caption(" Rate = ")
rateSliderReadout = wtext(text="60")

#Kinetic Energy wtext
kE_text = wtext(text=f"                     Kinetic Energy = ")
scene.append_to_caption("\n\n")

#Velocity GC1 Slider
scene.append_to_caption("   Velocity")
scene.append_to_caption("\n\n")
velocitySlider = slider(left=0, min=-20, max=20, step=.5, value=-10, bind=adjustVelocity, length=300)
scene.append_to_caption(" V GC1 = ")
velocitySliderReadout = wtext(text="-10")

#Total Energy wtext
tE_text = wtext(text=f"                 Total Energy = ")
scene.append_to_caption("\n\n")

#Velocity GC2 Slider
velocitySlider2 = slider(left=0, min=-20, max=20, step=.5, value=0, bind=adjustVelocity2, length=300)
scene.append_to_caption(" V GC2 = ")
velocitySlider2Readout = wtext(text="0")
scene.append_to_caption("\n\n")

#x-pos for GC1 Slider
scene.append_to_caption("   Position")
scene.append_to_caption("\n\n")
xSlider = slider(left=0, min=0, max=50, step=1, value=20, bind=adjustX, length=300)
scene.append_to_caption(" x Pos GC1 = ")
xSliderReadout = wtext(text="-20")
scene.append_to_caption("\n\n")

#y-pos for GC2 Slider
ySlider = slider(left=0, min=0, max=50, step=1, value=3, bind=adjustY, length=300)
scene.append_to_caption(" y Pos GC2 = ")
ySliderReadout = wtext(text="-3")
scene.append_to_caption("\n\n")
scene.append_to_caption("\n")

#Labels for GCs
GC1Label = label(pos=GC[0].pos, text="GC 1", xoffset=30, yoffset=40, space=5, height=15, border=5, box=False, opacity=0, color=color.yellow)
GC2Label = label(pos=GC[1].pos, text="GC 2", xoffset=-30, yoffset=-40, space=5, height=15, border=5, box=False, opacity=0, color=color.white)

#Variables for main loop
orbitCounter = 30        #Iterator for checking orbits
time = 0                #Time Iterator

while True:
    rate(rateSlider.value)          #Rate from rateSlider
    
    #Update scene center focus to follow GC2
    scene.center = GC[1].pos
    
    #Update Zoom
    scene.range = zoomSlider.value
    
    #Update label pos
    GC1Label.pos = GC[0].pos
    GC2Label.pos = GC[1].pos
        
    if running:
        
        #30 steps per 1 loop
        for i in range(30):
            singleStep()
                        
        #Update Energy wtexts                
        pE_text.text = "                  Potential Energy = {:.3f}".format(potentialE())
        kE_text.text = "                    Kinetic Energy = {:.3f}".format(kineticE())
        tE_text.text = "                Total Energy = {:.3f}".format(totalE())
        
        #Trail update
        curve1.append(pos=GC[0].pos)
        curve2.append(pos=GC[1].pos)
        
        #Checks orbits every 5 loops
        if orbitCounter == 30:
            starEnergy()
            capStar.plot(time, negTE)   #Plot on graph
            sInHalo.plot(time, inHalo)
            tStars.plot(time, 2000)
            orbitCounter = 0            #Reset counter
            
        orbitCounter += 1               #Iterate at end of loop
        time += 1