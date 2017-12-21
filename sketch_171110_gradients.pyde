from Saver import saves
N_PHI = 5000;N_save = 10; frameTrig = 5000

fc = 1; wd = int(1080/fc);ht = int(1080/fc)
nS = fc*.75e-3    #fc/400.   # noise/space step
nxoff = -(wd/2.)*nS + 1400   # noise offsets
nyoff = -(ht/2.)*nS + 2034
dt = .002  # time interval per update (s)
updates = 1  # n of updates per frame
tpf = dt*updates   # time per frame (s)
bgup = 1 #int(5/(tpf))  # update background every s or frames
lspan = tpf    # lifespan s or frames
dgrad = .09/nS  # index separation for calculating gradient
Urange = PVector(.447002,.663417)
hrange = (0,1)
maxF = 3.e0  # max force acting on particles  
maxG = .110038/(2**.5*nS*dgrad) # max gradient value
vini = PVector(1,1).mult(0./fc)   # initial vmags
R = 15./fc     # diameter of particle
maxL = R*1.5; minL = R  #length of vector line
wdoff = htoff = maxL   # space beyond canvas
np = 30   # n particles
mr = (2,2)  # mass range
pa = .15    # particle alpha value
# random range lims
xi = (-wdoff,wd+wdoff); yi = (-htoff,ht+htoff)
# rpos = (0, ((wd/2.)**2 + (ht/2.)**2)**.5)  #radial position
rpos = (wd*.1,wd*.2)
rposoff = PVector(wd/2.,ht/2.)   # radial position center
colorbg = False
Cd = .0    # drag coeff
rgbsft = maxL
#####
# z coord on points: k or constants
# tri = []
# phi = -PI/4
# rr = wd/4.
# for i in xrange(3):
#     th = TWO_PI*i/3
#     x = rr*cos(th+phi)+wd/2. 
#     y = rr*sin(th+phi)+ht/2.
#     tri.append(PVector(x,y,1))
# tri.append(PVector(wd/2.,ht/2.,-1))
# epts = tuple(tri)
sep = .25
# epts = [PVector(wd*.5,ht*.5,mr[0])]
epts = (PVector(0,ht/2.,1),PVector(wd,ht/2.,1))

def settings():
    size(wd,ht)

def setup():
    global bg,P,saving
    colorMode(HSB,1.)
    strokeWeight(1)
    noFill()
    noiseSeed(2456)
    noiseDetail(4)
    saving = saves(N_PHI,N_save)
    # saving.onClick()
    P = psyst(np,mr)
    bg = createImage(wd,ht,RGB)
    bg.loadPixels()
    umax = 0; umin = 1
    hmax = 0; hmin = 1
    gmax = 0; bmax = 0
    if colorbg:
        for ind in xrange(wd*ht):
            x = ind%wd
            y = ind/wd
            u = nice(x,y,0,0)
            h = Umap(u)
            if u > umax: umax = u; hmax = h
            elif u < umin: umin = u; hmin = h
            fx,fy = calcGrad(x,y)
            gmg = (fx**2+fy**2)**.5
            b = gmap(gmg)
            if gmg > gmax: gmax = gmg; bmax = b
            # h -= 0
            if not(0<=h<1): h -= floor(h)
            # if not(0<=b<1): b -= floor(b)
            bg.pixels[ind] = color(h,1,b) #color(h,1,b*.7+.3)
        print 'U', umin, umax
        print 'h', hmin,hmax
        print 'gmax',gmax, gmax*2**.5*nS*dgrad 
        print 'b', bmax
    else:
        for ind in xrange(wd*ht): 
            bg.pixels[ind] = color(1)
    bg.updatePixels()

def draw():
    global bg
    if frameCount%frameTrig==0: saving.onClick()
    blendMode(BLEND)
    background(1)
    image(bg,0,0)
    blendMode(SUBTRACT)
    P.update(updates)
    if frameCount%bgup == 0:
        bg = get(0,0,width,height)
    saving.save_frame()

###
class particle(object):
    def __init__(self,m,posi,veli,life,coff):
        self.posi = posi.copy()
        self.pos = posi.copy()
        self.veli = veli.copy()
        self.vel = veli.copy()
        self.grad = PVector(0,0)
        self.acc = PVector(0,0)
        self.m = m
        self.life = life
        self.lifespan = life
        self.coff = coff
        self.l = m*veli.mag()
    
    def update(self):
        # find gradient
        self.getGrad()
        # drag = self.vel.copy()
        # drag.setMag(self.vel.magSq()*Cd)
        # F = PVector.sub(self.grad,drag)
        self.acc = PVector.div(self.grad,self.m)
        self.vel.add(PVector.mult(self.acc,dt))
        self.pos.add(PVector.mult(self.vel,dt))
        # self.pos = PVector(mouseX,mouseY)
        self.oob(False,False)  # out of bounds calc
        self.acc = PVector(0,0)
        self.life -= dt
    
    def oob(self,die,edge):
        """ edge: bounce on edge"""
        out = False
        bounced = False
        if not(0 <= self.pos.x < width):
            if edge:
                self.vel = PVector(self.vel.x*-1,self.vel.y)
                bounced = True
            elif die:
                out = True
            else:
                self.pos.x -= width*abs(self.pos.x)/self.pos.x
        if not(0 <= self.pos.y < height):
            if edge:
                self.vel = PVector(self.vel.x,self.vel.y*-1)
                bounced = True
            elif die:
                out = True
            else:
                self.pos.y -= height*abs(self.pos.y)/self.pos.y
        if out: self.life = -9999
        elif bounced:
            v = PVector.mult(self.vel,dt)
            self.pos.add(v)
    
    def show(self,):
        b = gmap(self.grad.mag()/maxF)*rgbsft
        # b = int(.5+.5*cos(TWO_PI*frameCount/600))
        # h = Umap(nice(self.pos.x,self.pos.y,0,0))+self.coff
        # h = int(Umap(nice(self.pos.x,self.pos.y,0,0))*3)
        # h = ((h+b)%3)/3.
        h = 1/3.
        # if not(0<=h<1): h -= floor(h)
        L = map(self.grad.mag(),0,maxG,minL,maxL)
        gd = self.grad.copy(); gd.setMag(L)
        p = PVector.add(self.pos,gd)
        # alpha compensator; long lines less alpha
        aoff = map(self.grad.mag(),0,maxG,pa*.75,pa)
        # stroke(color(h,1,1,aoff))
        # point(self.pos.x,self.pos.y)
        # line(self.pos.x,self.pos.y,p.x,p.y)
        ### for tiling lines \/
        # gd = PVector(gd.y,-gd.x) # rotate 90CW
        gd.setMag(b)
        # xo = not(0<=self.pos.x<width)
        # yo = not(0<=self.pos.y<height)
        # pxo = not(0<=p.x<width)
        # pyo = not(0<=p.y<height)
        xo = (self.pos.x-gd.x) < 0 or (p.x+gd.x) < 0 or \
             (self.pos.x-gd.x) >= width or (p.x+gd.x) >= width
        yo = (self.pos.y-gd.y) < 0 or (p.y+gd.y) < 0 or \
             (self.pos.y-gd.y) >= height or (p.y+gd.y) >= height
        nl = False; npx = npy = 0
        if xo:
            npx = -width*sign(self.pos.x-width/2.)
            nl = True
        if yo:
            npy = -height*sign(self.pos.y-height/2.)
            nl = True
        if nl: aoff /= 2.
        for i in xrange(-1,2):
            stroke(color(h+i*1/3.,1,1,aoff))
            if nl:
                line(self.pos.x+npx+i*gd.x,self.pos.y+npy+i*gd.y,
                     p.x+npx+i*gd.x,p.y+npy+i*gd.y)
                # ellipse(self.pos.x+npx+i*gd.x,self.pos.y+npy+i*gd.y,
                #         4,4)
            line(self.pos.x+i*gd.x,self.pos.y+i*gd.y,
                 p.x+i*gd.x,p.y+i*gd.y)
            # ellipse(self.pos.x+i*gd.x,self.pos.y+i*gd.y,4,4)
        ### for tiling lines /\
        # ellipse(self.pos.x,self.pos.y,R,R)
        # vel vector
        # stroke(1,pa/2)
        # v = PVector.mult(self.vel,maxL/20)
        # if v.mag < minL: v.setMag(minL)
        # v.add(self.pos)
        # line(self.pos.x,self.pos.y,v.x,v.y)
    
    def getGrad(self):
        fx,fy = calcGrad(self.pos.x,self.pos.y)
        self.grad = PVector(fx,fy)
        self.grad.mult(maxF)
    
    def reset(self,orig):
        if orig: self.pos = self.posi.copy()
        else:
            # radial
            # theta = random(0,TWO_PI)
            # Rad = random(rpos[0],rpos[1])
            # self.pos = PVector(Rad*cos(theta) + rposoff.x,
            #                    Rad*sin(theta) + rposoff.y)
            # rectangular
            self.pos = PVector(random(xi[0],xi[1]),
                               random(yi[0],yi[1]))
        self.vel = self.veli.copy()
        self.acc = PVector(0,0)
        self.life = self.lifespan
        self.getGrad()

###
class psyst(object):
    def __init__(self,n,m):
        self.p = []
        for i in xrange(n):
            # rectangular
            # pos = PVector(wd*.5,
            #               ht*.75)
            pos = PVector(int(random(xi[0],xi[1])),
                          int(random(yi[0],yi[1])))
            # radial
            # theta = random(0,TWO_PI)
            # Rad = random(rpos[0],rpos[1])
            # pos = PVector(Rad*cos(theta) + rposoff.x,
            #               Rad*sin(theta) + rposoff.y)
            ## v perpendicular (mostly for single source)
            # v = PVector.sub(pos,rposoff)
            # sign = 0
            # while sign == 0:
            #     sign = random(-1,1)
            # sign = abs(sign)/sign
            # v.setMag(random(vini.x,vini.y)*sign)
            # v = PVector(v.y,-v.x)
            ## v2
            v = PVector(1,0).setMag(vini.x)
            ms = random(m[0],m[1])
            c = 0   # hshift from given by U func
            pr = particle(ms,pos,v,lspan,c)
            self.p.append(pr)
    
    def update(self,up):
        for i in xrange(len(self.p)-1,-1,-1):
            prnt = True
            for u in xrange(up):
                # print self.p[i].pos,self.p[i].grad
                self.p[i].update()
                if self.p[i].life <= 0:
                    self.p[i].reset(False)
                    # self.p.pop(i)
                    # print i,' died'
                    # prnt = False
                    break
            if prnt: self.p[i].show()
        
#####
def nice(x,y,dgx,dgy):
    # tiled noise
    xp = 1.*x/wd
    yp = 1.*y/ht
    h00 = noise((x+dgx)*nS+nxoff,(y+dgy)*nS+nyoff)
    h01 = noise((x+dgx)*nS+nxoff,(y+ht+dgy)*nS+nyoff)
    h10 = noise((x+wd+dgx)*nS+nxoff,(y+dgy)*nS+nyoff)
    h11 = noise((x+wd+dgx)*nS+nxoff,(y+ht+dgy)*nS+nyoff)
    h = xp*yp*h00 + xp*(1-yp)*h01 + \
        (1-xp)*yp*h10 + (1-xp)*(1-yp)*h11
    return h
    ### math pot
    # U = 0
    # for i in xrange(len(epts)):
    #     xr = (x-epts[i].x+dgx)
    #     yr = (y-epts[i].y+dgy)
    #     r = (fc*(xr*xr+yr*yr)**.5*nS + nS)
    #     # U += epts[i].z/r
    #     U -= epts[i].z/(r**2)
    # return U

def calcGrad(x,y):
    fx = (nice(x,y,dgrad,0)-nice(x,y,-dgrad,0))/(nS*2*dgrad)
    fy = (nice(x,y,0,dgrad)-nice(x,y,0,-dgrad))/(nS*2*dgrad)
    return -fx,-fy

def Umap(u):
    h = map(u,Urange.x,Urange.y,0,1)
    return map(h,0,1,hrange[0],hrange[1])

def gmap(G):
    return map(G,0,maxG,0,1)

def sign(s):
    if s >= 0: return 1
    else: return -1
    
def mouseClicked():
    saving.onClick()

def keyPressed():
    if key == 'f': print frameCount,frameCount*tpf