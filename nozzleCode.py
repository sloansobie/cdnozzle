from pandas import *
from math import *
from scipy.optimize import brentq
import numpy as np
from stl import mesh


# EITHER PROVIDE MASSFLOW RATE OR THRUST
def generateNozzle(
    CHAMBER_PRESSURE,
    CHAMBER_TEMPERATURE,
    DESIRED_THRUST,
    DESIRED_MASS_FLOW_RATE,
    ALTITUDE,
    GAMMA,
    GAS_CONSTANT,
    THROAT_RADIUS,
    CHAMBER_RADIUS):
    #THE BELOW LINES ARE USED TO FIND OUT THE ISENTROPIC RELATIONS
    if (11000  >  ALTITUDE)  and  (ALTITUDE  <  25000):
        T_0   =  -56.46
        P_0   =  1000  *  (22.65 * exp(1.73 - 0.000157 * ALTITUDE))
    elif (ALTITUDE  >=  25000):
        T_0   =  -131.21 + 0.00299 * ALTITUDE
        P_0   =  1000 * (2.488 * ((T_0 + 273.1) / 216.6) ** -11.388)
    else:
        T_0   =  15.04 - 0.00649 * ALTITUDE
        P_0   =  1000 * (101.29 * ((T_0 + 273.1) / 288.08) ** 5.256)


    _PR_    =    P_0 / CHAMBER_PRESSURE
    PR_2    =   (P_0 / CHAMBER_PRESSURE) ** ((GAMMA - 1) / GAMMA)
    T_T     =   (2 * GAMMA * GAS_CONSTANT * CHAMBER_TEMPERATURE) / (GAMMA - 1)
    P_T     =   ((2 / (GAMMA + 1)) ** (GAMMA / (GAMMA - 1))) * 2.068
    V_T     =    sqrt((2 * GAMMA * GAS_CONSTANT * CHAMBER_TEMPERATURE) / (GAMMA + 1))
    V_E     =    sqrt(T_T * (1 - PR_2))

    if  (DESIRED_MASS_FLOW_RATE  ==  0):
        DESIRED_MASS_FLOW_RATE  =  DESIRED_THRUST / V_E
    elif  (DESIRED_THRUST  ==  0):
        DESIRED_THRUST  =  DESIRED_MASS_FLOW_RATE / V_E

    T_E    =   CHAMBER_TEMPERATURE * (P_0 / CHAMBER_PRESSURE) ** ((GAMMA - 1) / GAMMA)
    A_E    =   sqrt(GAMMA * GAS_CONSTANT * T_E)
    M_e    =   V_E / A_E
    _RTOD_   =  180 / pi
    _DTOR_   =  pi / 180

    _A_     =  sqrt((GAMMA + 1) / (GAMMA - 1))
    _B_     =  (GAMMA - 1) / (GAMMA + 1)
    V__PM   =  lambda x: _A_ * atan(sqrt(_B_ * (x ** 2 - 1))) - atan(sqrt(x ** 2 - 1))


    T_MAX   = 0.5 * V__PM(M_e) * _RTOD_
    _DT_      = (90 - T_MAX) - round(90 - T_MAX)
    T_0,M,RR,LR,SL,P  =  [],[0.0000],[0.0000],[0.0000],[0.0000],[0.0000]
    T_0.append(_DT_*_DTOR_)
    n   =  T_MAX * 2

    for m in range(1, int(n) + 1):
        T_0.append((_DT_ + m) * _DTOR_)


        X_INT  = [ 1, 1.01 * M_e ]
        _FUNC_   = lambda x: T_0[m] - V__PM(x)
        M.append(brentq(_FUNC_, X_INT[0], X_INT[1]))
        P.append(0  +  THROAT_RADIUS * tan(T_0[m]))


        RR.append(-THROAT_RADIUS / P[m])


        LR.append(tan(T_0[m] + asin(1 / M[m])))
        SL.append(-RR[m])


    P.pop(0)
    l  =  len(P)
    for j in range(0,l):
        P1  =  [0,THROAT_RADIUS]
        P2  =  [P[j], 0]
    LR.pop(0)
    SL.pop(0)
    RR.pop(0)
    F   = RR[m - 1]
    x,y = [],[]
    for c in range(0,len(P) - 1):
        x.append((THROAT_RADIUS + SL[c] * P[c]) / (SL[c] - F))
        y.append(F * x[c] + THROAT_RADIUS)

    xD,yD,xC,yC,s,b = [],[],[],[],[],[]

    #Converging Section
    stepSize = 0.05 
    cNr = CHAMBER_RADIUS - THROAT_RADIUS
    #Before N
    t = 1*pi
    while t < 1.5*pi:
        xC.append(cNr*cos(t)+cNr)
        yC.append(1.5*cNr*sin(t)+cNr+THROAT_RADIUS)
        t += stepSize
    xCi = len(xC)
    yCi = len(yC)

    #After N
    while t < 2*pi:
        xC.append((0.382*THROAT_RADIUS)*cos(t)+(xC[xCi-1]-((0.382*THROAT_RADIUS)*cos(1.5*pi))))
        yC.append((0.382*THROAT_RADIUS)*sin(t)+(yC[yCi-1]-((0.382*THROAT_RADIUS)*sin(1.5*pi))))
        t += stepSize

    #Diverging Section
    _TM_    =   T_MAX  *  _DTOR_
    xD.append((THROAT_RADIUS + SL[0] * P[0]) / (SL[0] - tan(_TM_)))
    yD.append((tan(_TM_) * xD[0] + THROAT_RADIUS))


    _DTW_  =  tan(_TM_) / (len(P) - 1)
    s.append(tan(_TM_))
    b.append(THROAT_RADIUS)
    for k in range(1, len(P) - 1):
        s.append(tan(_TM_) - (k) * _DTW_)
        b.append(yD[k - 1] - s[k] * xD[k - 1])
        xD.append((b[k] + SL[k] * P[k]) / (SL[k] - s[k]))
        yD.append(s[k] * xD[k] + b[k])
        
    #YOU'LL GET A GRAPH WHICH REPRESENTS THE BELL NOZZLE CONTOUR

    xf  = (b[len(b) - 1] + SL[len(SL) - 1] * P[len(P) - 1]) / SL[len(SL) - 1]
    yf  =  b[len(b) - 1]
    xD = [0]  + xD
    yD = [THROAT_RADIUS] + yD
    RTHROAT  =  THROAT_RADIUS
    REXIT    =  yD[len(yD) - 1]
    AR       =  (RTHROAT / REXIT) ** 2

    #Slope of diverging nozzle
    sC,sD,sT = [],[],[]

    #Converging Nozzle derivative
    for i in range(len(xC)-1):
        dxC = xC[i+1] - xC[i]
        dyC = yC[i+1] - yC[i]
        if dyC == 0 or dxC == 0:
            dyC = 100
            dxC = 1 
        sC.append(dyC/dxC)

    #Diverging Nozzle derivative
    for i in range(len(xD)-1):
        dxD = xD[i+1] - xD[i]
        dyD = yD[i+1] - yD[i]
        if dyD == 0 or dxD == 0:
            dyD = 100
            dxD = 1 
        sD.append(dyD/dxD)

    for i in range (len(sC)-30, len(sC)):
        sT.append(sD[0] - sC[i])

    q = len(xC) - 30 + sT.index(min(sT, key=abs))
    del xC[q:] 
    del yC[q:]

    xT,yT,zT = [],[],[]
    for i in range(len(xC) + len(xD)):
        if i < len(xC):
            h = 0
            while h < 2*pi:
                xT.append(xC[i])
                yT.append(yC[i]*sin(h))
                zT.append(yC[i]*cos(h))
                h += stepSize
        else:
            h =0
            xValue = xD[i-len(xC)]-(xD[0]-xC[len(xC)-1])
            sRadius = yD[i-len(yC)]-(yD[0]-yC[len(yC)-1])
            while h <2*pi:
                xT.append(xValue)
                yT.append(sRadius*sin(h))
                zT.append(sRadius*cos(h))
                h += stepSize

    #making data into xyz Sarray
    vertices = np.array([[xT[0],yT[0],zT[0]]])
    for i in range(1, len(xT)):
        vertices = np.append(vertices,[[xT[i],yT[i],zT[i]]], axis=0)

    #adding faces
    faces = np.array([[125,0,126]])
    for i in range(len(xT)):
        if (i < 126):
            faces = np.append(faces,[[i,i+1,i+126]],axis=0)   
        else:
            if (i >= len(xT)-126):
                if (i == len(xT)-1):
                    break
                else:
                    faces = np.append(faces,[[i,i-125,i+1]],axis=0)
            else:
                faces = np.append(faces,[[i,i+1,i+126]],axis=0)
                faces = np.append(faces,[[i,i-125,i+1]],axis=0)

    shape = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
            for j in range(3):
                    shape.vectors[i][j] = vertices[f[j]]

    shape.save("generatedNozzle.stl")

