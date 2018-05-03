#!/usr/bin/env python3

"""
Shuttle Orbit Simulations
03/05/18
Sebastian Nickolls, University of Bristol Undergraduate Laboratory
Usage: /path/to/python3 /path/to/orbit.py <interval> <duration> <mass> <t0> <x0> <y0> <vx0> <vy0> <moon>
Intended for use on the command line, if this isn't working please pass the command line arguments through spyder using
Run > Configuration per file > Command line options:
If this is done incorrectly the script will return the usage message.
"""

import sys
import math

class RKApproximation():
    #calculate all sets of kvalues and return them as a vector
    def compute_kvals(data, pfunc, vfunc, interval):

        k0 = (0,0,0,0)

        #step through all the kvalues for each variable (4 per step, 16 total)
        k1 = RKApproximation.kstep(data, k0, pfunc, vfunc, 0)
        k2 = RKApproximation.kstep(data, k1, pfunc, vfunc, interval)
        k3 = RKApproximation.kstep(data, k2, pfunc, vfunc, interval)
        k4 = RKApproximation.kstep(data, k3, pfunc, vfunc, interval*2)
 
        kvector = (k1, k2, k3, k4)
 
        return kvector

    #calculate the next set of k values from the last set
    def kstep(data, kN, pfunc, vfunc, interval):
        intermediate = [0,0,0,0]
    
        for i in range(0, len(data)):
            intermediate[i] = data[i] + (kN[i] * interval/2)

        #calculate k value for each variable
        kx = pfunc(intermediate, 0)
        ky = pfunc(intermediate, 1)
        kvx = vfunc(intermediate, 0)
        kvy = vfunc(intermediate, 1)

        return (kx, ky, kvx, kvy)

    #sum up the kvalues in a kvector
    def sum_kvals(kvector, index, interval):
        i = index
        return (kvector[0][i] + 2*(kvector[1][i]+kvector[2][i]) + kvector[3][i]) * interval / 6

#Coordinate class has intrinsic positions, initial positions and a reset option
class Coordinate():
    def __init__(self, init_time, init_pos, init_vel):
        self.time = self.time0 = init_time
        self.x = self.x0 = init_pos[0]
        self.y = self.y0 = init_pos[1]
        self.vx = self.vx0 = init_vel[0]
        self.vy = self.vy0 = init_vel[1]

    def reset(self):
        self.time = self.time0
        self.x = self.x0
        self.y = self.y0
        self.vx = self.vx0
        self.vy = self.vy0

#Shuttle inherits vectors from Coordinate and allows updating the values with given position and velocity derivatives
class Shuttle(Coordinate):
    def __init__(self, mass, time0, pos0, vel0, p_derivative, v_derivative):
        super(Shuttle, self).__init__(time0, pos0, vel0)
        self.p_deriv = p_derivative
        self.v_deriv = v_derivative
        self.mass = mass

    def update(self, interval):
        data = (self.x, self.y, self.vx, self.vy)

        #retrieve kvectors for given interval
        kvect = RKApproximation.compute_kvals(data, self.p_deriv, self.v_deriv, interval)
		
        #update variables
        self.time += interval
        self.x += RKApproximation.sum_kvals(kvect, 0, interval)
        self.y += RKApproximation.sum_kvals(kvect, 1, interval)
        self.vx += RKApproximation.sum_kvals(kvect, 2, interval)
        self.vy += RKApproximation.sum_kvals(kvect, 3, interval)
    
        return 0

def usage():
    print("========================================= Shuttle Flight Simulator =========================================\n"+
          "============================================ Sebastian Nickolls ============================================\n"+
          "[*] Usage: /path/to/python3 /path/to/orbit.py <interval> <duration> <mass> <t0> <x0> <y0> <vx0> <vy0> <moon>\n"+
          "[*] Set <moon> to 1 to include the moon, 0 to ignore\n"
          "[*] Recommended use: python3 orbit.py <args> > data.csv\n"+
          "[*] For Windows: C:\\path\\to\\python3.exe orbit.py <args> > data.csv")
    return -1

def main():
    #constant variables
    G, Me, Mm, p_moon = 6.67E-11, 5.97E24, 7.35E22, 384E6
    moon_radius, earth_radius = 1737000, 6371000
    
    #position and velocity derivatives, with dimension index
    #choose whether include the moon or not
    def pfunc_earth(data, index):
        return data[2+index]
    
    def vfunc_earth(data, index):
        root = (data[0]**2 + data[1]**2)**(3/2)    
        return -1*G*Me*data[index]/root
    
    def pfunc_moon(data,index):
        return data[2+index]
    
    def vfunc_moon(data, index):
        earth = vfunc_earth(data, index)
        
        mroot = ( (data[0] - p_moon)**2 + data[1]**2) ** (3/2)
        
        if(index == 0):
            return earth - Mm*G*(data[0]-p_moon)/mroot
        else:
            return earth - Mm*G*(data[1])/mroot
     
    #script will not run if wrong amount of arguments given
    if len(sys.argv) != 10:
        return usage()
    
    #retrieve values from sys.argv (NOTE HEAVILY DUCKTYPED NO ERROR CHECKING)
    interval = float(sys.argv[1])
    duration = float(sys.argv[2])
    mass = float(sys.argv[3])
    t0 = float(sys.argv[4])
    x0 = float(sys.argv[5])
    y0 = float(sys.argv[6])
    vx0 = float(sys.argv[7])
    vy0 = float(sys.argv[8])
    moon = float(sys.argv[9])
    
    if(moon):
        pfunc = pfunc_moon
        vfunc = vfunc_moon
    else:
        pfunc = pfunc_earth
        vfunc = vfunc_earth
    try:
        S = Shuttle(mass, t0, (x0, y0), (vx0, vy0), pfunc, vfunc)
        
        #dump in csv format, pipe script into .csv file for best results
        print("time, x, y, vx, vy, KE, PE (w.r.t. Earth), E")
        while (S.time < duration):
            S.update(interval)
            
            #extra energy interpretations
            KE = 0.5 * S.mass * (S.vx**2 + S.vy**2)
            PE = -1*G*Me*S.mass / math.sqrt(S.x**2 + S.y**2)
            E = KE + PE
             
            print("%f, %f, %f, %f, %f, %f, %f, %f" % (S.time,S.x,S.y,S.vx, S.vy, KE, PE, E))
            dist_from_earth = math.sqrt(S.x**2 + S.y**2)
            
            #checking for crashes, return 1 if so
            if (dist_from_earth < earth_radius):
                return 1
            if (moon):
                dist_from_moon = math.sqrt((S.x-p_moon)**2 + S.y**2)
                if (dist_from_moon < moon_radius):
                    return 1

    except:
        print("[!] Caught an error in calculations, please check command line arguments!")
        return -1
            
    return 0

if __name__ == "__main__":
    main()
