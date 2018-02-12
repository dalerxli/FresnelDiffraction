#!/usr/bin/python3.6
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 09:10:37 2017

Fresnel Diffraction Approximation Tools for Bristol University Physics Laboratory.

@author: S Nickolls
"""

import numpy as np
import matplotlib.pyplot as plt
import time

class FresnelDiffractionSource():
    #initialize variable arguments on creation of Source object
    def __init__(self, wavelength, field_strength, ap_shape = [], ap_limits = []):
        self.wavelength = wavelength
        self.field_strength = field_strength
        self.ap_shape = ap_shape
        self.ap_limits = ap_limits
    
    #performs Simpsons integral approximation of function, over variables bounds bounds[] with N (even) approximation terms
    #input FresnelDiffractionSource, no of terms, function to integrate, integration limits as array, returns complex no.    
    def integration(self, N, function, bounds = []):
        if (N/2 != int(N/2)):
            return 'ERR'
        
        width = (bounds[1] - bounds[0]) / (N)
        
        integral = counter = 0
        
        x = bounds[0]
        integral = function(x)
        
        for i in range(1,N):
            x += width
            if (counter == 0):
                integral += 4 * function(x)
                counter = 1
            else:
                integral += 2 * function(x)
                counter = 0
        
        integral += function(bounds[1])        
        integral *= (width / 3)
        #print("<DEBUG> integral:", integral)
        
        return integral
    
    #performs double integral using Simpsons Approximation for a separable function f(x,y) = xfunction * yfunction
    #input FresnelDiffractionSource, no of terms, xfunction, yfunction, xlimits (array of functions of a single variable),
    #constant limits over the y axis, returns single complex no., the value of the double integral
    def double_integration(self, N, xfunction, yfunction, xlimits = [], ylimits = []):
        if (N/2 != int(N/2)):
            return 'ERR' 
        
        #wrapper function to adapt the limits for varying values of y
        def limit_function(y):
            x_lower = xlimits[0](y)
            x_upper = xlimits[1](y)
            limits = (x_lower, x_upper)
            
            return self.integration(N, xfunction, limits) * yfunction(y)

        integral = self.integration(N, limit_function, ylimits)
        
        return integral

    #calculates the complex exponential phase term, using values for wavelength, perpendicular distance from screen,
    #single dimensional coordinate on the screen, intention to integrate across aperture values ap_x
    def phase_function(self, distance, screen_x, ap_x):
        k = 2 * np.pi / self.wavelength
        return np.exp(k * 1j * (screen_x - ap_x)**2 / (2*distance))

    #produces a 1 dimensional plot of the Fresnel Diffraction Pattern at distance from the source, on a pointsxpoints grid
    #across the screen_limits. returns the value of plt.show(), which prints the graph. input initialized
    #FresnelDiffractionSource with non-zero attributes
    def single_interference(self, distance, points, screen_limits = []):
        xvals = np.zeros(points, dtype=np.complex)
        yvals = np.zeros(points, dtype=np.complex)
        
        width = (screen_limits[1] - screen_limits[0]) / points
        
        for i in range(0, points):
            
            screen_x = 0 + i * width
            
            #wrapper function for phase_function allowing variable screen_x input
            def phase(ap_x):
                return self.phase_function(distance, screen_x, ap_x)
            
            integral = self.integration(100, phase, [self.ap_limits[0], self.ap_limits[1]])
            
            xvals[i] = screen_x
            yvals[i] = abs(integral * integral.conjugate())
            #print(xvals[i], integral, yvals[i])
            
        plt.plot(xvals, yvals)
        
        return plt.show()
    
    #produces a 2D image of the fresnel diffraction at a distance from the source, on a pointsxpoints grid across 
    #the screen_limits. input initialized FresnelDiffractionSource with non-zero attributes, returns value of plt.show()
    def double_interference(self, distance, points, screen_limits = []):
        width = (screen_limits[1] - screen_limits[0]) / points
        
        ivals = np.zeros((points, points))   
        
        for i in range(0, points):
            screen_x = screen_limits[0] + i * width
            
            for j in range(0, points):  
                screen_y = screen_limits[0] + j * width
                
                #wrapper functions for phase_function allowing screen_x and screen_y to vary
                def xphase(ap_x):
                    return self.phase_function(distance, screen_x, ap_x)
                def yphase(ap_y):
                    return self.phase_function(distance, screen_y, ap_y)
                
                phase_integral = self.double_integration(50, xphase, yphase, self.ap_shape, self.ap_limits)
                phase_integral *= self.field_strength / (distance * self.wavelength)
                
                ivals[i,j] = (abs(phase_integral*phase_integral.conjugate()))**2 
           
        plt.imshow(ivals) 
        
        return plt.show()

#formatting functions because I'm used to C
def printf(string):
    return print(string, end='')
def selection():
    return input("Selection: ")

#main function handling user input and menu sequence, returns 0 on success.
def main():
    printf("Fresnel Diffraction Simulator by Sebastian Nickolls\nPhysics Undergraduate Laboratory, "+
           "University of Bristol 2018\n")
    choice = 'y'
    subchoice = ''
    
    printf("==========Menu==========\nEnter letter to continue.\na) graph 1d diffraction patterns\nb) image 2d "+
            "diffraction patterns.\nq) quit\n")
    
    while (choice != 'q'):
        choice = selection()
        
        if (choice == 'q'):
            return 0
         
        elif (choice == 'a' or choice == 'b'):
            printf("This script has default values for wavelength and electric field strength in standard units.\n")    
            printf("These values are 5E-7m, and 1N/C respectively.\n")
      
            if (choice == 'a'):
                printf("The light source pattern will be projected onto a 2E-4m length screen, resolution 1E-6m, "+
                       "from an aperture of roughly the same dimension.\n")
            if (choice == 'b'): 
                printf("The light source pattern will be projected onto a 2E-4mx2E-4m square screen, resolution 1E-6m, "+
                       "from an aperture of roughly the same dimensions.\n")
                        
            printf("The image will be generated and shown consecutively at distances z, 2z and 3z from the screen.\n"+
                   "The value z is default to the aperture shape, but may be changed in the source code. "+
                   "The resolution may also be change in the source.\n")
                
            screen_limits = [-1E-4, 1E-4]
            points = 75
            
            printf("Use default source values? [y/n]\n")
            subchoice = ''
            
            #set wavelength and field values
            while (subchoice != 'y' or 'n'):   
                subchoice = selection()
                
                if (subchoice == 'y'):
                    wavelength = 5E-7
                    field = 1
                    break
                elif (subchoice == 'n'):
                    printf("Enter wavelength value.")
                    wavelength = float(selection())
                    printf("Enter electric field value.")
                    field = float(selection())
                    break
                elif (subchoice != 'y' or subchoice != 'n'):
                    printf("Invalid Choice.\n")
            
            subchoice = 's'
            ap_limits = [-1E-4, 1E-4]
            z_interval = 5E-3
            
            if (choice == 'b'):
                subchoice = ''
                printf("Please choose aperture shape.\nc) circle\ns) square\nt) triangle (equilateral)\n")
                
                #set aperture shape
                while (subchoice != 'c' or subchoice != 's' or subchoice != 't'):
                    subchoice = selection()
                    
                    if (subchoice == 'c' or subchoice == 's'):
                        ap_limits = [-1E-4, 1E-4]
                        z_interval = 5E-3
                        break
                    if (subchoice == 't'):
                        ap_limits = [-0.5E-4,1E-4]
                        z_interval = 2.5E-3
                        break
                    else:
                        printf("Invalid Choice.\n")
                
                #defines the aperture lower limit based on the choice of shape
                def ap_lower(y, shape):
                    if (shape == 'c'): return -np.sqrt((1E-4)**2 - y**2) 
                    if (shape == 's'): return -1E-4 
                    if (shape == 't'): return (y - 1E-4) / np.sqrt(3) 
                    
                #defines the aperture upper limit based on the choice of shape
                def ap_upper(y, shape):
                    if (shape == 'c'): return np.sqrt((1E-4)**2 - y**2)
                    if (shape == 's'): return 1E-4
                    if (shape == 't'): return (1E-4 - y) / np.sqrt(3)
            
            #define aperture limits based on choice of shape
            def ap_low(y):
                return ap_lower(y, subchoice)
            def ap_high(y):
                return ap_upper(y, subchoice)
            
            #initialize a light source
            Source = FresnelDiffractionSource(wavelength, field, [ap_low, ap_high], ap_limits)
            
            #complete calculations
            for i in range(1,4):
                init_time = time.time()
                z = i * z_interval
                print("Distance from screen = %fm" % z)
                if (choice == 'a'): 
                    FresnelDiffractionSource.single_interference(Source, z, points, screen_limits)
                elif (choice == 'b'):
                    FresnelDiffractionSource.double_interference(Source, z, points, screen_limits)
                final_time = time.time()
                print("Time Elapsed: %fs\n" % (final_time - init_time))
            
            printf("Execution complete.\n")
            
            break
        
        else:
            printf("Invalid Choice.\n")
    
    printf("Exiting...\n")
    return 0

#execute script
if (main() != 0):
    print("Something went horribly wrong!")
            
