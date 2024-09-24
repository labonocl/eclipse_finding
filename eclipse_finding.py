# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 10:42:21 2024

@author: lbnc
"""
import numpy as np
from skyfield.api import load

def rotate(axis, angle, vector):
    # https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    u = axis/np.linalg.norm(axis)
    s, c = np.sin(np.radians(angle)), np.cos(np.radians(angle))
    x = (u[0]*u[0]*(1 - c) + c)*vector[0] + (u[0]*u[1]*(1 - c) - u[2]*s)*vector[1] + (u[0]*u[2]*(1 - c) + u[1]*s)*vector[2]
    y = (u[0]*u[1]*(1 - c) + u[2]*s)*vector[0] + (u[1]*u[1]*(1 - c) + c)*vector[1] + (u[1]*u[2]*(1 - c) - u[0]*s)*vector[2]
    z = (u[0]*u[2]*(1 - c) - u[1]*s)*vector[0] + (u[1]*u[2]*(1 - c) + u[0]*s)*vector[1] + (u[2]*u[2]*(1 - c) + c)*vector[2]
    return np.array([x, y, z])

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def cross(a, b):
    return np.array([a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]])

def find_eclipse(eclipsee_position, eclipsee_R, 
                 eclipser_position, eclipser_R, 
                 screen_position, screen_R):
    # position of screen
    screen_position_km = screen_position.position.km # in km
    # position of eclipsee
    eclipsee_position_km = screen_position_km + eclipsee_position.position.km # in km
    # position of eclipser
    eclipser_position_km = screen_position_km + eclipser_position.position.km # in km
    # distance eclipsee -> eclipser
    d = np.linalg.norm(eclipser_position_km - eclipsee_position_km) # in km
    # distance eclipser -> screen
    d2 = np.linalg.norm(screen_position_km - eclipser_position_km) # in km
    # angles
    # minimum angle of ray that is boundary/edge for umbra
    a1 = np.arctan2(eclipsee_R - eclipser_R, d)
    # maximum angle of ray that is boundary/edge for penumbra
    a2 = np.arctan2(eclipsee_R + eclipser_R, d)
    # maximum angle of ray that is boundary/edge for umbra
    a3 = np.arctan2(eclipsee_R + eclipser_R - d2*np.tan(a1), d + d2)
    # shadow region diameter
    shadow_R = eclipser_R - d2*np.tan(a1)
    # vector from eclipsee center to eclipser center
    along_dir = eclipser_position_km - eclipsee_position_km
    # vector perpendicular to along_dir
    perp_dir = cross(eclipsee_position_km, along_dir)
    # normalize
    perp_dir /= np.linalg.norm(perp_dir)
    # radius vector of eclipsee
    eclipsee_radvec = eclipsee_R*perp_dir
    # radius vector of eclipser
    eclipser_radvec = eclipser_R*perp_dir
    # points along surface of eclipsee and eclipser
    # assuming circular bodies
    eclipsee_rim = np.array([eclipsee_position_km + rotate(along_dir, a, eclipsee_radvec) for a in range(0, 360, 10)])
    eclipser_rim = np.array([eclipser_position_km + rotate(along_dir, a, eclipser_radvec) for a in range(0, 360, 10)])
    # for each combination of points on the rims, 
    # calculate whether the connecting line intersects with observer
    # again assuming a circular observer
    intercpt_pos = np.zeros((eclipsee_rim.shape[0], eclipser_rim.shape[0], 3))
    intercpt_ang = np.zeros((eclipsee_rim.shape[0], eclipser_rim.shape[0]))
    intercpt_cat = np.zeros((eclipsee_rim.shape[0], eclipser_rim.shape[0]))
    # fudge factor
    eps = 1e-4
    for x,eclipsee_pt in enumerate(eclipsee_rim):
        for y,eclipser_pt in enumerate(eclipser_rim):
            # https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection#Calculation_using_vectors_in_3D
            # origin of line (eclipsee rim)
            o = eclipsee_pt
            # direction of line (eclipsee rim to eclipser rim)
            u = eclipser_pt - eclipsee_pt
            # center of sphere (screen)
            t = screen_position_km
            # radius of sphere (screen)
            r = screen_R
            # coefficients
            a = dot(u, u)
            b = 2*dot(u, o-t)
            c = dot(o-t, o-t) - r**2
            # doesn't hit the screen
            if b**2 - 4*a*c < 0:
                pass
            # Bingo!
            else:
                # distance from eclipsee rim to intersection point on screen
                # only take the minus here, because that is the shortest
                d = (-b - np.sqrt(b**2 - 4*a*c))/2/a
                # intersection position in screen-centric coordinates
                intercpt_pos[x,y,:] = o + d*u
                # find category of intercept, umbra or penumbra
                intercpt_ang[x,y] = np.arccos(dot(u, along_dir)/(np.linalg.norm(u)*np.linalg.norm(along_dir)))
                # check for antumbra
                if a1 > 0 and eclipser_R/np.tan(a1) < d2:
                    intercpt_cat[x,y] = 3
                elif intercpt_ang[x,y] >= (1-eps)*a1 and intercpt_ang[x,y] < a3:
                    # umbra
                    intercpt_cat[x,y] = 2
                elif intercpt_ang[x,y] >= a3 and intercpt_ang[x,y] <= (1+eps)*a2:
                    # penumbra
                    intercpt_cat[x,y] = 1
                else:
                    # let's hope we never get here
                    intercpt_cat[x,y] = 0
    return intercpt_pos, intercpt_cat


ts = load.timescale()
t0 = ts.utc(2024, 4, 8, 18, 0, 0)
bodies = load('de421.bsp')
earth, sun, moon = bodies["earth"], bodies["sun"], bodies["moon"]

# body A, the eclipsee
eclipsee, eclipsee_R = sun, 696340
# body B, the eclipser
eclipser, eclipser_R = moon, 1737.4
# body C, the screen
screen, screen_R = earth, 6371

# barycentric position of screen
screen_position = screen.at(t0)
# eclipser position seen from screen
eclipser_position = screen_position.observe(eclipser).apparent()
# eclipsee position seen from screen
eclipsee_position = screen_position.observe(eclipsee).apparent()

intercpt_pos, intercpt_cat = find_eclipse(eclipsee_position, eclipsee_R, 
                                          eclipser_position, eclipser_R, 
                                          screen_position, screen_R)
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7,7),dpi=150)
labels = ["penumbra", "umbra"]
colors = ["tab:gray", "k"]
for i in [1,2]:
    idx = intercpt_cat == i
    ax.scatter(intercpt_pos[idx,0], intercpt_pos[idx,1], c=colors[i-1], label=labels[i-1])
plt.axis("equal")
plt.legend()
