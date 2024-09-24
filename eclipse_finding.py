# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 10:42:21 2024

@author: lbnc
"""
import numpy as np
from skyfield.api import load
from skyfield.functions import dots, length_of, angle_between

def rotate(axis, angle, vector):
    # https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    # Can this be more elegant? Sure!
    # Do I care? No!
    s, c = np.sin(np.radians(angle)), np.cos(np.radians(angle))
    x = (axis[0]*axis[0]*(1 - c) + c)*vector[0] + \
        (axis[0]*axis[1]*(1 - c) - axis[2]*s)*vector[1] + \
            (axis[0]*axis[2]*(1 - c) + axis[1]*s)*vector[2]
    y = (axis[0]*axis[1]*(1 - c) + axis[2]*s)*vector[0] + \
        (axis[1]*axis[1]*(1 - c) + c)*vector[1] + \
            (axis[1]*axis[2]*(1 - c) - axis[0]*s)*vector[2]
    z = (axis[0]*axis[2]*(1 - c) - axis[1]*s)*vector[0] + \
        (axis[1]*axis[2]*(1 - c) + axis[0]*s)*vector[1] + \
            (axis[2]*axis[2]*(1 - c) + c)*vector[2]
    return np.array([x, y, z])

def cross(a, b):
    return np.array([a[1]*b[2] - a[2]*b[1], 
                     a[2]*b[0] - a[0]*b[2], 
                     a[0]*b[1] - a[1]*b[0]])

def find_eclipse(eclipsee_position, eclipsee_R, 
                 eclipser_position, eclipser_R, 
                 screen_position, screen_R):
    # position of screen
    screen_position_km = screen_position.position.km # in km
    # position of eclipsee
    eclipsee_position_km = screen_position_km + eclipsee_position.position.km # in km
    # position of eclipser
    eclipser_position_km = screen_position_km + eclipser_position.position.km # in km
    # intersection of centerline
    o = eclipsee_position_km
    # direction of line (eclipsee rim to eclipser rim)
    u = eclipser_position_km - eclipsee_position_km
    # center of sphere (screen)
    t = screen_position_km
    # radius of sphere (screen)
    r = screen_R
    # coefficients
    a = dots(u, u)
    b = 2*dots(u, o-t)
    c = dots(o-t, o-t) - r**2
    # doesn't hit the screen
    if b**2 - 4*a*c < 0:
        return None, None, None
    # Bingo!
    else:
        s = (-b - np.sqrt(b**2 - 4*a*c))/2/a
        # intersection position in screen-centric coordinates
        intersct_cpos = o + s*u
    # distance eclipsee -> eclipser
    d = length_of(eclipser_position_km - eclipsee_position_km) # in km
    # distance eclipser -> screen
    d2 = length_of(screen_position_km - eclipser_position_km) # in km
    # angles
    # minimum angle of ray that is boundary/edge for umbra
    a1 = np.arctan2(eclipsee_R - eclipser_R, d)
    # maximum angle of ray that is boundary/edge for penumbra
    a2 = np.arctan2(eclipsee_R + eclipser_R, d)
    # maximum angle of ray that is boundary/edge for umbra
    a3 = np.arctan2(eclipsee_R + eclipser_R - d2*np.tan(a1), d + d2)
    # shadow region radius
    shadow_R = eclipser_R - d2*np.tan(a1)
    # vector from eclipsee center to eclipser center
    along_dir = eclipser_position_km - eclipsee_position_km
    # normalize
    along_dir /= length_of(along_dir)
    # vector perpendicular to along_dir
    perp_dir = cross(eclipsee_position_km, along_dir)
    # normalize
    perp_dir /= length_of(perp_dir)
    # radius vector of eclipsee
    eclipsee_radvec = eclipsee_R*perp_dir
    # radius vector of eclipser
    eclipser_radvec = eclipser_R*perp_dir
    # points along surface of eclipsee and eclipser
    # assuming circular bodies
    eclipsee_rim = eclipsee_position_km + \
        np.array([rotate(along_dir, a, eclipsee_radvec) for a in range(0, 360, 20)])
    eclipser_rim = eclipser_position_km + \
        np.array([rotate(along_dir, a, eclipser_radvec) for a in range(0, 360, 20)])
    # for each combination of points on the rims, 
    # calculate whether the connecting line intersects with observer
    # again assuming a circular observer
    intersct_pos = np.zeros((eclipsee_rim.shape[0], eclipser_rim.shape[0], 3))
    intersct_ang = np.zeros((eclipsee_rim.shape[0], eclipser_rim.shape[0]))
    intersct_cat = np.zeros((eclipsee_rim.shape[0], eclipser_rim.shape[0]))
    # fudge factor
    eps = 1e-6
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
            a = dots(u, u)
            b = 2*dots(u, o-t)
            c = dots(o-t, o-t) - r**2
            # doesn't hit the screen
            if b**2 - 4*a*c < 0:
                pass
            # Bingo!
            else:
                # distance from eclipsee rim to intersection point on screen
                # only take the minus here, because that is the shortest
                # I think
                s = (-b - np.sqrt(b**2 - 4*a*c))/2/a
                # intersection position in screen-centric coordinates
                intersct_pos[x,y,:] = o + s*u
                # find category of intercept, umbra or penumbra
                intersct_ang[x,y] = angle_between(u,along_dir)
                # check for antumbra
                if a1 > 0 and eclipser_R/np.tan(a1) < d2:
                    intersct_cat[x,y] = 3
                elif intersct_ang[x,y] >= (1-eps)*a1 and intersct_ang[x,y] < a3:
                    # umbra
                    intersct_cat[x,y] = 2
                elif intersct_ang[x,y] >= a3 and intersct_ang[x,y] <= (1+eps)*a2:
                    # penumbra
                    intersct_cat[x,y] = 1
                else:
                    # let's hope we never get here
                    intersct_cat[x,y] = 0
    return intersct_cpos, intersct_pos, intersct_cat


ts = load.timescale()
t0 = ts.utc(2024, 4, 8, 18, 0, 0)
#t0 = ts.utc(2017, 8, 21, 18, 30, 0)
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

# intersection coordinates in B(I)CRS coordinates
# because everything is relative to the screen position, 
# which is in ICRF coords
intersct_cpos, intersct_pos, intersct_cat = find_eclipse(
    eclipsee_position, eclipsee_R, 
    eclipser_position, eclipser_R, 
    screen_position, screen_R)
# plot result in ICRF
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(7,7),dpi=150)
labels = ["penumbra", "umbra"]
colors = ["tab:gray", "k"]
for i in [1,2]:
    idx = intersct_cat == i
    ax.scatter(intersct_pos[idx,0], intersct_pos[idx,1], 
               c=colors[i-1], label=labels[i-1])
ax.scatter(intersct_cpos[0], intersct_cpos[1], 
           c="tab:orange", label="center")
plt.axis("equal")
plt.legend()

# convert to geocentric latitude/longitude
from skyfield.positionlib import Geocentric
from skyfield.api import wgs84

# geocentric ECEF position of centerline intersection
cposition = Geocentric(
    intersct_cpos - screen_position.position.km, 
    velocity_au_per_d=[0,0,0], 
    t=t0, 
    center=screen_position.target)
cp = wgs84.subpoint(cposition)
# This is slightly off compared to 
# https://eclipse.gsfc.nasa.gov/SEpath/SEpath2001/SE2024Apr08Tpath.html
# Not sure why
print(cp.latitude.degrees, cp.longitude.degrees)

# all intersections of the traced rimlines
g = np.empty_like(intersct_cat, dtype=object)
for x in range(intersct_pos.shape[0]):
    for y in range(intersct_pos.shape[1]):
        position = Geocentric(
            intersct_pos[x,y,:] - screen_position.position.km, 
            velocity_au_per_d=[0,0,0], 
            t=t0, 
            center=screen_position.target)
        g[x,y] = wgs84.subpoint(position)
# This seems belony, as the umbra extends by more then 10deg, 
# from about 14degN to over 25degN
# from shadow_R above (and NASA) the diameter of the shadow region should be
# about 100 km, i.e., about 1 deg in latitude
for gp in g[intersct_cat == 2]:
    print(gp.latitude.degrees, gp.longitude.degrees)
