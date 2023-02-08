# Rock Matching Functions
# Maddie Schwarz 1/26/23
# Code modified from: https://github.com/DREAMS-lab/Rock-Matching
#%%
import warnings
warnings.filterwarnings('ignore')
import geopandas as gpd
import math
import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt
from matplotlib.patches import ConnectionPatch
import os
from scipy.spatial.distance import euclidean
import shapely
from shapely.geometry import LineString
from shapely import affinity

#FUNCTIONS

## GEOMETRIC
def get_ecc(area_polygon): # calc eccentricity
    points = np.array(area_polygon)

    small_latwise = np.min(points[points[:, 0] == np.min(points[:, 0])], 0)
    small_lonwise = np.min(points[points[:, 1] == np.min(points[:, 1])], 0)
    big_latwise = np.max(points[points[:, 0] == np.max(points[:, 0])], 0)
    big_lonwise = np.max(points[points[:, 1] == np.max(points[:, 1])], 0)
    distance_lat = euclidean(big_latwise, small_latwise)
    distance_lon = euclidean(big_lonwise, small_lonwise)
    if distance_lat >= distance_lon:
        major_axis_length = distance_lat
        minor_axis_length = distance_lon
    else:
        major_axis_length = distance_lon
        minor_axis_length = distance_lat
    a = major_axis_length/2
    b = minor_axis_length/2
    ecc = np.sqrt(np.square(a)-np.square(b))/a
    return ecc

def pp_compactness(geom): # Compactness Polsby-Popper
    p = geom.length
    a = geom.area    
    return (4*pi*a)/(p*p)

def convexity(geom): # Convexity Polsby-Popper
    p = geom.convex_hull.length
    a = geom.length    
    return (p/a)

def solidity(geom): # Calc solidity (Polsby-Popper)
    convex_area = geom.convex_hull.area
    area = geom.area    
    return (area/convex_area)

def getAngle(geom): # Primary axis orientation
    g = geom
    a = g.minimum_rotated_rectangle
    l = a.boundary
    coords = [c for c in l.coords]
    segments = [LineString([a, b]) for a, b in zip(coords,coords[1:])]
    longest_segment = max(segments, key=lambda x: x.length)

    p1, p2 = [c for c in longest_segment.coords]
    angle = math.degrees(math.atan2(p2[1]-p1[1], p2[0]-p1[0]))
    
    return angle



## SIMILARITY MEASURES
# functions sourced from this code: https://www.analyticsvidhya.com/blog/2021/06/nlp-answer-retrieval-from-document-using-ts-ss-similarity-python/#:~:text=TS%2DSS%20computes%20the%20similarity,the%20difference%20between%20their%20magnitudes.

def VectorSize(vec) : # calculates the length of a vector (ref: https://mathinsight.org/vectors_cartesian_coordinates_2d_3d#:~:text=Using%20the%20Pythagorean%20Theorem%2C%20we,21%2Ba22.)
    return math.sqrt(sum(math.pow(v,2) for v in vec))

def InnerProduct(vec1, vec2) : #inner product is a generalization of the dot product; way to multiply vectors together to get scalar result
    return sum(v1*v2 for v1,v2 in zip(vec1,vec2))

def Cosine(vec1, vec2) : 
    result = InnerProduct(vec1,vec2) / (VectorSize(vec1) * VectorSize(vec2))
    return result

def Theta(vec1, vec2) :# inverse of Cosine is our angle theta
    CosTrun = float(f'{Cosine(vec1,vec2):.8f}') # temp fix (1/19/23): truncate calculation of Cos() to prevent floating point err w/ arcos()
    return math.acos(CosTrun) + math.radians(10)

def Euclidean(vec1, vec2) :# calc euclid. dist between 2 vectors
    return math.sqrt(sum(math.pow((v1-v2),2) for v1,v2 in zip(vec1, vec2)))

def Triangle(vec1, vec2) :
    theta = math.radians(Theta(vec1,vec2))
    return (VectorSize(vec1) * VectorSize(vec2) * math.sin(theta)) / 2

def Magnitude_Difference(vec1, vec2) :
    return abs(VectorSize(vec1) - VectorSize(vec2))

def Sector(vec1, vec2) :
    ED = Euclidean(vec1, vec2)
    MD = Magnitude_Difference(vec1, vec2)
    theta = Theta(vec1, vec2)
    return math.pi * math.pow((ED+MD),2) * theta/360

def TS_SS(vec1, vec2) :
    return Triangle(vec1, vec2) * Sector(vec1, vec2)

# Visualization
def get_rock_and_line(method, df1, df2):
    df1['matching_rock'] = None
    for i in range(len(df1)):
        df1['matching_rock'][i] = df2['centroid'][df1['nearest_matching_polygon'+method][i]]
    df1['line'] = df1.apply(lambda x: LineString([x['centroid'], x['matching_rock']]), axis=1)
    return df1 
        #if df1['nearest_matching_polygon'+method][i] != None:
            #df1['matching_rock'][i] = df2['centroid'][df1['nearest_matching_polygon'+method][i]]
        #else:
           # continue
        #df1['line'] = df1.apply(lambda x: LineString([x['centroid'], x['matching_rock']]), axis=1)
    #return df1 

# 1/26/23: Note-- this function is designed to check
# synthetic rock datasets, not real dataset (where rock order in df may differ)
def check_matches(method,df1,df2):
    correct_m = 0
    incorrect_m = 0
    print('checking matches:')
    for i in range(len(df1)):
        print(i, df1['nearest_matching_polygon'+method][i])
        if df1['nearest_matching_polygon'+method][i] != i:
            #print('incorrect match')
            incorrect_m = incorrect_m + 1
        else:
            #print('correct')
            correct_m = correct_m + 1
    print('Similarity measure used: ', method)
    print('Number of rocks in TimeFrame 1: ', len(df1))
    print('Number of rocks in TimeFrame 2: ', len(df2))
    print('Correct matches:', correct_m)
    print('Incorrect matches:', incorrect_m)

    #%%
#vec1 = [0]
#vec2 = [3]

#print(TS_SS(vec1, vec2) )

# %%
#print(math.radians(Theta(vec1,vec2)))

# %%
