# Calculate the UVW velocities for an incoming asymptotic trajectory from 
# its velocity vector: direction (alpha and delta) and module (vinf)
#
# alpha, delta: ra and dec of direction of approximation (J2000) 
# vinf: velocity of approximation in km/s (positive if approating the Solar system)
#
# U, V, W: output given in km/s
#
#
# Author: Felipe Almeida-Fernandes
# email: felipefer42@gmail.com
#
# DON'T PANIC!
#
# Date: March, 2023
#
#
# Explanation:
# Johnson & Soderblom 1987: Formulas to convert velocity vector in equatorial to UVW
#
# For the incoming trajectory, pmra = pmdec = 0 and rho = -vinf
# Equation (1) becomes:
#
# [ U ]               [ cos(delta) cos(alpha) ]
# [ V ] = -vinf x T x [ cos(delta) sin(alpha) ]
# [ W ]               [ sin(delta)            ]

import numpy as np

# Obtaining the matrix T
# Coordinates of the NGP from Karim & Mamajek (2017)

t0 = 123 * np.pi/180
aNGP = 192.25 * np.pi/180
dNGP = 27.4 * np.pi/180

T_t0   = np.array([[+np.cos(t0), +np.sin(t0),  0],
                   [+np.sin(t0), -np.cos(t0),  0],
                   [          0,           0, +1]])

T_dNGP = np.array([[-np.sin(dNGP),  0, +np.cos(dNGP)],
                   [            0, -1,             0],
                   [+np.cos(dNGP),  0, +np.sin(dNGP)]])

T_aNGP = np.array([[+np.cos(aNGP), +np.sin(aNGP),  0],
                   [+np.sin(aNGP), -np.cos(aNGP),  0],
                   [            0,             0, +1]])

T = np.matmul(T_t0, T_dNGP)
T = np.matmul(T, T_aNGP)


# Get each element to optmize operations

T11 = T[0][0]
T12 = T[0][1]
T13 = T[0][2]

T21 = T[1][0]
T22 = T[1][1]
T23 = T[1][2]

T31 = T[2][0]
T32 = T[2][1]
T33 = T[2][2]


def radecvinf_to_UVW(ra, dec, vinf):
    """
    Estimates the U, V, W galactic velocities for an incoming asymptotic trajectory
    defined by the incoming direction: ra, dec (J2000); and velocity: vinf (km/s)
    
    Inputs can also be arrays (and so will be the outputs)
    
    Parameters
    ----------
    ra : float
        Aproaching right ascension in degrees (J2000)
    dec : float
        Aproaching declination in degrees (J2000)
    vinf : float
        Aproaching velocity, before arriving in the Solar System (km/s)
        Positive for incoming objects
    
    Returns
    -------
    float, float, float
        right-handed U, V, W velocities (km/s). 
        U is positive towards the galactic center
        

    Example
    -------
    
    for 1I/'Oumuamua
    >>> U, V, W = radecvinf_to_UVW(ra = 279.804, dec = 33.997, vinf = 26.3209)
    """
    
    ra *= np.pi/180
    dec *= np.pi/180

    cosdcosa = np.cos(dec) * np.cos(ra)
    cosdsina = np.cos(dec) * np.sin(ra)
    sind = np.sin(dec)
    
    U = -vinf * (T11*cosdcosa + T12*cosdsina + T13*sind)
    V = -vinf * (T21*cosdcosa + T22*cosdsina + T23*sind)
    W = -vinf * (T31*cosdcosa + T32*cosdsina + T33*sind)
    
    return U, V, W

