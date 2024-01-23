import numpy as np

# dot product of many vectors in a single array, where the first dimension is the component
def many_mags(a):
    return np.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])

def many_dots(a,b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def many_crosses(a,b):
    xs = a[1]*b[2] - a[2]*b[1]
    ys = a[2]*b[0] - a[0]*b[2]
    zs = a[0]*b[1] - a[1]*b[0]
    return np.stack((xs,ys,zs),axis=0)

def many_ijk_to_xyz(s_ijk, i_xyz, j_xyz, k_xyz):
    sii = np.stack((s_ijk[0]*i_xyz[0],s_ijk[0]*i_xyz[1],s_ijk[0]*i_xyz[2]),axis=0)
    sjj = np.stack((s_ijk[1]*j_xyz[0],s_ijk[1]*j_xyz[1],s_ijk[1]*j_xyz[2]),axis=0)
    skk = np.stack((s_ijk[2]*k_xyz[0],s_ijk[2]*k_xyz[1],s_ijk[2]*k_xyz[2]),axis=0)
    s_xyz = sii+sjj+skk
    s_xyz /= many_mags(s_xyz)
    return s_xyz

# returns orbit normal unit vector for orbit of body at index index
def calc_orbit_normal(r_xyz=None,v_xyz=None,incl=np.radians(20.),outer=False):
    if outer:
        l = np.array([0.,-np.sin(incl),np.cos(incl)])
        l_hat = l/np.sqrt(np.dot(l,l))
        l_hats = np.array([[l_hat[j] for i in range(np.shape(r_xyz)[1])] for j in range(3)])
    else:
        l = many_crosses(v_xyz,r_xyz)
        l_hats = l/many_mags(l)
    return l_hats

# get theta, assuming l_hat is in z direction
def get_theta(s_ijk,i_xyz,j_xyz,k_xyz):
    s_x = s_ijk[0]*i_xyz[0] + s_ijk[1]*j_xyz[0] + s_ijk[2]*k_xyz[0]
    s_y = s_ijk[0]*i_xyz[1] + s_ijk[1]*j_xyz[1] + s_ijk[2]*k_xyz[1]
    s_z = s_ijk[0]*i_xyz[2] + s_ijk[1]*j_xyz[2] + s_ijk[2]*k_xyz[2]
    s_xyz = np.stack((s_x,s_y,s_z),axis=0)
    s_xyz /= many_mags(s_xyz)
    return np.degrees(np.arccos(s_xyz[2]))

# returns (obliquity, phi) of body at index 1 in radians
def get_theta_phi(s_ijk,i_xyz,j_xyz,k_xyz,r,v):
    # calculate theta
    s_x = s_ijk[0]*i_xyz[0] + s_ijk[1]*j_xyz[0] + s_ijk[2]*k_xyz[0]
    s_y = s_ijk[0]*i_xyz[1] + s_ijk[1]*j_xyz[1] + s_ijk[2]*k_xyz[1]
    s_z = s_ijk[0]*i_xyz[2] + s_ijk[1]*j_xyz[2] + s_ijk[2]*k_xyz[2]
    s_xyz = np.stack((s_x,s_y,s_z),axis=0)
    s_xyz /= many_mags(s_xyz)

    l_hat = calc_orbit_normal(r_xyz=r,v_xyz=v) # orbit normal of triaxial planet
    theta = np.arccos(many_dots(l_hat,s_xyz))

    # calculate phi
    l_p_hat = calc_orbit_normal(r_xyz=r,outer=True) # orbit normal of perturbing planet
    y = many_crosses(l_p_hat, l_hat) # unrelated to y basis unit vector
    y_hat = y/many_mags(y)
    x = many_crosses(y_hat, l_hat) # unrelated to x basis unit vector
    x_hat = x/many_mags(x)
    phi = np.arctan2(many_dots(s_xyz,y_hat),many_dots(s_xyz,x_hat))

    # range from 0 to 360
    phi += (phi < 0).astype(int) * 2*np.pi

    return (theta, phi)

# returns angle in rad
def get_psi(rs, iss, js):
    psi = np.arccos(many_dots(rs,iss))
    # range from 0 to 360
    # psi *= ((2 * (many_dots(rs,js) > 0).astype(int)) - 1) # add sign to [0,pi] angle
    # psi += (psi < 0).astype(int) * 2*np.pi # range from 0 to 180
    return psi

# psi is angle between r and {j cross s, or k cross s if s = j}
# psi is positive if in direction of spin, negative if otherwise
def get_psi_v2(rs,i_xyz,j_xyz,k_xyz,s_ijk):
    s_x = s_ijk[0]*i_xyz[0] + s_ijk[1]*j_xyz[0] + s_ijk[2]*k_xyz[0]
    s_y = s_ijk[0]*i_xyz[1] + s_ijk[1]*j_xyz[1] + s_ijk[2]*k_xyz[1]
    s_z = s_ijk[0]*i_xyz[2] + s_ijk[1]*j_xyz[2] + s_ijk[2]*k_xyz[2]
    s_xyz = np.stack((s_x,s_y,s_z),axis=0)
    s_xyz /= many_mags(s_xyz)

    long_axes = many_crosses(j_xyz,s_xyz)
    long_axes += (long_axes == 0.).astype(int) * many_crosses(k_xyz,s_xyz)
    long_axes /= many_mags(long_axes)
    signs = ((many_dots(many_crosses(rs,long_axes), s_xyz) > 0).astype(int)*2) - 1
    proto_psis = np.cos(many_dots(long_axes,rs))
    psis = signs * proto_psis
    return psis

# returns angle in rad
def get_beta(ss):
    proto_beta = np.arccos(ss[2])
    beta = (proto_beta*(proto_beta <= np.pi/2).astype(int)) + ((np.pi - proto_beta)*(proto_beta > np.pi/2).astype(int))
    return beta

def get_theta_kl(ks,r,v):
    l_hat = calc_orbit_normal(r_xyz=r,v_xyz=v)
    theta_kl = np.arccos(many_dots(ks,l_hat))
    return theta_kl

def get_theta_primes(triaxial_bool,ss,r,v,i_xyz,j_xyz,k_xyz):
    Re = 4.263e-5 # radius of Earth in AU
    Me = 3.003e-6 # mass of Earth in solar masses
    R_p = 2.*Re # radius of inner planet
    M_p = 4.*Me # mass of inner planet

    if triaxial_bool:
        moment2 = 1.e-5 # 1e-1 # (Ij - Ii) / Ii, < moment3
    else:
        moment2 = 0.
    moment3 = 1.e-3 # 2e-1 # (Ik - Ii) / Ii, > moment2
    k = 0.331
    Ii = k*M_p*R_p*R_p
    Ij = Ii*(1+moment2)
    Ik = Ii*(1+moment3)

    vec = np.stack((Ii*ss[0],Ij*ss[1],Ik*ss[2]),axis=0)
    vec_xyz = np.stack((vec[0]*i_xyz[0] + vec[1]*j_xyz[0] + vec[2]*k_xyz[0],vec[0]*i_xyz[1] + vec[1]*j_xyz[1] + vec[2]*k_xyz[1],vec[0]*i_xyz[2] + vec[1]*j_xyz[2] + vec[2]*k_xyz[2]),axis=0)
    vec_xyz /= many_mags(vec_xyz)
    return many_dots(vec_xyz,calc_orbit_normal(r_xyz=r,v_xyz=v))
