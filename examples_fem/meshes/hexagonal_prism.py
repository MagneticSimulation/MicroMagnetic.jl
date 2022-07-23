import numpy as np

#compute the rotation_axis, typically axis_from and axis_to have a acute angle
def compute_rotation_axis(axis_from, axis_to):
    axis = np.cross(axis_from, axis_to)
    if np.sqrt(sum(axis**2)) == 0.0:
        return axis_from
    return axis / np.sqrt(sum(axis**2))


#rotate the vector v around the rotation axis with theta degree.
def rotate_to(rotation_axis, v, theta):
    matrix = np.zeros((3,3))
    m = rotation_axis
    st, ct = np.sin(theta), np.cos(theta)
    matrix[0,0] = ct+(1-ct)*m[0]*m[0]
    matrix[0,1] = (1-ct)*m[0]*m[1] - st*m[2]
    matrix[0,2] = (1-ct)*m[0]*m[2] + st*m[1]
    matrix[1,0] = (1-ct)*m[1]*m[0] + st*m[2]
    matrix[1,1] = ct+(1-ct)*m[1]*m[1]
    matrix[1,2] = (1-ct)*m[1]*m[2] - st*m[0]
    matrix[2,0] = (1-ct)*m[0]*m[2] - st*m[1]
    matrix[2,1] = (1-ct)*m[1]*m[2] + st*m[0]
    matrix[2,2] = ct+(1-ct)*m[2]*m[2]
    return np.dot(matrix,v)

#rotate the given vector when the z axis is rotated to `axis_to`.
def rotate(vector, theta=0, phi=0, gamma=0):
    nx, ny, nz = np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta) 
    axis = compute_rotation_axis([0,0,1], [nx, ny, nz])
    theta = np.arccos(nz)
    vector = rotate_to(axis, vector, theta) #rotate_to (nx, ny, nz)
    vector = rotate_to([nx, ny, nz], vector, gamma) #rotate around (nx, ny, nz)
    return vector

# the center of the hexagonal prism is (x0,y0,z0)
# the height of the hexagonal prism is h
# the radius of the hexagonal prism is R
# the normal direction of the hexagonal prism is first rotated to (theta, phi), i.e.,
# n = (cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta))
# then the hexagonal prism is rotated gamma around its normal axis
def hexagonal_prism(x0,y0,z0, h=1, R=2, theta=0, phi=0, gamma=0, maxh=0.1, id=0):

    def single_surface_info(dp):
        dv = rotate(dp, theta, phi, gamma)
        x, y, z = x0 + dv[0], y0 + dv[1], z0 + dv[2]
        return "plane(%g, %g, %g; %g, %g, %g)\n"%(x, y, z, dv[0], dv[1], dv[2]) 

    name = "hp_%d"%id
    output = "solid %s = "%name
    output += single_surface_info([0, 0, h/2.0])  #top layer:
    output += "and " + single_surface_info([0, 0, -h/2.0]) #bottom layer:
    output += "and " + single_surface_info([0, -R*np.sqrt(3)/2, 0])  #front
    output += "and " + single_surface_info([0, R*np.sqrt(3)/2, 0]) #back
    output += "and " + single_surface_info([-R*3/4.0, -R*np.sqrt(3)/4, 0]) #left front
    output += "and " + single_surface_info([-R*3/4.0, R*np.sqrt(3)/4, 0]) #left back
    output += "and " + single_surface_info([R*3/4.0, -R*np.sqrt(3)/4, 0]) #right front
    output += "and " + single_surface_info([R*3/4.0, R*np.sqrt(3)/4, 0])[:-1]+";\n" #right back

    output += "tlo %s -maxh=%g;\n\n"%(name, maxh)

    nx, ny, nz = np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta), np.cos(theta) 
    print("normal: ", [nx, ny, nz])

    return output


def main(output="hex.geo"):
    f = open(output, "w")
    f.write("algebraic3d \n\n")

    R = 21/2.0
    h = 7.0
    for i in range(2):
        x = i*25
        theta = (np.pi/2)*(2*np.random.random()-1)
        phi = 2*np.pi*np.random.random()
        gamma = np.pi/3*np.random.random()
        p = hexagonal_prism(x, 0, 0, h=h, R=R, theta=theta, phi=phi, gamma=gamma, maxh=2.0, id=i)
        f.write(p)
    f.close()  

main()

