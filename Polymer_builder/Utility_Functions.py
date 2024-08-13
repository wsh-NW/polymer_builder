import numpy as np
import matplotlib.pyplot as plt

########################################################################################################################
""" MATH FUNCTIONS"""

def unit_v(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def generate_random_unit_vector() -> np.array:
    """
    generates random unit vector (r=1) with phi and theta
    Source: https://stackoverflow.com/questions/5408276/sampling-uniformly-distributed-random-points-inside-a-spherical-volume
    :return: random unit vector
    """
    phi = np.random.uniform(0, 2 * np.pi)  # phi = random(0,2pi)
    costheta = np.random.uniform(-1, 1)  # costheta = random(-1,1)
    theta = np.arccos(costheta)
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])


def unit_basis_3pts(left_center_right: list[np.ndarray]) -> list[np.ndarray]:
    """
    creates unit vectors from an isoceles triangle:
    [k: points upwards out of triangle, w: orthogonal to C-C-C, v: crossproduct(k, w)]
    ** assumes isosceles triangle
    :param list[np.ndarray] left_center_right: 3 carbons, input from add_H function
    :return: list[np.ndarray] [upwards unit vector on center carbon, vector orthogonal to C-C-C plane, vector in ccc]
    """
    c1, c2, c3 = left_center_right
    c2c1 = c2 - c1
    c2c3 = c2 - c3
    k = unit_v(c2c1 + c2c3)  # axis to rotate H's off of
    w = unit_v(np.cross(c2c1, c2c3))  # axis orthogonal to c-c-c plane
    v = np.cross(w, k)  # axis completing rank 3 matrix
    return [k, w, v]


def rodrigues_rotation(angle: float, vector: np.ndarray, axis: np.ndarray) -> np.ndarray:
    """
    rotates about vector about axis <angle> degrees CCW
    :param angle: given in degrees
    :param vector: vector to rotate
    :param axis: vertical axis to rotate about
    :return: rotated vector
    Source: https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    Information: https://www.youtube.com/watch?v=q-ESzg03mQc&list=PLpzmRsG7u_gr0FO12cBWj-15_e0yqQQ1U&index=5
    """
    return np.dot(rotation_matrix(axis, angle), vector)


def rotation_matrix(axis: np.ndarray, theta: float):
    """
    Return the rotation matrix associated with counterclockwise rotation about the given axis by theta degrees.
    """
    axis = unit_v(axis)
    theta = theta * np.pi / 180  # converting to radians
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def convert_xyz(list_of_points: list[np.ndarray]):
    array3d = np.array(list_of_points)
    x, y, z = array3d[:, 0], array3d[:, 1], array3d[:, 2]
    return x, y, z

def angle_2vecs(vecs: list[np.ndarray]):
    return np.arccos(np.clip(np.dot(unit_v(vecs[0]), unit_v(vecs[1])), -1.0, 1.0)) * 180/np.pi


########################################################################################################################
""" PLOTTING FUNCTIONS """

def plot_carbon_positions(carbon_list):
    """
    returns n list of carbon positions: array[x,y,z]
    :return: list of 3D np.arrays
    """
    array3d = np.array(carbon_list)
    x_label, y_label, z_label = "x (A)", "y (A)", "z (A)"
    x, y, z = array3d[:, 0], array3d[:, 1], array3d[:, 2]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, marker='o')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    plt.title("carbon backbone")
    plt.show()


def plot_Polymer_CH2(carbon_list, hydrogen_list):
    """
    plots carbons and hydrogens
    :return: none
    """
    if not hydrogen_list:
        print("carbon backbone not complete, hydrogens not added")

    # carbon backbone plotted
    x_label, y_label, z_label = "x (A)", "y (A)", "z (A)"
    x, y, z = convert_xyz(carbon_list)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, marker='o')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    plt.title("CH2 Backbone")

    # hydrogens plotting (units of CH2 and CH3 connected)
    h_index = 3  # start on
    for i, c in enumerate(carbon_list):
        if i == 0:
            for h in hydrogen_list[0:3]:
                x, y, z = convert_xyz([h] + [c])
                plt.plot(x, y, z, marker="o")
        elif i == len(carbon_list) - 1:
            for h in hydrogen_list[-3:]:
                x, y, z = convert_xyz([h] + [c])
                plt.plot(x, y, z, marker="o")
        else:
            x, y, z = convert_xyz([hydrogen_list[h_index], c, hydrogen_list[h_index + 1]])
            plt.plot(x, y, z, marker="o")
            h_index += 2

    plt.show()

def plot_multiple_polymers(cn, hn, hydrogen_plotting=True):
    """
    plots carbons and hydrogens
    :return: none
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x_label, y_label, z_label = "x (A)", "y (A)", "z (A)"
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_zlabel(z_label)
    plt.title("Polymer Melt")

    #iterate over polymers
    for n in range(len(cn)):
        carbon_list = cn[n]
        hydrogen_list = hn[n]
        #plotting carbon backbone
        x, y, z = convert_xyz(carbon_list)
        ax.plot(x, y, z, marker='o')

        if hydrogen_plotting:
            # hydrogens plotting (units of CH2 and CH3 connected)
            h_index = 3  # start on
            for i, c in enumerate(carbon_list):
                if i == 0:
                    for h in hydrogen_list[0:3]:
                        x, y, z = convert_xyz([h] + [c])
                        plt.plot(x, y, z, marker="o", color="black", alpha=0.05)
                elif i == len(carbon_list) - 1:
                    for h in hydrogen_list[-3:]:
                        x, y, z = convert_xyz([h] + [c])
                        plt.plot(x, y, z, marker="o", color="black", alpha=0.05)
                else:
                    x, y, z = convert_xyz([hydrogen_list[h_index], c, hydrogen_list[h_index + 1]])
                    plt.plot(x, y, z, marker="o", color="black", alpha=0.05)
                    h_index += 2

    plt.show()