def calc_q_vectors(p1, p2, p3, p4):
    """Function to calculate q vectors"""
    import numpy as np
    # Calculate coordinates for vectors q1, q2 and q3
    q1 = np.subtract(p2,p1) # b-a
    q2 = np.subtract(p3,p2) # c-b
    q3 = np.subtract(p4,p3) # d-c
    return q1,q2,q3

def calc_cross_vectors(q1,q2,q3):
    """Function to calculate cross vectors"""
    import numpy as np
    # Calculate cross vectors
    q1_x_q2 = np.cross(q1,q2)
    q2_x_q3 = np.cross(q2,q3)
    return q1_x_q2, q2_x_q3

def calc_normals(q1_x_q2,q2_x_q3):
    """Function to calculate normal vectors to planes"""
    import numpy as np
    # Calculate normal vectors
    n1 = q1_x_q2/np.sqrt(np.dot(q1_x_q2,q1_x_q2))
    n2 = q2_x_q3/np.sqrt(np.dot(q2_x_q3,q2_x_q3))
    return n1,n2

def calc_orthogonal_unit_vectors(n2,q2):
    """Function to calculate orthogonal unit vectors"""
    import numpy as np
    # Calculate unit vectors
    u1 = n2
    u3 = q2/(np.sqrt(np.dot(q2,q2)))
    u2 = np.cross(u3,u1)
    return u1,u2,u3

def calc_dihedral_angle(n1,u1,u2,u3):
    """Function to calculate dihedral angle"""
    import numpy as np
    import math
    # Calculate cosine and sine
    cos_theta = np.dot(n1,u1)
    sin_theta = np.dot(n1,u2)
    # Calculate theta
    theta = -math.atan2(sin_theta,cos_theta)    # it is different from Fortran math.atan2(y,x)
    theta_deg = np.degrees(theta)
    # Show results
    #print("theta (rad) = %8.3f"%theta)
    final_angle = float("%8.3f"%theta_deg)
    return final_angle

def get_dihedrals(a, b, c, d):
    q1, q2, q3 = calc_q_vectors(a, b, c, d)
    q1_x_q2, q2_x_q3 = calc_cross_vectors(q1, q2, q3)
    n1, n2 = calc_normals(q1_x_q2, q2_x_q3)
    u1, u2, u3 = calc_orthogonal_unit_vectors(n2, q2)
    final_angle = calc_dihedral_angle(n1, u1, u2, u3)
    return final_angle

def get_psi_phi(atoms, u, res, id):

    atom_list = []
    entry_number = []
    i = 0

    for atom in atoms:
        current_atom = atom.split(' ')
        fca = list(filter(lambda x: x != "", current_atom))
        if 'ATOM' == fca[0]:
            if len(fca) == 12 and int(fca[id]) == int(res)-1:
                if fca[2] == 'C':
                    atom_list.append(list(map(float, fca[5:8])))
                    entry_number.append(i)
            if len(fca) == 12 and int(fca[id]) == int(res):
                if fca[2] == 'N' or fca[2] == 'CA' or fca[2] == 'C':
                    atom_list.append(list(map(float, fca[5:8])))
                    entry_number.append(i)
            if len(fca) == 12 and int(fca[id]) == int(res)+1:
                if fca[2] == 'N':
                    atom_list.append(list(map(float, fca[5:8])))
                    entry_number.append(i)
            i += 1

    dihedral_list = []



    for ts in u.trajectory:
        pos = ts.positions
        temp_list = []
        for index in entry_number:
            temp_list.append(pos[index].tolist())
        dihedral_list.append(temp_list)

    phi = []
    psi = []

    for sub in dihedral_list:
        phi_temp = [sub[0], sub[1], sub[2], sub[3]]
        psi_temp = [sub[1], sub[2], sub[3], sub[4]]
        phi.append(phi_temp)
        psi.append(psi_temp)

    return phi, psi


def prepare_capped_dihedrals(path, u, res, id):

    atoms = []

    with open(path) as file:
        for line in list(file):
            if not len(line.split()) == 0:
                if line.split()[0] == 'ATOM':
                    atoms.append(line)

    phi, psi = get_psi_phi(atoms, u, res, id)

    phi_list = []
    psi_list = []
    i = 0
    angle_list = [phi ,psi]

    for angle in angle_list:
        for sublist in angle:
            p1 = sublist[0]
            p2 = sublist[1]
            p3 = sublist[2]
            p4 = sublist[3]
            finale_angle = get_dihedrals(p1, p2, p3, p4)
            if i == 0:
                phi_list.append(finale_angle)
            else:
                psi_list.append(finale_angle)
        i += 1

    return phi_list, psi_list


