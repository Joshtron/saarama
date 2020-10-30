def prepare_all_dihedrals(path, u, res, id, start, end):

    #This is the core of torsion angle calculation where everything comes together-
    #atoms is a list that contains all lines of the original topology file.

    atoms = []

    with open(path) as file:
        for line in list(file):
            if not len(line.split()) == 0:
                if line.split()[0] == 'ATOM':
                    atoms.append(line)

    #This function gets the positions of needed atoms and their coordinates from the trajectory file
    #and returns them as well as the residue names.

    phi, psi, res_name_list = get_psi_phi(atoms, u, res, id, start, end)

    phi_list = []
    psi_list = []
    j = 0

    #This actually calculates the torsion angles with help of the coordinates.

    for i in range(0, (end-start)):
        angle_list = [phi[i], psi[i]]
        phi_temp = []
        psi_temp = []
        for angle in angle_list:
            for sublist in angle:
                p1 = sublist[0]
                p2 = sublist[1]
                p3 = sublist[2]
                p4 = sublist[3]
                final_angle = get_dihedrals(p1, p2, p3, p4)
                if j == 0 or (j % 2) == 0:
                    phi_temp.append(final_angle)
                else:
                    psi_temp.append(final_angle)
            j += 1




        phi_list.append(phi_temp)
        psi_list.append(psi_temp)

    return phi_list, psi_list, res_name_list

def get_psi_phi(atoms, u, res, id, start, end):

    atom_list = []
    entry_number = []
    pdb_list = []
    res_name_list = []

    #This gives us the residue names that are used for plotting later

    for atom in atoms:
        current_atom = atom.split(' ')
        fca = list(filter(lambda x: x != "", current_atom))
        pdb_list.append(fca)

    i = start

    for atom in pdb_list:
        if int(atom[4]) == i and i <= end:
            res_name_list.append(atom[3])
            i += 1

    #This atrocity is actually quite useful and helps you by looping through all atoms
    #that are going to be plotted. Basically this extracts all positions that are needed in order
    #to calculate torsion angles. Those positions of important atoms gets saved in atom_list and entry_number
    #so they can get looked up in the trajectory file to extract the coordinates.

    for i in range(start, end):
        temp_atom = []
        temp_entry = []
        j = 0
        for atom in atoms:
            current_atom = atom.split(' ')
            fca = list(filter(lambda x: x != "", current_atom))
            if 'ATOM' == fca[0]:
                if len(fca) == 12 and int(fca[id]) == int(i)-1:
                    if fca[2] == 'C':
                        temp_atom.append(list(map(float, fca[5:8])))
                        temp_entry.append(j)
                if len(fca) == 12 and int(fca[id]) == int(i):
                    if fca[2] == 'N' or fca[2] == 'CA' or fca[2] == 'C':
                        temp_atom.append(list(map(float, fca[5:8])))
                        temp_entry.append(j)
                if len(fca) == 12 and int(fca[id]) == int(i)+1:
                    if fca[2] == 'N':
                        temp_atom.append(list(map(float, fca[5:8])))
                        temp_entry.append(j)
                j += 1

        atom_list.append(temp_atom)
        entry_number.append(temp_entry)

    dihedral_list = []

    #Coordinates get extracted here

    for i in range(0, (end-start)):
        temp_dihedrals = []
        for ts in u.trajectory:
            pos = ts.positions
            temp_list = []
            for index in entry_number[i]:
                temp_list.append(pos[index].tolist())
            temp_dihedrals.append(temp_list)
        dihedral_list.append(temp_dihedrals)

    phi = []
    psi = []

    #We get the coordinates from 5 atoms where atom 1-4 forms the phi angles and atom 2-5 froms the
    #psi angle. coordinates get splitted here.

    for i in range(0, (end-start)):
        phi_temp = []
        psi_temp = []
        for sub in dihedral_list[i]:
            phi_sub = [sub[0], sub[1], sub[2], sub[3]]
            psi_sub = [sub[1], sub[2], sub[3], sub[4]]
            phi_temp.append(phi_sub)
            psi_temp.append(psi_sub)
        phi.append(phi_temp)
        psi.append(psi_temp)

    return phi, psi, res_name_list

#Not really elegant but this calculates torsion angles the old fashioned way and it is somehow guranteed that it works
#Trust me, I tested this a lot and validated most of my angles with vmd.

def get_dihedrals(a, b, c, d):
    q1, q2, q3 = calc_q_vectors(a, b, c, d)
    q1_x_q2, q2_x_q3 = calc_cross_vectors(q1, q2, q3)
    n1, n2 = calc_normals(q1_x_q2, q2_x_q3)
    u1, u2, u3 = calc_orthogonal_unit_vectors(n2, q2)
    final_angle = calc_dihedral_angle(n1, u1, u2, u3)
    return final_angle

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