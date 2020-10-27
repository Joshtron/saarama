def prepare_top_dihedrals(top):
    from Bio import PDB
    import math

    structure = PDB.PDBParser().get_structure('input_structure', top)

    phi_gen = []
    psi_gen = []
    phi_pre = []
    psi_pre = []
    phi_pro = []
    psi_pro = []
    phi_gly = []
    psi_gly = []

    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            for poly_index, poly in enumerate(polypeptides):
                phi_psi = poly.get_phi_psi_list()
                for res_index, residue in enumerate(poly):
                    res_name = "{}".format(residue.resname)
                    phi, psi = phi_psi[res_index]
                    if phi and psi:
                        if str(poly[res_index + 1].resname) == "PRO":
                            phi_pre.append(math.degrees(phi))
                            psi_pre.append(math.degrees(psi))
                        elif res_name == "PRO":
                            phi_pro.append(math.degrees(phi))
                            psi_pro.append(math.degrees(psi))
                        elif res_name == "GLY":
                            phi_gly.append(math.degrees(phi))
                            psi_gly.append(math.degrees(psi))
                        else:
                            phi_gen.append(math.degrees(phi))
                            psi_gen.append(math.degrees(psi))

    return phi_gen, psi_gen, phi_pro, psi_pro, phi_gly, psi_gly, phi_pre, psi_pre










    #Previous attempt but for some reason MDanalysis calculates torsion angles wrong when the .pdb file
    #is too long or maybe if a chain is discontinous. I have not seen a pattern in this.
    #I switched to the dihedral method that was used in PYRAMA which seems to be a better solution.
    #I still wanted to leave this here in case I can fix this, which means that there's one import less

    '''
    from MDAnalysis.analysis.dihedrals import Ramachandran
    r_general = u.select_atoms("backbone and segid B and resname VAL PHE ALA LYS ARG CYS GLU LEU MET HIS TYR TRP SER ASN GLN THR ASP ILE")
    r_pro = u.select_atoms("resname PRO")
    r_gly = u.select_atoms("resname GLY")
    R_general = Ramachandran(r_general).run()
    R_pro = Ramachandran(r_pro).run()
    R_gly = Ramachandran(r_gly).run()

    for atom in u.select_atoms("backbone"):
        print(atom)

    phi_general = []
    psi_general = []

    for line in list(R_general.angles):
        for entry in line:
            splitted_entry = entry.tolist()
            phi_general.append(splitted_entry[0])
            psi_general.append(splitted_entry[1])

    phi_pro = []
    psi_pro = []

    for line in list(R_pro.angles):
        for entry in line:
            splitted_entry = entry.tolist()
            phi_pro.append(splitted_entry[0])
            psi_pro.append(splitted_entry[1])

    phi_gly = []
    psi_gly = []

    for line in list(R_gly.angles):
        for entry in line:
            splitted_entry = entry.tolist()
            phi_gly.append(splitted_entry[0])
            psi_gly.append(splitted_entry[1])
    '''

    '''

    res_counter = 2
    pro_list = []
    gly_list = []
    pre_list = []

    with open(top) as pdb:
        for atom in pdb:
            current_atom = atom.split(' ')
            filtered_current_atom = list(filter(lambda x: x != "", current_atom))
            if 'ATOM' in filtered_current_atom[0] and int(filtered_current_atom[int(id)]) == res_counter and res_counter < len(phi)+2:
                print(filtered_current_atom)
                if filtered_current_atom[3] == 'GLY':
                    gly_list.append(res_counter)
                elif filtered_current_atom[3] == 'PRO':
                    pro_list.append(res_counter)
                elif filtered_current_atom[3] == 'PRO' and res_counter > 2:
                    pre_list.append(res_counter - 1)
                res_counter += 1

    print(len(phi))
    '''



