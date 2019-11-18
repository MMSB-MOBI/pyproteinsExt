import ccmap as c_ccmap

def ccmap(pdbRec, pdbLig, dist):
    int_map = c_ccmap.duals([(pdbRec.atomDictorize, pdbLig.atomDictorize)], dist)
    lig_residues = pdbLig.getResID
    rec_residues = pdbRec.getResID
    ligResCount = len(lig_residues)
    indexes=[(int(i/ligResCount), i% ligResCount) for i in int_map[0]]

    tr_ccmap = [(rec_residues[i[0]], lig_residues[i[1]]) for i in indexes]
    return tr_ccmap