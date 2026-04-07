"""Prepare PDB structures for fpocket: remove HETATM, waters, keep ATOM + END."""
import os
import sys

def clean_pdb(input_path, output_path):
    """Strip non-protein atoms from PDB for pocket detection."""
    with open(input_path) as f:
        lines = f.readlines()
    
    kept = []
    for line in lines:
        if line.startswith('ATOM'):
            kept.append(line)
        elif line.startswith('END'):
            kept.append(line)
            break
    
    # For multi-model files, keep only first model
    cleaned = []
    for line in kept:
        if line.startswith('ENDMDL'):
            cleaned.append(line)
            break
        cleaned.append(line)
    
    if not any(l.startswith('END') for l in cleaned):
        cleaned.append('END\n')
    
    with open(output_path, 'w') as f:
        f.writelines(cleaned)
    
    n_atoms = sum(1 for l in cleaned if l.startswith('ATOM'))
    return n_atoms

if __name__ == '__main__':
    proj = sys.argv[1] if len(sys.argv) > 1 else '.'
    pdb_dir = os.path.join(proj, 'data', 'pdb')
    af_dir = os.path.join(proj, 'data', 'alphafold')
    out_dir = os.path.join(proj, 'data', 'cleaned')
    os.makedirs(out_dir, exist_ok=True)
    
    structures = {
        # Positive controls
        'ADAR2_crystal': os.path.join(pdb_dir, '1ZY7.pdb'),
        'Pds5B_crystal': os.path.join(pdb_dir, '5HDT.pdb'),
        'HDAC1_crystal': os.path.join(pdb_dir, '5ICN.pdb'),
        'HDAC3_crystal': os.path.join(pdb_dir, '4A69.pdb'),
        # AlphaFold
        'ADAR2_alphafold': os.path.join(af_dir, 'AF-P78563-F1.pdb'),
        # Negative controls (PH domains)
        'PLCd1_PH': os.path.join(pdb_dir, '1MAI.pdb'),
        'Btk_PH':   os.path.join(pdb_dir, '1BTK.pdb'),
        'DAPP1_PH': os.path.join(pdb_dir, '1FAO.pdb'),
        'Grp1_PH':  os.path.join(pdb_dir, '1FGY.pdb'),
    }
    
    for name, path in structures.items():
        if os.path.exists(path):
            out = os.path.join(out_dir, f'{name}.pdb')
            n = clean_pdb(path, out)
            print(f'  {name}: {n} atoms')
        else:
            print(f'  {name}: FILE NOT FOUND — {path}')
