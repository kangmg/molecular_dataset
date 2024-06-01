from rdkit import Chem
from rdkit.Chem import AllChem
from chemspipy import ChemSpider


def smiles2mol(smiles:str, N_conformers:int)->str:
  """
    get mol format block from smiles
  """
  # smiles to mol block
  mol = Chem.MolFromSmiles(smiles)
  mol = Chem.AddHs(mol)

  # sanitize
  Chem.SanitizeMol(mol)

  # generate 3D embedded conformers
  num_conformers = N_conformers
  params = AllChem.ETKDGv3() # embedding algorithm
  params.randomSeed = 42

  conformer_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers)

  energies = []
  for conf_id in conformer_ids:
    # MMFF94 힘장 설정
    ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)

    # force field optimize
    ff.Minimize()
    energy = ff.CalcEnergy()
    energies.append((conf_id, energy))

  # get lowest conformer
  lowest_energy_idx, _ = min(energies, key=lambda x: x[-1])

  # mol format block
  mol_block = Chem.MolToMolBlock(mol, confId=lowest_energy_idx)

  return mol_block


def get_properties(_API_KEY_:str, CSID:int|str)->dict:
    """
    Description
    -----------
    CSID (ChemSpider ID)로부터 name, smiles, formula를 얻은 후 반환하는 함수

    Parameters
    ----------
      - _API_KEY_ : chemspider API 키
      - CSID : chemspider ID 값
    
    Returns
    -------
      - properties (dict) : name, smiles, formula를 포함한 dict
    """
    # get properties using chemspider API
    cs = ChemSpider(_API_KEY_)
    compound = cs.get_compound(CSID)
    # parsing
    try:
      NAME = compound.common_name
    except KeyError:
      NAME = ''
    SMILES = compound.smiles
    FORMULA = compound.molecular_formula
    return {"CSID"      :   CSID,
            "NAME"      :   NAME,
            "SMILES"    :   SMILES,
            "FORMULA"   :   FORMULA}


def get_charge(smiles:str)->int:
  """
  Description
  -----------
  estimates net charge from smiles

  Parameters
  ----------
  smiles(str) : simplified molecular-input line-entry system

  Returns
  -------
  net charge (int) : estimated net charge
  """
  try:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.ComputeGasteigerCharges(mol)
    _charges = [atom.GetDoubleProp('_GasteigerCharge') for atom in mol.GetAtoms()]
    return round(sum(_charges))
  except:
    print(f"[Missing Value Warning] rdkit - ComputeGasteigerCharges may not supports to compute formal charge of {smiles}")
    return None