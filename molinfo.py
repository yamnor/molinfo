import streamlit as st

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import SimilarityMaps, rdMolDraw2D

import py3Dmol
from stmol import showmol

import urllib.parse
from urllib.request import urlopen

def smi2mol(smi):
  mol = Chem.MolFromSmiles(smi)
  if mol is not None:
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, maxIters = 200)
    return mol
  else:
    return None

def show_2dview(smi):
  mol = Chem.MolFromSmiles(smi)
  if mol is not None:
    col = st.columns(3)
    col[0].write(' ')
    col[1].image(Draw.MolToImage(mol))
    col[2].write(' ')
  else:
    st.error('Try again.')

def show_3dview(smi):
  #viewsize = (300, 800)
  mol = smi2mol(smi)
  if mol is not None:
    viewer = py3Dmol.view(height = viewsize[0], width = viewsize[1])
    molblock = Chem.MolToMolBlock(mol)
    viewer.addModel(molblock, 'mol')
    viewer.setStyle({'stick':{}})
    viewer.zoomTo()
    viewer.spin('y', 1)
    #showmol(viewer, height = viewsize[0], width = viewsize[1])
    showmol(viewer)
    st.balloons()
  else:
    st.error('Try again.')

def get_alogps(smi):
  url = 'http://www.vcclab.org/web/alogps/calc?' + urllib.parse.urlencode({'SMILES': smi})
  txt = urlopen(url).read()
  dat = txt[44:-20].split()
  return {'logP' : float(dat[0]), 'logS' : float(dat[1])}

def show_properties(smi):
  mol = Chem.MolFromSmiles(smi)  
  alogps = get_alogps(smi)
  if mol is not None:
    col = st.columns(6)
    col[0].metric(label = "log P",      value = '{:.2f}'.format(alogps['logP']))
    col[1].metric(label = "log S",      value = '{:.2f}'.format(alogps['logS']))
    col[2].metric(label = "tPSA",       value = '{:.1f}'.format(Descriptors.TPSA(mol)))
    col[3].metric(label = "H-bond Accep", value = Descriptors.NumHAcceptors(mol))
    col[4].metric(label = "H-bond Donor", value = Descriptors.NumHDonors(mol))
    col[5].metric(label = "MW",         value = '{:.1f}'.format(Descriptors.MolWt(mol)))
  else:
    st.error('Try again.')

def main():

  st.set_page_config(
    page_title = 'molinfo')

  st.title("molinfo")
  
  molecule_dic = {
    'Enter SMILES' : '',
    'acetylsalicylic acid (Aspirin)' : 'CC(=O)Oc1ccccc1C(=O)O',
    'methyl salicylate' : 'O=C(OC)c1ccccc1O',
    'zanamivir (Relenza)' : 'CC(=O)NC1C(C=C(OC1C(C(CO)O)O)C(=O)O)N=C(N)N',
    'oseltamivir (Tamiflu)': 'CCC(CC)OC1C=C(CC(C1NC(=O)C)N)C(=O)OCC'}

  select_molecule = st.selectbox(
    "Select molecule",
    tuple(molecule_dic.keys()))

  if select_molecule == "Enter SMILES":

    smi = st.text_input('SMILES:', 'CC(=O)Oc1ccccc1C(=O)O')
    st.markdown(
      """
      [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
      is a specification for describing the chemical structure of molecules using short strings.
      """)

  else:

    smi = molecule_dic[select_molecule]

  if st.button("😊"):

    st.markdown("---")
    show_properties(smi)
    st.markdown("---")
    show_2dview(smi)
    st.markdown("---")
    show_3dview(smi)
    st.markdown("---")

if __name__ == "__main__":
    main()
