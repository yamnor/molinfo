import streamlit as st

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors

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
  viewsize = (340, 340)
  mol = smi2mol(smi)  
  if mol is not None:
    col = st.columns(3)
    col[0].write(' ')
    with col[1]:
      viewer = py3Dmol.view(height = viewsize[0], width = viewsize[1])
      molblock = Chem.MolToMolBlock(mol)
      viewer.addModel(molblock, 'mol')
      viewer.setStyle({'stick':{}})
      viewer.zoomTo()
      viewer.spin('y', 1)
      viewer.setBackgroundColor('black')
      showmol(viewer, height = viewsize[0], width = viewsize[1])
    col[2].write(' ')
  else:
    st.error('Try again.')

def get_alogps(smi):
  url = 'http://www.vcclab.org/web/alogps/calc?' + urllib.parse.urlencode({'SMILES': smi})
  txt = urlopen(url).read()
  if 'logP' in str(txt):
    dat = txt[44:-20].split()
    return {'logP' : float(dat[0]), 'logS' : float(dat[1])}
  else:
    return {'logP' : 0.0, 'logS' : 0.0}

def show_properties(smi):
  mol = Chem.MolFromSmiles(smi)  
  alogps = get_alogps(smi)
  if mol is not None:
    col = st.columns(6)
    col[0].metric(label = "logP", value = '{:.2f}'.format(alogps['logP']))
    col[1].metric(label = "logS", value = '{:.2f}'.format(alogps['logS']))
    col[2].metric(label = "tPSA", value = '{:.1f}'.format(Descriptors.TPSA(mol)))
    col[3].metric(label = "HBAc",  value = Descriptors.NumHAcceptors(mol))
    col[4].metric(label = "HBDo",  value = Descriptors.NumHDonors(mol))
    col[5].metric(label = "RotB",  value = Descriptors.NumRotatableBonds(mol))
    '''
    `logP` : [オクタノール・水分配係数](https://w.wiki/62VE)
    `logS` : 水に対する溶解度
    `tPSA` : [トポロジカル極性表面積](https://w.wiki/62VA)（**t**opological **P**olar **S**urface **A**rea）
    `HBAc` : 水素結合アクセプター（**H**ydrogen **B**ond **Ac**ceptors）の数
    `HBDo` : 水素結合ドナー（**H**ydrogen **B**ond **Do**nors）の数
    `RotB` : 回転できる共有結合（**Rot**atable **B**onds）の数
    '''
  else:
    st.error('Try again.')

def show_info(smi):
  st.balloons()
  show_2dview(smi)
  st.markdown("---")
  with st.spinner("Please wait..."):
    show_properties(smi)
  st.markdown("---")
  show_3dview(smi)

def main():

  st.set_page_config(
    page_title = 'molinfo',
    initial_sidebar_state = 'collapsed')

  st.title("molinfo")
  
  molecule_dic = {
    'Enter SMILES' : '',
    'acetylsalicylic acid (Aspirin)' : 'CC(=O)Oc1ccccc1C(=O)O',
    'methyl salicylate' : 'O=C(OC)c1ccccc1O',
    'zanamivir (Relenza)' : 'CC(=O)NC1C(C=C(OC1C(C(CO)O)O)C(=O)O)N=C(N)N',
    'oseltamivir (Tamiflu)': 'CCC(CC)OC1C=C(CC(C1NC(=O)C)N)C(=O)OCC',
    'benzene hexachloride (BHC)' : 'C1(C(C(C(C(C1Cl)Cl)Cl)Cl)Cl)Cl',
    'pentachlorobiphenyl (PCB)' : 'C1=CC(=C(C=C1C2=CC(=C(C(=C2)Cl)Cl)Cl)Cl)Cl',
    }

  select_molecule = st.selectbox(
    "Select molecule",
    tuple(molecule_dic.keys()))

  if select_molecule == "Enter SMILES":
    smi = st.text_input('SMILES:', 'CC(=O)Oc1ccccc1C(=O)O')
    if st.button("😊"):
      show_info(smi)
  else:
    smi = molecule_dic[select_molecule]
    show_info(smi)

  with st.sidebar:

    st.subheader('Notice')
    st.markdown(
      """
      * This [streamlit](https://streamlit.io/) app uses [RDKit](https://www.rdkit.org),
        [py3Dmol](https://pypi.org/project/py3Dmol/), and [stmol](https://github.com/napoles-uach/stmol) libraries.
      * `log P` and `log S` are obtained from [the ALOGPS homepage](http://www.vcclab.org/lab/alogps/).
      """
    )

if __name__ == "__main__":
    main()
