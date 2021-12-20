#this contains the fragmentation function - produces RDKit fragments and atoms
import rdkit
import pandas as pd
import numpy as np
from rdkit import Chem
import rdkit.Chem.Fragments
from rdkit.Chem import rdqueries    

df1 = pd.DataFrame()
df2 = pd.DataFrame()
   
def fragments(mol):
    Al_COO = []
    Al_OH = []
    Al_OH_noTert = []
    ArN = []
    Ar_COO = []
    Ar_N = []
    Ar_NH = []
    Ar_OH = []
    COO = []
    COO2 = []
    C_O = []
    C_O_noCOO = []
    C_S = []
    HOCCN = []
    Imine = []
    NH0 = []
    NH1 = []
    NH2 = []
    N_O = []
    Ndealkylation1 = []
    Ndealkylation2 = []
    Nhpyrrole = []
    SH = []
    aldehyde = []
    alkyl_carbamate = []
    alkyl_halide = []
    allylic_oxid = []
    amide = []
    amidine = []
    aniline = []
    aryl_methyl = []
    azide = []
    azo = []
    barbitur = []
    benzene = []
    benzodiazepine = []
    bicyclic = []
    diazo = []
    dihydropyridine = []
    epoxide = []
    ester = []
    ether = []
    furan = []
    guanido = []
    halogen = []
    hdrzine = []
    hdrzone = []
    imidazole = []
    imide = []
    isocyan = []
    isothiocyan = []
    ketone = []
    ketone_Topliss = []
    lactam = []
    lactone = []
    methoxy = []
    morpholine = []
    nitrile = []
    nitro = []
    nitro_arom = []
    nitro_arom_nonortho = []
    nitroso = []
    oxazole = []
    oxime = []
    para_hydroxylation = []
    phenol = []
    phenol_noOrthoHbond = []
    phos_acid = []
    phos_ester = []
    piperdine = []
    piperzine = []
    priamide = []
    prisulfonamd = []
    pyridine = []
    quatN = []
    sulfide = []
    sulfonamd = []
    sulfone = []
    term_acetylene = []
    tetrazole = []
    thiazole = []
    thiocyan = []
    thiophene = []
    unbrch_alkane = []
    urea = []

    Al_COO.append(rdkit.Chem.Fragments.fr_Al_COO(mol))
    Al_OH.append(rdkit.Chem.Fragments.fr_Al_OH(mol))
    Al_OH_noTert.append(rdkit.Chem.Fragments.fr_Al_OH_noTert(mol))
    ArN.append(rdkit.Chem.Fragments.fr_ArN(mol))
    Ar_COO.append(rdkit.Chem.Fragments.fr_Ar_COO(mol))
    Ar_N.append(rdkit.Chem.Fragments.fr_Ar_N(mol))
    Ar_NH.append(rdkit.Chem.Fragments.fr_Ar_NH(mol))
    Ar_OH.append(rdkit.Chem.Fragments.fr_Ar_OH(mol))
    COO.append(rdkit.Chem.Fragments.fr_COO(mol))
    COO2.append(rdkit.Chem.Fragments.fr_COO2(mol))
    C_O.append(rdkit.Chem.Fragments.fr_C_O(mol))
    C_O_noCOO.append(rdkit.Chem.Fragments.fr_C_O_noCOO(mol))
    C_S.append(rdkit.Chem.Fragments.fr_C_S(mol))
    HOCCN.append(rdkit.Chem.Fragments.fr_HOCCN(mol))
    Imine.append(rdkit.Chem.Fragments.fr_Imine(mol))
    NH0.append(rdkit.Chem.Fragments.fr_NH0(mol))
    NH1.append(rdkit.Chem.Fragments.fr_NH1(mol))
    NH2.append(rdkit.Chem.Fragments.fr_NH2(mol))
    N_O.append(rdkit.Chem.Fragments.fr_N_O(mol))
    Ndealkylation1.append(rdkit.Chem.Fragments.fr_Ndealkylation1(mol))
    Ndealkylation2.append(rdkit.Chem.Fragments.fr_Ndealkylation2(mol))
    Nhpyrrole.append(rdkit.Chem.Fragments.fr_Nhpyrrole(mol))
    SH.append(rdkit.Chem.Fragments.fr_SH(mol))
    aldehyde.append(rdkit.Chem.Fragments.fr_aldehyde(mol))
    alkyl_carbamate.append(rdkit.Chem.Fragments.fr_alkyl_carbamate(mol))
    alkyl_halide.append(rdkit.Chem.Fragments.fr_alkyl_halide(mol))
    allylic_oxid.append(rdkit.Chem.Fragments.fr_allylic_oxid(mol))
    amide.append(rdkit.Chem.Fragments.fr_amide(mol))
    amidine.append(rdkit.Chem.Fragments.fr_amidine(mol))
    aniline.append(rdkit.Chem.Fragments.fr_aniline(mol))
    aryl_methyl.append(rdkit.Chem.Fragments.fr_aryl_methyl(mol))
    azide.append(rdkit.Chem.Fragments.fr_azide(mol))
    azo.append(rdkit.Chem.Fragments.fr_azo(mol))
    barbitur.append(rdkit.Chem.Fragments.fr_barbitur(mol))
    benzene.append(rdkit.Chem.Fragments.fr_benzene(mol))
    benzodiazepine.append(rdkit.Chem.Fragments.fr_benzodiazepine(mol))
    bicyclic.append(rdkit.Chem.Fragments.fr_bicyclic(mol))
    diazo.append(rdkit.Chem.Fragments.fr_diazo(mol))
    dihydropyridine.append(rdkit.Chem.Fragments.fr_dihydropyridine(mol))
    epoxide.append(rdkit.Chem.Fragments.fr_epoxide(mol))
    ester.append(rdkit.Chem.Fragments.fr_ester(mol))
    ether.append(rdkit.Chem.Fragments.fr_ether(mol))
    furan.append(rdkit.Chem.Fragments.fr_furan(mol))
    guanido.append(rdkit.Chem.Fragments.fr_guanido(mol))
    halogen.append(rdkit.Chem.Fragments.fr_halogen(mol))
    hdrzine.append(rdkit.Chem.Fragments.fr_hdrzine(mol))
    hdrzone.append(rdkit.Chem.Fragments.fr_hdrzone(mol))
    imidazole.append(rdkit.Chem.Fragments.fr_imidazole(mol))
    imide.append(rdkit.Chem.Fragments.fr_imide(mol))
    isocyan.append(rdkit.Chem.Fragments.fr_isocyan(mol))
    isothiocyan.append(rdkit.Chem.Fragments.fr_isothiocyan(mol))
    ketone.append(rdkit.Chem.Fragments.fr_ketone(mol))
    ketone_Topliss.append(rdkit.Chem.Fragments.fr_ketone_Topliss(mol))
    lactam.append(rdkit.Chem.Fragments.fr_lactam(mol))
    lactone.append(rdkit.Chem.Fragments.fr_lactone(mol))
    methoxy.append(rdkit.Chem.Fragments.fr_methoxy(mol))
    morpholine.append(rdkit.Chem.Fragments.fr_morpholine(mol))
    nitrile.append(rdkit.Chem.Fragments.fr_nitrile(mol))
    nitro.append(rdkit.Chem.Fragments.fr_nitro(mol))
    nitro_arom.append(rdkit.Chem.Fragments.fr_nitro_arom(mol))
    nitro_arom_nonortho.append(rdkit.Chem.Fragments.fr_nitro_arom_nonortho(mol))
    nitroso.append(rdkit.Chem.Fragments.fr_nitroso(mol))
    oxazole.append(rdkit.Chem.Fragments.fr_oxazole(mol))
    oxime.append(rdkit.Chem.Fragments.fr_oxime(mol))
    para_hydroxylation.append(rdkit.Chem.Fragments.fr_para_hydroxylation(mol))
    phenol.append(rdkit.Chem.Fragments.fr_phenol(mol))
    phenol_noOrthoHbond.append(rdkit.Chem.Fragments.fr_phenol_noOrthoHbond(mol))
    phos_acid.append(rdkit.Chem.Fragments.fr_phos_acid(mol))
    phos_ester.append(rdkit.Chem.Fragments.fr_phos_ester(mol))
    piperdine.append(rdkit.Chem.Fragments.fr_piperdine(mol))
    piperzine.append(rdkit.Chem.Fragments.fr_piperzine(mol))
    priamide.append(rdkit.Chem.Fragments.fr_priamide(mol))
    prisulfonamd.append(rdkit.Chem.Fragments.fr_prisulfonamd(mol))
    pyridine.append(rdkit.Chem.Fragments.fr_pyridine(mol))
    quatN.append(rdkit.Chem.Fragments.fr_quatN(mol))
    sulfide.append(rdkit.Chem.Fragments.fr_sulfide(mol))
    sulfonamd.append(rdkit.Chem.Fragments.fr_sulfonamd(mol))
    sulfone.append(rdkit.Chem.Fragments.fr_sulfone(mol))
    term_acetylene.append(rdkit.Chem.Fragments.fr_term_acetylene(mol))
    tetrazole.append(rdkit.Chem.Fragments.fr_tetrazole(mol))
    thiazole.append(rdkit.Chem.Fragments.fr_thiazole(mol))
    thiocyan.append(rdkit.Chem.Fragments.fr_thiocyan(mol))
    thiophene.append(rdkit.Chem.Fragments.fr_thiophene(mol))
    unbrch_alkane.append(rdkit.Chem.Fragments.fr_unbrch_alkane(mol))
    urea.append(rdkit.Chem.Fragments.fr_urea(mol))

    df1['Al_COO'] = Al_COO
    df1['Al_OH'] = Al_OH
    df1['Al_OH_noTert'] = Al_OH_noTert
    df1['ArN'] = ArN
    df1['Ar_COO'] = Ar_COO
    df1['Ar_N'] = Ar_N
    df1['Ar_NH'] = Ar_NH
    df1['Ar_OH'] = Ar_OH
    df1['COO'] = COO
    df1['COO2'] = COO2
    df1['C_O'] = C_O
    df1['C_O_noCOO'] = C_O_noCOO
    df1['C_S'] = C_S
    df1['HOCCN'] = HOCCN
    df1['Imine'] = Imine
    df1['NH0'] =NH0
    df1['NH1'] = NH1
    df1['NH2'] = NH2
    df1['N_O'] = N_O
    df1['Ndealkylation1'] = Ndealkylation1
    df1['Ndealkylation2'] = Ndealkylation2
    df1['Nhpyrrole'] = Nhpyrrole
    df1['SH'] = SH
    df1['aldehyde'] = aldehyde
    df1['alkyl_carbamate'] = alkyl_carbamate
    df1['alkyl_halide'] = alkyl_halide
    df1['allylic_oxid'] = allylic_oxid
    df1['amide'] = amide
    df1['amidine'] = amidine
    df1['aniline'] = aniline
    df1['aryl_methyl'] = aryl_methyl
    df1['azide'] = azide
    df1['azo'] = azo
    df1['barbitur'] = barbitur
    df1['benzene'] = benzene
    df1['benzodiazepine'] = benzodiazepine
    df1['bicyclic'] = bicyclic
    df1['diazo'] = diazo
    df1['dihydropyridine'] = dihydropyridine
    df1['epoxide'] = epoxide
    df1['ester'] = ester
    df1['ether'] = ether
    df1['furan'] = furan
    df1['guanido'] = guanido
    df1['halogen'] = halogen
    df1['hdrzine'] = hdrzine
    df1['hdrzone'] = hdrzone
    df1['imidazole'] = imidazole
    df1['imide'] = imide
    df1['isocyan'] = isocyan
    df1['isothiocyan'] = isothiocyan
    df1['ketone'] = ketone
    df1['ketone_Topliss'] = ketone_Topliss
    df1['lactam'] = lactam
    df1['lactone'] = lactone
    df1['methoxy'] = methoxy
    df1['morpholine'] = morpholine
    df1['nitrile'] = nitrile
    df1['nitro'] = nitro
    df1['nitro_arom'] = nitro_arom
    df1['nitro_arom_nonortho'] = nitro_arom_nonortho
    df1['nitroso'] = nitroso
    df1['oxazole'] = oxazole
    df1['oxime'] = oxime
    df1['para_hydroxylation'] = para_hydroxylation
    df1['phenol'] = phenol
    df1['phenol_noOrthoHbond'] = phenol_noOrthoHbond
    df1['phos_acid'] = phos_acid
    df1['phos_ester'] = phos_ester
    df1['piperdine'] = piperdine
    df1['piperzine'] = piperzine
    df1['priamide'] = priamide
    df1['prisulfonamd'] = prisulfonamd
    df1['pyridine'] = pyridine
    df1['quatN'] = quatN
    df1['sulfide'] = sulfide
    df1['sulfonamd'] = sulfonamd
    df1['sulfone'] = sulfone
    df1['term_acetylene'] =term_acetylene
    df1['tetrazole'] = tetrazole
    df1['thiazole'] = thiazole
    df1['thiocyan'] = thiocyan
    df1['thiophene'] = thiophene
    df1['unbrch_alkane'] = unbrch_alkane
    df1['urea'] = urea
    
    return (df1)
    
    
    
def atom_numbers(mol):    
    H = rdqueries.AtomNumEqualsQueryAtom(1)
    C = rdqueries.AtomNumEqualsQueryAtom(6)
    N = rdqueries.AtomNumEqualsQueryAtom(7)
    O = rdqueries.AtomNumEqualsQueryAtom(8)
    F = rdqueries.AtomNumEqualsQueryAtom(9)
    P = rdqueries.AtomNumEqualsQueryAtom(15)
    S = rdqueries.AtomNumEqualsQueryAtom(16)
    Cl = rdqueries.AtomNumEqualsQueryAtom(17)
    I = rdqueries.AtomNumEqualsQueryAtom(53)
    Br = rdqueries.AtomNumEqualsQueryAtom(35)
    
    nH = []
    nC = []
    nN = []
    nO = []
    nF = []
    nP = []
    nS = []
    nCl =[]
    nI = []
    nBr =[]
    
    nH.append(len(mol.GetAtomsMatchingQuery(H)))
    nC.append(len(mol.GetAtomsMatchingQuery(C)))
    nN.append(len(mol.GetAtomsMatchingQuery(N)))
    nO.append(len(mol.GetAtomsMatchingQuery(O)))
    nF.append(len(mol.GetAtomsMatchingQuery(F)))
    nP.append(len(mol.GetAtomsMatchingQuery(P)))
    nS.append(len(mol.GetAtomsMatchingQuery(S)))
    nCl.append(len(mol.GetAtomsMatchingQuery(Cl)))
    nI.append(len(mol.GetAtomsMatchingQuery(I)))
    nBr.append(len(mol.GetAtomsMatchingQuery(Br)))
    
    df2['H'] = nH
    df2['C'] = nC
    df2['N'] = nN
    df2['O'] = nO
    df2['F'] = nF
    df2['P'] = nP
    df2['S'] = nS
    df2['Cl'] = nCl
    df2['I'] = nI
    df2['Br'] = nBr
    
    return (df2)

