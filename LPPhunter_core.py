# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

import ConfigParser
import time

import pandas as pd

from rdkit import Chem

from LibLPPhunter.XIC import XIC
from LibLPPhunter.MSMS import MSMS
from LibLPPhunter.Plot import Spectra_Ploter

from LibLPPhunter.EncodeChecker import check_encode
from LibLPPhunter.ExactMassCalc import Elem2Mass

from LibLPPhunter.SDFparser import MolReader


print('Start --->')
config = ConfigParser.ConfigParser()
config.read('config.ini')
st_time = time.clock()

infile_type = config.get('inputfile', 'filetype')
ms1_precision = config.get('inputfile', 'ms1_precision')
msn_precision = config.get('inputfile', 'msn_precision')

output_folder = config.get('output', 'output_folder')

sdf_f = 'newsdf.sdf'

rt_range = [10.0, 30.0]

# mz2get = 778.560
# mz2get_lst = [778.560, 776.544, 774.529, 806.591]

mz2get_lst = []
sdf_obj = Chem.SDMolSupplier(sdf_f)
mzcalc = Elem2Mass()
sdf_dct = {'hmdb_id': [], 'formula': [], 'pr_mz': [],
           'hg_smi': [], 'hg_mz': [],
           'sn1_smi': [], 'sn1_mz': [],
           'sn2_smi': [], 'sn2_mz': []}

print('Reading SDF--->')
molreader = MolReader()

for _mol in sdf_obj:
    _hmdb_id = _mol.GetProp('HMDB_ID')
    _formula = _mol.GetProp('CHEMICAL_FORMULA')
    _formula_COOH = ''.join([_formula, 'COOH'])
    _formula_COOH = mzcalc.get_elem(_formula_COOH)
    _mz_COOH = mzcalc.get_mass(_formula_COOH)

    _pl_dct = molreader.get_pl_backbone(_mol)

    sdf_dct['hmdb_id'].append(_hmdb_id)
    sdf_dct['formula'].append(_formula)
    sdf_dct['pr_mz'].append(_mz_COOH)
    sdf_dct['hg_smi'].append(_pl_dct['hg_smi'])
    sdf_dct['sn1_smi'].append(_pl_dct['sn1_smi'])
    sdf_dct['sn2_smi'].append(_pl_dct['sn2_smi'])
    sdf_dct['hg_mz'].append(_pl_dct['hg_mz'])
    sdf_dct['sn1_mz'].append(_pl_dct['sn1_mz'])
    sdf_dct['sn2_mz'].append(_pl_dct['sn2_mz'])

    if _mz_COOH not in mz2get_lst:
        mz2get_lst.append(_mz_COOH)
    else:
        pass

sdf_df = pd.DataFrame(data=sdf_dct)

mz2get_lst.sort()

print len(mz2get_lst), ' different m/z need to find===>'


# for mz2get in mz2get_lst:
if infile_type.lower() == 'mzml':
    infile_name = config.get('inputfile', 'filename')
    print infile_name
    encode_typ = check_encode(infile_name)
    print 'assume to be:', encode_typ, 'encoded'

    xic_spec = XIC(infile_name, encode_typ, ms1_precision, msn_precision)

    rt_dct, xic_dct, ms_spectra_dct = xic_spec.extract_mz(mz2get_lst, rt_range)
    #
    # print 'main peaks in range', rt_lst
    # print xic_dct.keys()
    # print ms_spectra_dct.keys()
    # for x in xic_dct.keys():
    #     print type(xic_dct[x])
    #     print xic_dct[x].head()
    # for y in rt_dct.keys():
    #     print(rt_dct[y])
    #     break
    # for x in ms_spectra_dct.keys():
    #     print type(ms_spectra_dct[x])

    msms_spec = MSMS(infile_name, encode_typ, ms1_precision, msn_precision)
    msms_spectra_dct = msms_spec.get_ms2(mz2get_lst, rt_dct)

    # print msms_spectra_dct.keys()
    # for x in msms_spectra_dct.keys():
    #     print msms_spectra_dct[x]

    spec_plt = Spectra_Ploter()
    spec_plt.plot_all(mz2get_lst, sdf_df, xic_dct, ms_spectra_dct, msms_spectra_dct, path=output_folder)

# if infile_type.lower() == 'mzml.gz':
#     infile_name = config.get('inputfile', 'filename')
#     print infile_name
#
#     xic_spec = XIC(infile_name)
#     xic_spec.extract_mz(778.560, ppm=20)

else:
    print 'No input mzML'

ed_time = time.clock() - st_time
print ed_time
