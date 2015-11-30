# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

import ConfigParser
import time

from rdkit import Chem

from LibLPPhunter.XIC import XIC
from LibLPPhunter.MSMS import MSMS
from LibLPPhunter.Plot import Spectra_Ploter

from LibLPPhunter.EncodeChecker import check_encode
from LibLPPhunter.ExactMassCalc import Elem2Mass


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
for _sdf in sdf_obj:
    _formula = _sdf.GetProp('CHEMICAL_FORMULA')
    _formula_COOH = ''.join([_formula, 'COOH'])
    _formula_COOH = mzcalc.get_elem(_formula_COOH)
    _mz_COOH = mzcalc.get_mass(_formula_COOH)
    if _mz_COOH not in mz2get_lst:
        mz2get_lst.append(_mz_COOH)
    else:
        pass

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
    spec_plt.plot_all(mz2get_lst, xic_dct, ms_spectra_dct, msms_spectra_dct, path=output_folder)

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
