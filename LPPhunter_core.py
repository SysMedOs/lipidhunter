# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig
# The software is currently  under development and is not ready to be released.
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

import ConfigParser
import time

from rdkit import Chem

from LibLPPhunter.xic import XIC
from LibLPPhunter.msms import MSMS
from LibLPPhunter.plot import Spectra_Ploter

from LibLPPhunter.encode_checker import check_encode
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


for mz2get in mz2get_lst:
    if infile_type.lower() == 'mzml':
        infile_name = config.get('inputfile', 'filename')
        print infile_name
        encode_typ = check_encode(infile_name)
        print 'assume to be:', encode_typ, 'encoded'

        xic_spec = XIC(infile_name, encode_typ, ms1_precision, msn_precision)

        rt_lst, xic_df, ms_spectra_dct = xic_spec.extract_mz(mz2get, rt_range)

        print 'main peaks in range', rt_lst
        print xic_df.tail(5)
        msms_spec = MSMS(infile_name, encode_typ, ms1_precision, msn_precision)
        msms_spectra_dct = msms_spec.get_ms2(mz2get, rt_lst)

        spec_plt = Spectra_Ploter()
        spec_plt.plot_all(mz2get, xic_df, ms_spectra_dct, msms_spectra_dct, path=output_folder)

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
