# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.
# For more info please contact zhixu.ni@uni-leipzig.de

import pymzml


class XIC(object):

    def __init__(self, mzml):

        self.mzml_obj = pymzml.run.Reader(mzml, MS1_Precision=20e-6, MSn_Precision=200e-6)


    def find_mz(self, mz):

        timeDependentIntensities = []
        for spectrum in self.mzml_obj:
            print spectrum.keys()
            if spectrum['MS:1000511'] == 1:  # 'MS:1000511' == 'ms level'
                print spectrum
                for _mz, _i in spectrum.peaks:
                    print _mz, _i

                break
        #         matchList = spectrum.hasPeak(mz)
        #     if matchList != []:
        #         for mz, I in matchList:
        #             timeDependentIntensities.append([spectrum['scan time'], I, mz])
        # for rt, i, mz in timeDependentIntensities:
        #     print('{0:5.3f} {1:13.4f}       {2:10}'.format(rt, i, mz))
