# -*- coding: utf-8 -*-
# Copyright 2014-16 Zhixu Ni, AG Bioanalytik,BBZ,University of Leipzig          
# The software is currently  under development and is not ready to be released. 
# A suitable license will be chosen before the official release.               
# For more info please contact zhixu.ni@uni-leipzig.de

import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw


class MolReader(object):

    def __init__(self):
        self.pc_hg_smi = r'OP(=O)([O-])OCC[N+](C)(C)C'
        self.pc_hg_mol = Chem.MolFromSmiles(self.pc_hg_smi)
        AllChem.Compute2DCoords(self.pc_hg_mol)

    def get_pl_backbone(self, usr_mol):

        # Get std_smi
        AllChem.Compute2DCoords(usr_mol)

        std_smi = Chem.MolToSmiles(usr_mol)

        # make regular expressions to check PL substructures

        checker_1_std = re.compile(r'(?P<sn1>\S*C[(][=]O[)]O)(\S*)(?P<sn2>OC[(][=]O[)][^P^p]*)([)]CO\S*)')
        # std_smi = r'CCCCCCCCCCCCCCCC(=O)OC[C](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC'
        checker_2_std = re.compile(r'(?P<sn1>\S*C[(][=]O[)]O)([COPNSHcopnsh\[\]\(\)\+\-\=]*)'
                                   r'(?P<sn2>OC[(][=]O[)][COHcoh\[\]\(\)\+\-\=]*)')
        # std_smi = r'CCCCCC[C]=[C]CCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCC)COP(=O)([O-])OCC[N+](C)(C)C'
        checker_3_std = re.compile(r'(?P<sn1>\S*C[(][=]O[)]O)(C[(]C)(?P<sn2>OC[(][=]O[)]\S*)([)]COP\S*)')
        # std_smi = r'CCCCCCCC[C]=[C]CCCCCCCC(=O)OC[C](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC[C]=[C]CCCCCCCC'
        checker_4_std = re.compile(r'(?P<sn1>\S*C[(][=]O[)]O)([COPNSHcopnsh\[\]\(\)\+\-\=]*)(?P<sn2>OC[(][=]O[)]\S*)')
        # std_smi = r'CCCCCC[C]=[C]CCCCCCCC(=O)OC(COC(=O)CCCCCCCCCCCCCCC)COP(=O)([O-])OCC[N+](C)(C)C'
        checker_lst = [checker_1_std, checker_2_std, checker_3_std, checker_4_std]
        try:

            for _checker in checker_lst:
                # print('checker', checker_lst.index(_checker) + 1)
                _pl_dct = self._re_checker(_checker, std_smi, usr_mol)
                if len(_pl_dct.keys()) > 0:
                    _h_mass = 1.007825
                    _pl_dct['sn1_mz'] = Descriptors.ExactMolWt(Chem.MolFromSmiles(_pl_dct['sn1_smi'])) - _h_mass
                    _pl_dct['sn2_mz'] = Descriptors.ExactMolWt(Chem.MolFromSmiles(_pl_dct['sn2_smi'])) - _h_mass
                    _pl_dct['hg_mz'] = Descriptors.ExactMolWt(Chem.MolFromSmiles(_pl_dct['hg_smi'])) - _h_mass
                    _pl_dct['rest_mz'] = Descriptors.ExactMolWt(Chem.MolFromSmiles(_pl_dct['rest_smi'])) - _h_mass
                    break
                else:
                    pass
        except:
            print('Not pass --->', std_smi)
            _pl_dct = {'sn1_smi': '', 'sn2_smi': '', 'rest_smi': '', 'hg_smi': '',
                       'sn1_mz': '', 'sn2_mz': '', 'rest_mz': '', 'hg_mz': ''}

        return _pl_dct

    def _re_checker(self, usr_patt, usr_smi, usr_mol):

        _pl_dct = {}

        # checker_std = re.compile(r'(?P<sn1>\S*C[(][=]O[)]O)(\S*)(?P<sn2>OC[(][=]O[)]\S*)')
        checker = re.match(usr_patt, usr_smi)

        if checker:

            pl_dct = checker.groupdict()
            _sn1 = pl_dct['sn1']
            _sn2 = pl_dct['sn2']

            _rest_lst = list(checker.groups())

            rebuild_str = ''.join(_rest_lst)

            if rebuild_str == usr_smi:

                _rest_lst.remove(_sn1)
                _rest_lst.remove(_sn2)
                _rest_str_pre = ''.join(_rest_lst)
                # print 'sn1', _sn1
                # print 'sn2', _sn2
                # print '_rest_lst', _rest_lst
                # print _rest_str_pre

                # remove '[C]' from smiles
                _bad_C_std = re.compile(r'[\[]C[\]]')
                _rest_C_str = _bad_C_std.split(_rest_str_pre)
                if len(_rest_C_str) > 1:
                    _rest_str = 'C'.join(_rest_C_str)
                    # print '[C] replaced--->'
                else:
                    _rest_str = _rest_str_pre
                    # print 'No [C]'
                # print _rest_str

                _rest_mol = Chem.MolFromSmiles(_rest_str)
                AllChem.Compute2DCoords(_rest_mol)

                if usr_mol.HasSubstructMatch(self.pc_hg_mol):
                    # _match_hg = usr_mol.GetSubstructMatch(_rest_mol)
                    if usr_mol.HasSubstructMatch(_rest_mol):
                        _match_hg = usr_mol.GetSubstructMatch(self.pc_hg_mol)
                        # _match = usr_mol.GetSubstructMatch(_rest_mol)

                        _pl_dct['sn1_smi'] = _sn1
                        _pl_dct['sn2_smi'] = _sn2
                        _pl_dct['rest_smi'] = _rest_str
                        _pl_dct['hg_smi'] = self.pc_hg_smi

                        # save structure as images
                        # _hmdb_id = usr_mol.GetProp('HMDB_ID')
                        # _img_fp = _hmdb_id + '.png'
                        #
                        # img = Draw.MolToImage(usr_mol, highlightAtoms=_match_hg, size=(600, 600))
                        # img.save(_img_fp)

                        return _pl_dct

                    else:
                        print '_rest do not contain HG'
                        return _pl_dct
                        # _rest_img = Draw.MolToImage(_rest_mol)
                        # _rest_img.save(_img_fp)
                else:
                    # _rest do not contain HG
                    print '_rest do not contain HG'
                    return _pl_dct
            else:

                print 'less groups'
                return _pl_dct

        else:
            # print 'usr_smi not fit', usr_smi
            return _pl_dct


# usr_sdf = r'newsdf.sdf'
# sdf_obj = Chem.SDMolSupplier(usr_sdf)
# sdf = SDFreader()
# i = 1
# for m in sdf_obj:
#     sdf.get_pl_backbone(m)
#     print i
#     i += 1
#     # break



