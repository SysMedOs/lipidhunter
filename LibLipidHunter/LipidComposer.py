# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2019  SysMedOs_team @ AG Bioanalytik, University of Leipzig:
# SysMedOs_team: Zhixu Ni, Georgia Angelidou, Mike Lange, Maria Fedorova
# LipidHunter is Dual-licensed
#     For academic and non-commercial use: `GPLv2 License` Please read more information by the following link:
#         [The GNU General Public License version 2] (https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
#     For commercial use:
#         please contact the SysMedOs_team by email.
# Please cite our publication in an appropriate form.
# Ni, Zhixu, Georgia Angelidou, Mike Lange, Ralf Hoffmann, and Maria Fedorova.
# "LipidHunter identifies phospholipids by high-throughput processing of LC-MS and shotgun lipidomics datasets."
# Analytical Chemistry (2017).
# DOI: 10.1021/acs.analchem.7b01126
#
# For more info please contact:
#     Developer Zhixu Ni zhixu.ni@uni-leipzig.de
#     Developer Georgia Angelidou georgia.angelidou@uni-leipzig.de

import itertools
from typing import List

import pandas as pd

from LibLipidHunter.LipidNomenclature import NameParserFA
from LibLipidHunter.AbbrElemCalc import ElemCalc
from LibLipidHunter.ParallelFunc import ppm_window_para


class LipidComposer:
    """
    Lipid composer is designed to generate all possible combinations of given lipid class based on the
    FA white list from user settings. and calculate all information required for MS/MS identification
    including elemental composition, m/z of precursor and fragments ion under selected charge status
    as fragment or from neutral losses.
    """

    def __init__(self):
        """
        Default parameters for the elemental compositions
        """
        pa_hg_elem = {"C": 0, "H": 3, "O": 4, "P": 1, "N": 0}
        pc_hg_elem = {"C": 5, "H": 14, "O": 4, "P": 1, "N": 1}
        pe_hg_elem = {"C": 2, "H": 8, "O": 4, "P": 1, "N": 1}
        pg_hg_elem = {"C": 3, "H": 9, "O": 6, "P": 1, "N": 0}
        pi_hg_elem = {"C": 6, "H": 13, "O": 9, "P": 1, "N": 0}
        pip_hg_elem = {"C": 6, "H": 14, "O": 12, "P": 2, "N": 0}
        ps_hg_elem = {"C": 3, "H": 8, "O": 6, "P": 1, "N": 1}
        tg_hg_elem = {"C": 0, "H": 0, "O": 0, "P": 0, "N": 0}
        # Todo: add SM corresponding HG sections around here

        self.lipid_hg_lst = ["PA", "PC", "PE", "PG", "PI", "PS", "PIP", "TG"]

        self.lipid_hg_elem_dct = {
            "PA": pa_hg_elem,
            "PC": pc_hg_elem,
            "PE": pe_hg_elem,
            "PG": pg_hg_elem,
            "PI": pi_hg_elem,
            "PS": ps_hg_elem,
            "PIP": pip_hg_elem,
            "TG": tg_hg_elem,
            "DG": tg_hg_elem,
        }

        self.glycerol_bone_elem_dct = {"C": 3, "H": 2}
        self.link_o_elem_dct = {"O": -1, "H": 2}
        self.link_p_elem_dct = {"O": -1}
        self.elem_dct = {
            "H": [1.0078250321, 0.999885],
            "D": [2.0141017780, 0.0001157],
            "C": [12.0, 0.9893],
            "N": [14.0030740052, 0.99632],
            "O": [15.9949146221, 0.99757],
            "Na": [22.98976967, 1.0],
            "P": [30.97376151, 1.0],
            "S": [31.97207069, 0.9493],
            "K": [38.9637069, 0.932581],
        }

    @staticmethod
    def calc_fa_df(lipid_class: str, fa_df: pd.DataFrame) -> List[List[str]]:
        """
        Read list of FA abbreviations from FA white list for lipid enumeration later
        Args:
            lipid_class(str): lipid class abbreviation
            fa_df(pd.DataFrame): FA white list in pd.DataFrame

        Returns:
            sn_units_lst(List[List[str]]): list of FA for each FA position
        """

        sn_units_lst = []

        header_lst = fa_df.columns.values.tolist()

        if lipid_class in ["PA", "PC", "PE", "PG", "PI", "PS"]:
            if "PL" in header_lst and "FATTYACID" in header_lst:
                pl_fa1_df = fa_df.query('PL == "T" and FA1 == "T"')
                pl_fa2_df = fa_df.query('PL == "T" and FA2 == "T"')

                pl_fa1_lst = pl_fa1_df["FATTYACID"].tolist()
                pl_fa2_lst = pl_fa2_df["FATTYACID"].tolist()

                sn_units_lst = [pl_fa1_lst, pl_fa2_lst]

        elif lipid_class in ["LPA", "LPC", "LPE", "LPG", "LPI", "LPS"]:
            if "LPL" in header_lst and "FATTYACID" in header_lst:
                pl_fa1_df = fa_df.query('LPL == "T" and FA1 == "T"')
                pl_fa1_lst = pl_fa1_df["FATTYACID"].tolist()
                sn_units_lst = [pl_fa1_lst, [""]]
            elif "PL" in header_lst and "FATTYACID" in header_lst:
                pl_fa1_df = fa_df.query('PL == "T" and FA1 == "T"')
                pl_fa1_lst = pl_fa1_df["FATTYACID"].tolist()
                sn_units_lst = [pl_fa1_lst, [""]]

        elif lipid_class in ["TG"]:
            if "TG" in header_lst and "FATTYACID" in header_lst:
                tg_fa1_df = fa_df.query('TG == "T" and FA1 == "T"')
                tg_fa2_df = fa_df.query('TG == "T" and FA2 == "T"')
                tg_fa3_df = fa_df.query('TG == "T" and FA3 == "T"')

                tg_fa1_lst = tg_fa1_df["FATTYACID"].tolist()
                tg_fa2_lst = tg_fa2_df["FATTYACID"].tolist()
                tg_fa3_lst = tg_fa3_df["FATTYACID"].tolist()

                sn_units_lst = [tg_fa1_lst, tg_fa2_lst, tg_fa3_lst]
        elif lipid_class in ["DG"]:
            if "DG" in header_lst and "FATTYACID" in header_lst:
                dg_fa1_df = fa_df.query('TG == "T" and FA1 == "T"')
                dg_fa2_df = fa_df.query('TG == "T" and FA2 == "T"')

                dg_fa1_lst = dg_fa1_df["FATTYACID"].tolist()
                dg_fa2_lst = dg_fa2_df["FATTYACID"].tolist()

                sn_units_lst = [dg_fa1_lst, dg_fa2_lst]
        elif lipid_class in ["CL"]:
            if "CL" in header_lst and "FATTYACID" in header_lst:
                cl_fa1_df = fa_df.query('CL == "T" and FA1 == "T"')
                cl_fa2_df = fa_df.query('CL == "T" and FA2 == "T"')
                cl_fa3_df = fa_df.query('CL == "T" and FA3 == "T"')
                cl_fa4_df = fa_df.query('CL == "T" and FA4 == "T"')

                cl_fa1_lst = cl_fa1_df["FATTYACID"].tolist()
                cl_fa2_lst = cl_fa2_df["FATTYACID"].tolist()
                cl_fa3_lst = cl_fa3_df["FATTYACID"].tolist()
                cl_fa4_lst = cl_fa4_df["FATTYACID"].tolist()

                sn_units_lst = [cl_fa1_lst, cl_fa2_lst, cl_fa3_lst, cl_fa4_lst]

        return sn_units_lst

    def calc_fa_query(
        self, lipid_class: str, fa_whitelist: str, ms2_ppm: int = 100
    ) -> pd.DataFrame:
        """
        Prepare all query strings for the MS/MS identification
        build all pandas query code in advance to save the time during identification
        Args:
            lipid_class(str): Lipid class abbreviation
            fa_whitelist(str): the ile path to FA white list
            ms2_ppm(int): ppm at MS/MS level, set default to 100 for some qTOF instrument

        Returns (pd.DataFrame): The calculated information including m/z and query code in DataFrame

        """

        usr_fa_df = pd.read_excel(fa_whitelist)
        usr_fa_df = usr_fa_df.fillna(value="F")
        tmp_columns = usr_fa_df.columns.tolist()
        usr_fa_df.columns = usr_fa_df.columns.str.upper()

        if lipid_class in ["PL", "PA", "PC", "PE", "PG", "PI", "PS", "SM"]:
            if lipid_class in tmp_columns:
                pass
            elif "PL" in tmp_columns:
                pass
            else:
                return False
        elif lipid_class in ["LPL", "LPA", "LPC", "LPE", "LPG", "LPI", "LPS"]:
            if lipid_class in tmp_columns:
                pass
            elif "LPL" in tmp_columns or "PL" in tmp_columns:
                pass
            else:
                return False
        elif lipid_class in ["TG", "DG"]:
            if lipid_class in tmp_columns:
                pass
            else:
                return False
        else:
            return False

        sn_units_lst = self.calc_fa_df(
            lipid_class, usr_fa_df
        )  # Return a list with list of the FA for each sn position
        fa_abbr_lst = []
        # For PL lem(sn_units_lst) = 2 and for TG len(sn_units_lst) = 3
        for _s in sn_units_lst:
            fa_abbr_lst.extend(_s)  # Compine all the FA in one list
        fa_abbr_lst = sorted(list(set(fa_abbr_lst)))

        abbr_parser = NameParserFA()
        elem_calc = ElemCalc()
        usr_fa_dct = {}
        for _fa_abbr in fa_abbr_lst:
            if _fa_abbr:
                _fa_info_dct = abbr_parser.get_fa_info(
                    _fa_abbr
                )  # Calculate all the information for each FA
                _lipid_formula, _lipid_elem_dct = elem_calc.get_formula(
                    _fa_abbr
                )  # get the elemental composition of FA
                # add the Abbr, formula and the exact mass in the dictionary
                _fa_info_dct["ABBR"] = _fa_abbr
                _fa_info_dct["FORMULA"] = _lipid_formula
                _fa_info_dct["EXACTMASS"] = elem_calc.get_exactmass(
                    _lipid_elem_dct
                )  # Calc. the exact mass for each FA
                usr_fa_dct[_fa_abbr] = _fa_info_dct

        usr_fa_df = pd.DataFrame(
            usr_fa_dct
        ).T.copy()  # put all the info for the FA in a dataframe
        # usr_fa_df.is_copy = False

        # create the queries for the FA fragments and MG
        for _fa_ion in ["[FA-H]-", "[FA-H2O-H]-", "[FA-H2O+H]+"]:
            usr_fa_df["%s_MZ_LOW" % _fa_ion] = ppm_window_para(
                usr_fa_df["%s_MZ" % _fa_ion].values.tolist(), ms2_ppm * -1
            )
            usr_fa_df["%s_MZ_HIGH" % _fa_ion] = ppm_window_para(
                usr_fa_df["%s_MZ" % _fa_ion].values.tolist(), ms2_ppm
            )
            usr_fa_df["%s_Q" % _fa_ion] = (
                usr_fa_df["%s_MZ_LOW" % _fa_ion].astype(str)
                + " <= mz <= "
                + usr_fa_df["%s_MZ_HIGH" % _fa_ion].astype(str)
            )

        # More specific fragments for PL
        if lipid_class in ["PA", "PC", "PE", "PG", "PI", "PS", "PIP"]:

            # there alreadz the backbone of the stucture eg. for PE C5H10NO4P
            # for the [LPE(16:0)-H]- missing the exact mass of the FA
            # for the [LPE(16:0)-H2O-H]- we have the loss of water which already calculated in the structure
            # thats why addition of the FA but without of the water
            lyso_type_dct = {
                "[L%s-H]-" % lipid_class: "EXACTMASS",
                "[L%s-H2O-H]-" % lipid_class: "[FA-H2O]_MZ",
            }

            # backbone creation for the different PL
            lyso_base_elem_dct = self.lipid_hg_elem_dct[lipid_class]
            for _e in list(self.glycerol_bone_elem_dct.keys()):
                lyso_base_elem_dct[_e] += self.glycerol_bone_elem_dct[_e]

            # the element here is with no Hydroxyl on fa1 and fa2, here a [M-H]- is already considered
            lyso_base_mz = (
                elem_calc.get_exactmass(lyso_base_elem_dct)
                + 1.0078250321
                + 15.9949146221
            )
            # Todo: SM section
            if lipid_class in ["PC", "SM"]:
                lyso_base_mz -= (
                    12.0 + 2 * 1.0078250321
                )  # LPC loss one -CH3 from HG (one H already remove above)

            for _lyso_ion in list(lyso_type_dct.keys()):

                if lipid_class in ["PC", "SM"]:

                    if lyso_type_dct[_lyso_ion] == "EXACTMASS":
                        usr_fa_df["%s_ABBR" % _lyso_ion] = (
                            "[L"
                            + lipid_class
                            + "("
                            + usr_fa_df["ABBR"].str.strip("FA")
                            + ")-CH3]-"
                        )
                    elif lyso_type_dct[_lyso_ion] == "[FA-H2O]_MZ":
                        usr_fa_df["%s_ABBR" % _lyso_ion] = (
                            "[L"
                            + lipid_class
                            + "("
                            + usr_fa_df["ABBR"].str.strip("FA")
                            + ")-H2O-CH3]-"
                        )
                    else:
                        usr_fa_df["%s_ABBR" % _lyso_ion] = "ERROR"
                else:
                    if lyso_type_dct[_lyso_ion] == "EXACTMASS":
                        usr_fa_df["%s_ABBR" % _lyso_ion] = (
                            "[L"
                            + lipid_class
                            + "("
                            + usr_fa_df["ABBR"].str.strip("FA")
                            + ")-H]-"
                        )
                    elif lyso_type_dct[_lyso_ion] == "[FA-H2O]_MZ":
                        usr_fa_df["%s_ABBR" % _lyso_ion] = (
                            "[L"
                            + lipid_class
                            + "("
                            + usr_fa_df["ABBR"].str.strip("FA")
                            + ")-H2O-H]-"
                        )
                    else:
                        usr_fa_df["%s_ABBR" % _lyso_ion] = "ERROR"

                usr_fa_df["%s_MZ" % _lyso_ion] = (
                    lyso_base_mz + usr_fa_df[lyso_type_dct[_lyso_ion]]
                )
                usr_fa_df["%s_MZ_LOW" % _lyso_ion] = ppm_window_para(
                    usr_fa_df["%s_MZ" % _lyso_ion].values.tolist(), ms2_ppm * -1
                )
                usr_fa_df["%s_MZ_HIGH" % _lyso_ion] = ppm_window_para(
                    usr_fa_df["%s_MZ" % _lyso_ion].values.tolist(), ms2_ppm
                )
                usr_fa_df["%s_Q" % _lyso_ion] = (
                    usr_fa_df["%s_MZ_LOW" % _lyso_ion].astype(str)
                    + " <= mz <= "
                    + usr_fa_df["%s_MZ_HIGH" % _lyso_ion].astype(str)
                )
        elif lipid_class in ["TG"]:
            # Cannot calculate the theoritical m/z values of the DG fragments sinc we calculate the loss of a FA
            # and we dont know the cobination of the remaining 2
            # TODO(georgia.angelidou@uni-leipzig.de): create the section for theuniue fragments when there is TG
            mg_type_dct = {"[MG-H2O+H]+": "EXACTMASS"}
            mg_base_elem_dct = self.lipid_hg_elem_dct[lipid_class]

            for _e in self.glycerol_bone_elem_dct.keys():
                mg_base_elem_dct[_e] += self.glycerol_bone_elem_dct[_e]

            # Calculate the rest of monoglycerol after the neutral loss of 2 FA in protonated form (one without Water)
            mg_base_elem_dct = (
                elem_calc.get_exactmass(mg_base_elem_dct)
                + 15.9949146221
                + (3 * 1.0078250321)
            )
            for _mg_ion in mg_type_dct.keys():
                usr_fa_df["%s_ABBR" % _mg_ion] = (
                    "[MG(" + usr_fa_df["ABBR"].str.strip("FA") + ")-H2O+H]+"
                )
                usr_fa_df["%s_MZ" % _mg_ion] = (
                    mg_base_elem_dct + usr_fa_df[mg_type_dct[_mg_ion]]
                )
                usr_fa_df["%s_MZ_LOW" % _mg_ion] = ppm_window_para(
                    usr_fa_df["%s_MZ" % _mg_ion].values.tolist(), ms2_ppm * -1
                )
                usr_fa_df["%s_MZ_HIGH" % _mg_ion] = ppm_window_para(
                    usr_fa_df["%s_MZ" % _mg_ion].values.tolist(), ms2_ppm
                )
                usr_fa_df["%s_Q" % _mg_ion] = (
                    usr_fa_df["%s_MZ_LOW" % _mg_ion].astype(str)
                    + " <= mz <= "
                    + usr_fa_df["%s_MZ_HIGH" % _mg_ion].astype(str)
                )

        elif lipid_class in ["DG"]:
            mg_type_dct = {"[MG-H2O+H]+": "EXACTMASS"}
            mg_base_elem_dct = self.lipid_hg_elem_dct[lipid_class]

            for _e in self.glycerol_bone_elem_dct.keys():
                mg_base_elem_dct[_e] += self.glycerol_bone_elem_dct[_e]

            # Calculate the rest of monoglycerol after the neutral loss of 2 FA in protonated form (one without Water)
            mg_base_elem_dct = (
                elem_calc.get_exactmass(mg_base_elem_dct)
                + 15.9949146221
                + (3 * 1.0078250321)
            )
            for _mg_ion in mg_type_dct.keys():
                usr_fa_df["%s_ABBR" % _mg_ion] = (
                    "[MG(" + usr_fa_df["ABBR"].str.strip("FA") + ")-H2O+H]+"
                )
                usr_fa_df["%s_MZ" % _mg_ion] = (
                    mg_base_elem_dct + usr_fa_df[mg_type_dct[_mg_ion]]
                )
                usr_fa_df["%s_MZ_LOW" % _mg_ion] = ppm_window_para(
                    usr_fa_df["%s_MZ" % _mg_ion].values.tolist(), ms2_ppm * -1
                )
                usr_fa_df["%s_MZ_HIGH" % _mg_ion] = ppm_window_para(
                    usr_fa_df["%s_MZ" % _mg_ion].values.tolist(), ms2_ppm
                )
                usr_fa_df["%s_Q" % _mg_ion] = (
                    usr_fa_df["%s_MZ_LOW" % _mg_ion].astype(str)
                    + " <= mz <= "
                    + usr_fa_df["%s_MZ_HIGH" % _mg_ion].astype(str)
                )
        else:
            # Todo: SM section like above
            pass
        return usr_fa_df

    def gen_all_comb(
        self, lipid_class: str, usr_fa_df: pd.DataFrame, position: bool = False
    ) -> dict:

        """
        Generate all possible FA combination of selected lipid class in discrete form without mirrored duplicates

        Args:
            lipid_class(str): Lipid class abbreviation
            usr_fa_df(pd.DataFrame): FA white list in DataFrame
            position(bool): set as False by default to consider discrete form without position specific isomers

        Returns(dict): all combinations stored in dict

        """

        fa_units_lst = self.calc_fa_df(lipid_class, usr_fa_df)

        if (
            lipid_class in ["PA", "PC", "PE", "PG", "PI", "PS", "DG", "SM"]
            and len(fa_units_lst) == 2
        ):
            fa_comb_lst = list(itertools.product(fa_units_lst[0], fa_units_lst[1]))
            fa_df_header_lst = ["FA1", "FA2"]
            # lipid_template = '{}'
        elif (
            lipid_class in ["LPA", "LPC", "LPE", "LPG", "LPI", "LPS"]
            and len(fa_units_lst) == 2
        ):
            fa_comb_lst = list(itertools.product(fa_units_lst[0], fa_units_lst[1]))
            fa_df_header_lst = ["FA1", "FA2"]
            # lipid_template = '{}'
        elif lipid_class == "TG" and len(fa_units_lst) == 3:
            fa_comb_lst = list(
                itertools.product(fa_units_lst[0], fa_units_lst[1], fa_units_lst[2])
            )
            fa_df_header_lst = ["FA1", "FA2", "FA3"]
        elif lipid_class == "CL" and len(fa_units_lst) == 4:
            fa_comb_lst = list(
                itertools.product(
                    fa_units_lst[0], fa_units_lst[1], fa_units_lst[2], fa_units_lst[3]
                )
            )
            fa_df_header_lst = ["FA1", "FA2", "FA3", "FA4"]
        else:
            fa_comb_lst = []
            fa_df_header_lst = []

        fa_combo_df = pd.DataFrame(data=fa_comb_lst, columns=fa_df_header_lst)

        fa_combo_df["CLASS"] = lipid_class

        # Todo: check SM section below:
        if lipid_class in ["PA", "PC", "PE", "PG", "PI", "PS", "DG", "SM"]:

            fa_combo_link_df = fa_combo_df.copy()
            # fa_combo_link_df.is_copy = False
            fa_combo_link_df["LINK"] = fa_combo_link_df["FA1"].str[0:2]
            fa_link_df = fa_combo_link_df[fa_combo_link_df["LINK"] == "FA"].copy()

            # fa_link_df.is_copy = False
            fa_link_df.drop(["LINK", "CLASS"], axis=1, inplace=True)
            # fa_link_df.values.argsort(kind='mergesort')
            # fa_link_df.drop(columns=['CLASS'], inplace=True)
            fa_link_df.values.sort(kind="mergesort")  # safe sort by numpy
            fa_link_df["CLASS"] = lipid_class
            fa_link_df["DISCRETE_ABBR"] = (
                fa_link_df["CLASS"]
                + "("
                + fa_link_df["FA1"].str.strip("FA")
                + "_"
                + fa_link_df["FA2"].str.strip("FA")
                + ")"
            )
            fa_link_df.sort_values(by="DISCRETE_ABBR", inplace=True)

            if lipid_class in ["PC", "PE"]:
                op_link_df = fa_combo_link_df[
                    (fa_combo_link_df["LINK"] == "O-")
                    | (fa_combo_link_df["LINK"] == "P-")
                ].copy()
                if not op_link_df.empty:
                    # op_link_df.is_copy = False
                    op_link_df.drop(["LINK"], axis=1, inplace=True)
                    op_link_df["DISCRETE_ABBR"] = (
                        op_link_df["CLASS"]
                        + "("
                        + op_link_df["FA1"].str.strip("FA")
                        + "_"
                        + op_link_df["FA2"].str.strip("FA")
                        + ")"
                    )
                    op_link_df.sort_values(by="DISCRETE_ABBR", inplace=True)

                    fa_combo_df = fa_link_df.append(op_link_df)
                    del op_link_df
            else:
                fa_combo_df = fa_link_df

            del fa_combo_link_df
            del fa_link_df
            print(
                "[INFO] --> Number of predicted lipids (exact position): ",
                fa_combo_df.shape[0],
            )

        elif lipid_class in ["LPA", "LPC", "LPE", "LPG", "LPI", "LPS"]:

            fa_combo_link_df = fa_combo_df.copy()
            # fa_combo_link_df.is_copy = False

            fa_combo_link_df["LINK"] = fa_combo_link_df["FA1"].str[0:2]
            fa_link_df = fa_combo_link_df[fa_combo_link_df["LINK"] == "FA"].copy()

            # fa_link_df.is_copy = False
            fa_link_df.drop(["LINK", "CLASS"], axis=1, inplace=True)
            # fa_link_df.values.argsort(kind='mergesort')
            # fa_link_df.drop(columns=['CLASS'], inplace=True)
            # fa_link_df.values.sort(kind='mergesort')  # safe sort by numpy
            fa_link_df["CLASS"] = lipid_class

            fa_link_df["DISCRETE_ABBR"] = (
                fa_link_df["CLASS"] + "(" + fa_link_df["FA1"].str.strip("FA") + ")"
            )
            fa_link_df.sort_values(by="DISCRETE_ABBR", inplace=True)

            if lipid_class in ["LPC", "LPE"]:
                op_link_df = fa_combo_link_df[
                    (fa_combo_link_df["LINK"] == "O-")
                    | (fa_combo_link_df["LINK"] == "P-")
                ].copy()
                if not op_link_df.empty:
                    # op_link_df.is_copy = False
                    op_link_df.drop(["LINK"], axis=1, inplace=True)
                    op_link_df["DISCRETE_ABBR"] = (
                        op_link_df["CLASS"]
                        + "("
                        + op_link_df["FA1"].str.strip("FA")
                        + ")"
                    )
                    op_link_df.sort_values(by="DISCRETE_ABBR", inplace=True)

                    fa_combo_df = fa_link_df.append(op_link_df)
                    del op_link_df
            else:
                fa_combo_df = fa_link_df

            del fa_combo_link_df
            del fa_link_df
            print(
                "[INFO] --> Number of predicted lipids (exact position): ",
                fa_combo_df.shape[0],
            )

        elif lipid_class in ["TG"]:
            fa_combo_link_df = fa_combo_df.copy()
            # fa_combo_link_df.is_copy = False
            fa_combo_link_df["LINK"] = fa_combo_link_df["FA1"].str[0:2]
            fa_link_df = fa_combo_link_df[fa_combo_link_df["LINK"] == "FA"]

            # fa_link_df.is_copy = False
            fa_link_df.drop(["LINK"], axis=1, inplace=True)
            fa_link_df.values.sort(kind="mergesort")  # safe sort by numpy
            fa_link_df["DISCRETE_ABBR"] = (
                fa_link_df["CLASS"]
                + "("
                + fa_link_df["FA1"].str.strip("FA")
                + "_"
                + fa_link_df["FA2"].str.strip("FA")
                + "_"
                + fa_link_df["FA3"].str.strip("FA")
                + ")"
            )
            fa_link_df.sort_values(by="DISCRETE_ABBR", inplace=True)
            op_link_df = fa_combo_link_df[
                (fa_combo_link_df["LINK"] == "O-") | (fa_combo_link_df["LINK"] == "P-")
            ].copy()
            if not op_link_df.empty:
                # op_link_df.is_copy = False
                op_link_df.drop(["LINK"], axis=1, inplace=True)
                op_link_df["DISCRETE_ABBR"] = (
                    op_link_df["CLASS"]
                    + "("
                    + op_link_df["FA1"].str.strip("FA")
                    + "_"
                    + op_link_df["FA2"].str.strip("FA")
                    + "_"
                    + op_link_df["FA3"].str.strip("FA")
                    + ")"
                )
                op_link_df.sort_values(by="DISCRETE_ABBR", inplace=True)
                fa_combo_df = fa_link_df.append(op_link_df)
                del op_link_df
            else:
                fa_combo_df = fa_link_df

            del fa_link_df
            del fa_combo_link_df
            print(
                "[INFO] --> Number of predicted lipids (exact position): ",
                fa_combo_df.shape[0],
            )

        elif lipid_class in ["CL"]:
            fa_combo_df.values.sort(kind="mergesort")  # safe sort by numpy
            print(
                "[INFO] --> Number of predicted lipids (exact position): ",
                fa_combo_df.shape[0],
            )
            fa_combo_df["DISCRETE_ABBR"] = (
                fa_combo_df["CLASS"]
                + "("
                + fa_combo_df["FA1"].str.strip("FA")
                + "_"
                + fa_combo_df["FA2"].str.strip("FA")
                + "_"
                + fa_combo_df["FA3"].str.strip("FA")
                + "_"
                + fa_combo_df["FA4"].str.strip("FA")
                + ")"
            )
            fa_combo_df.sort_values(by="DISCRETE_ABBR", inplace=True)
            print(
                "[INFO] --> Number of predicted lipids (exact position): ",
                fa_combo_df.shape[0],
            )
        else:
            fa_combo_df["DISCRETE_ABBR"] = ""
            print("[WARNING] !!! Number of predicted lipids (exact position): 0")

        if position is False:
            print("[INFO] --> Use discrete form for identification ...")
            fa_combo_lite_df = fa_combo_df.drop_duplicates(
                subset=["DISCRETE_ABBR"], keep="first"
            ).copy()
            print(
                "[INFO] --> Number of predicted lipids (discrete form): ",
                fa_combo_lite_df.shape[0],
            )
        else:
            fa_combo_lite_df = fa_combo_df.copy()

        # fa_combo_lite_df.is_copy = False
        fa_combo_lite_df["idx"] = fa_combo_lite_df["DISCRETE_ABBR"]
        fa_combo_lite_df.set_index("idx", drop=True, inplace=True)

        lipid_comb_dct = fa_combo_lite_df.to_dict(orient="index")

        return lipid_comb_dct

    @staticmethod
    def calc_fragments(lipid_dct, charge="", ms2_ppm=100):

        # m_formula = lipid_dct['FORMULA']
        m_exactmass = lipid_dct["EXACTMASS"]
        m_class = lipid_dct["CLASS"]

        h_exactmass = 1.0078250321
        na_exactmass = 22.98976967
        nh3_exactmass = 3 * 1.0078250321 + 14.0030740052
        ch3_exactmass = 12.0 + 3 * 1.0078250321
        nl_water = 2 * 1.0078250321 + 15.9949146221
        gly_mg_base_exactmass = 3 * 12.0 + 5 * 1.0078250321 + 15.9949146221
        fa1_abbr = lipid_dct["FA1"].strip("FA")
        fa2_abbr = lipid_dct["FA2"].strip("FA")

        fa1_exactmass = lipid_dct["FA1_EXACTMASS"]
        if "FA2_EXACTMASS" in list(lipid_dct.keys()):
            fa2_exactmass = lipid_dct["FA2_EXACTMASS"]

        if m_class in ["PA", "PE", "PG", "PI", "PS", "PIP", "PL"]:

            lyso_str = "L" + m_class
            # create the abbreviation name for the Lyso fragments eg. LPE(18:0)-H]-_ABBR
            # without the loss of water
            lipid_dct["[LPL(FA1)-H]-_ABBR"] = "[%s(%s)-H]-" % (lyso_str, fa1_abbr)
            lipid_dct["[LPL(FA2)-H]-_ABBR"] = "[%s(%s)-H]-" % (lyso_str, fa2_abbr)
            # with the loss of water
            lipid_dct["[LPL(FA1)-H2O-H]-_ABBR"] = "[%s(%s)-H2O-H]-" % (
                lyso_str,
                fa1_abbr,
            )
            lipid_dct["[LPL(FA2)-H2O-H]-_ABBR"] = "[%s(%s)-H2O-H]-" % (
                lyso_str,
                fa2_abbr,
            )

            # calculation of the exact mass for the different lyso fragments
            lipid_dct["[LPL(FA1)-H]-_MZ"] = round(
                m_exactmass - (fa2_exactmass - nl_water) - h_exactmass, 6
            )
            lipid_dct["[LPL(FA2)-H]-_MZ"] = round(
                m_exactmass - (fa1_exactmass - nl_water) - h_exactmass, 6
            )
            lipid_dct["[LPL(FA1)-H2O-H]-_MZ"] = round(
                m_exactmass - fa2_exactmass - h_exactmass, 6
            )
            lipid_dct["[LPL(FA2)-H2O-H]-_MZ"] = round(
                m_exactmass - fa1_exactmass - h_exactmass, 6
            )

        # Todo: Add SM section here and check
        elif m_class in ["PC", "SM"]:

            lyso_str = "L" + m_class
            # The abbr. here is not exactly correct due to the compatibility issues with ranks core calc functions
            lipid_dct["[LPL(FA1)-H]-_ABBR"] = "[%s(%s)-CH3]-" % (lyso_str, fa1_abbr)
            lipid_dct["[LPL(FA2)-H]-_ABBR"] = "[%s(%s)-CH3]-" % (lyso_str, fa2_abbr)
            lipid_dct["[LPL(FA1)-H2O-H]-_ABBR"] = "[%s(%s)-H2O-CH3]-" % (
                lyso_str,
                fa1_abbr,
            )
            lipid_dct["[LPL(FA2)-H2O-H]-_ABBR"] = "[%s(%s)-H2O-CH3]-" % (
                lyso_str,
                fa2_abbr,
            )

            lipid_dct["[LPL(FA1)-H]-_MZ"] = round(
                m_exactmass - (fa2_exactmass - nl_water) - ch3_exactmass, 6
            )
            lipid_dct["[LPL(FA2)-H]-_MZ"] = round(
                m_exactmass - (fa1_exactmass - nl_water) - ch3_exactmass, 6
            )
            lipid_dct["[LPL(FA1)-H2O-H]-_MZ"] = round(
                m_exactmass - fa2_exactmass - ch3_exactmass, 6
            )
            lipid_dct["[LPL(FA2)-H2O-H]-_MZ"] = round(
                m_exactmass - fa1_exactmass - ch3_exactmass, 6
            )

        elif m_class in ["LPA", "LPE", "LPG", "LPI", "LPS", "LPIP", "LPL"]:
            pass

        elif m_class in ["LPC"]:
            pass

        elif m_class in ["TG"]:
            # Here maybe should get the DG fragments
            # TODO(georgia.angelidou@uni-leipzig.de): create the section for theuniue fragments when there is TG
            #   Missing the fragments for the sodium adduct

            # The different frgments for triacylglycerol names when neutral loss of the FA
            # Take the correspond information of the 3 FA
            fa3_abbr = lipid_dct["FA3"].strip("FA")
            fa3_exactmass = lipid_dct["FA3_EXACTMASS"]
            dg_str = "M"
            if charge in ["[M+Na]+"]:
                fa1_Na_exactmass = lipid_dct["FA1_[FA-H+Na]_MZ"]
                fa2_Na_exactmass = lipid_dct["FA2_[FA-H+Na]_MZ"]
                fa3_Na_exactmass = lipid_dct["FA3_[FA-H+Na]_MZ"]

                lipid_dct["[M-(FA1)+Na]+_ABBR"] = "[%s-FA%s+Na]+" % (dg_str, fa1_abbr)
                lipid_dct["[M-(FA2)+Na]+_ABBR"] = "[%s-FA%s+Na]+" % (dg_str, fa2_abbr)
                lipid_dct["[M-(FA3)+Na]+_ABBR"] = "[%s-FA%s+Na]+" % (dg_str, fa3_abbr)

                lipid_dct["[M-(FA1-H+Na)+H]+_ABBR"] = "[%s-(FA%s-H+Na)+H]+" % (
                    dg_str,
                    fa1_abbr,
                )
                lipid_dct["[M-(FA2-H+Na)+H]+_ABBR"] = "[%s-(FA%s-H+Na)+H]+" % (
                    dg_str,
                    fa2_abbr,
                )
                lipid_dct["[M-(FA3-H+Na)+H]+_ABBR"] = "[%s-(FA%s-H+Na)+H]+" % (
                    dg_str,
                    fa3_abbr,
                )

                lipid_dct["[M-(FA1)+Na]+_MZ"] = round(
                    m_exactmass - fa1_exactmass + na_exactmass, 6
                )
                lipid_dct["[M-(FA2)+Na]+_MZ"] = round(
                    m_exactmass - fa2_exactmass + na_exactmass, 6
                )
                lipid_dct["[M-(FA3)+Na]+_MZ"] = round(
                    m_exactmass - fa3_exactmass + na_exactmass, 6
                )

                lipid_dct["[M-(FA1-H+Na)+H]+_MZ"] = round(
                    m_exactmass - fa1_Na_exactmass + na_exactmass, 6
                )
                lipid_dct["[M-(FA2-H+Na)+H]+_MZ"] = round(
                    m_exactmass - fa2_Na_exactmass + na_exactmass, 6
                )
                lipid_dct["[M-(FA3-H+Na)+H]+_MZ"] = round(
                    m_exactmass - fa3_Na_exactmass + na_exactmass, 6
                )

                lipid_dct["[M-(FA1)+Na]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - (fa1_exactmass) + na_exactmass), ms2_ppm * -1
                )
                lipid_dct["[M-(FA1)+Na]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - (fa1_exactmass) + na_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA1)+Na]+_Q"] = (
                    lipid_dct["[M-(FA1)+Na]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA1)+Na]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[M-(FA2)+Na]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - (fa2_exactmass) + na_exactmass), ms2_ppm * -1
                )
                lipid_dct["[M-(FA2)+Na]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - (fa2_exactmass) + na_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA2)+Na]+_Q"] = (
                    lipid_dct["[M-(FA2)+Na]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA2)+Na]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[M-(FA3)+Na]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - (fa3_exactmass) + na_exactmass), ms2_ppm * -1
                )
                lipid_dct["[M-(FA3)+Na]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - (fa3_exactmass) + na_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA3)+Na]+_Q"] = (
                    lipid_dct["[M-(FA3)+Na]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA3)+Na]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[M-(FA1-H+Na)+H]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - (fa1_Na_exactmass) + na_exactmass), ms2_ppm * -1
                )
                lipid_dct["[M-(FA1-H+Na)+H]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - (fa1_Na_exactmass) + na_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA1-H+Na)+H]+_Q"] = (
                    lipid_dct["[M-(FA1-H+Na)+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA1-H+Na)+H]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[M-(FA2-H+Na)+H]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - (fa2_Na_exactmass) + na_exactmass), ms2_ppm * -1
                )
                lipid_dct["[M-(FA2-H+Na)+H]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - (fa2_Na_exactmass) + na_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA2-H+Na)+H]+_Q"] = (
                    lipid_dct["[M-(FA2-H+Na)+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA2-H+Na)+H]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[M-(FA3-H+Na)+H]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - (fa3_Na_exactmass) + na_exactmass), ms2_ppm * -1
                )
                lipid_dct["[M-(FA3-H+Na)+H]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - (fa3_Na_exactmass) + na_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA3-H+Na)+H]+_Q"] = (
                    lipid_dct["[M-(FA3-H+Na)+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA3-H+Na)+H]+_MZ_HIGH"].astype(str)
                )
            else:
                # Neutral loss of a FA with a water
                lipid_dct["[M-(FA1)+H]+_ABBR"] = "[%s-FA%s+H]+" % (dg_str, fa1_abbr)
                lipid_dct["[M-(FA2)+H]+_ABBR"] = "[%s-FA%s+H]+" % (dg_str, fa2_abbr)
                lipid_dct["[M-(FA3)+H]+_ABBR"] = "[%s-FA%s+H]+" % (dg_str, fa3_abbr)
                # Neutral loss of a FA minus a water
                lipid_dct["[M-(FA1-H2O)+H]+_ABBR"] = "[%s-(FA%s-H2O)+H]+" % (
                    dg_str,
                    fa1_abbr,
                )
                lipid_dct["[M-(FA2-H2O)+H]+_ABBR"] = "[%s-(FA%s-H2O)+H]+" % (
                    dg_str,
                    fa2_abbr,
                )
                lipid_dct["[M-(FA3-H2O)+H]+_ABBR"] = "[%s-(FA%s-H2O)+H]+" % (
                    dg_str,
                    fa3_abbr,
                )

                lipid_dct["[M-(FA1)+H]+_MZ"] = round(
                    m_exactmass - fa1_exactmass + h_exactmass, 6
                )
                lipid_dct["[M-(FA2)+H]+_MZ"] = round(
                    m_exactmass - fa2_exactmass + h_exactmass, 6
                )
                lipid_dct["[M-(FA3)+H]+_MZ"] = round(
                    m_exactmass - fa3_exactmass + h_exactmass, 6
                )

                lipid_dct["[M-(FA1-H2O)+H]+_MZ"] = round(
                    m_exactmass - (fa1_exactmass - nl_water) + h_exactmass, 6
                )
                lipid_dct["[M-(FA2-H2O)+H]+_MZ"] = round(
                    m_exactmass - (fa2_exactmass - nl_water) + h_exactmass, 6
                )
                lipid_dct["[M-(FA3-H2O)+H]+_MZ"] = round(
                    m_exactmass - (fa3_exactmass - nl_water) + h_exactmass, 6
                )

                lipid_dct["[MG(FA1)-H2O+H]+_MZ_LOW"] = ppm_window_para(
                    (fa1_exactmass + gly_mg_base_exactmass), ms2_ppm * -1
                )
                lipid_dct["[MG(FA1)-H2O+H]+_MZ_HIGH"] = ppm_window_para(
                    (fa1_exactmass + gly_mg_base_exactmass), ms2_ppm
                )
                lipid_dct["[MG(FA1)-H2O+H]+_Q"] = (
                    lipid_dct["[MG(FA1)-H2O+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[MG(FA1)-H2O+H]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[MG(FA2)-H2O+H]+_MZ_LOW"] = ppm_window_para(
                    (fa2_exactmass + gly_mg_base_exactmass), ms2_ppm * -1
                )
                lipid_dct["[MG(FA2)-H2O+H]+_MZ_HIGH"] = ppm_window_para(
                    (fa2_exactmass + gly_mg_base_exactmass), ms2_ppm
                )
                lipid_dct["[MG(FA2)-H2O+H]+_Q"] = (
                    lipid_dct["[MG(FA2)-H2O+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[MG(FA2)-H2O+H]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[MG(FA3)-H2O+H]+_MZ_LOW"] = ppm_window_para(
                    (fa3_exactmass + gly_mg_base_exactmass), ms2_ppm * -1
                )
                lipid_dct["[MG(FA3)-H2O+H]+_MZ_HIGH"] = ppm_window_para(
                    (fa3_exactmass + gly_mg_base_exactmass), ms2_ppm
                )
                lipid_dct["[MG(FA3)-H2O+H]+_Q"] = (
                    lipid_dct["[MG(FA3)-H2O+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[MG(FA3)-H2O+H]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[M-(FA1)+H]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - fa1_exactmass + h_exactmass), ms2_ppm * -1
                )
                lipid_dct["[M-(FA1)+H]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - fa1_exactmass + h_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA1)+H]+_Q"] = (
                    lipid_dct["[M-(FA1)+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA1)+H]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[M-(FA2)+H]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - fa2_exactmass + h_exactmass), ms2_ppm * -1
                )
                lipid_dct["[M-(FA2)+H]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - fa2_exactmass + h_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA2)+H]+_Q"] = (
                    lipid_dct["[M-(FA2)+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA2)+H]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[M-(FA3)+H]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - fa3_exactmass + h_exactmass), ms2_ppm * -1
                )
                lipid_dct["[M-(FA3)+H]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - fa3_exactmass + h_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA3)+H]+_Q"] = (
                    lipid_dct["[M-(FA3)+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA3)+H]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[M-(FA1-H2O)+H]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - (fa1_exactmass - nl_water) + h_exactmass),
                    ms2_ppm * -1,
                )
                lipid_dct["[M-(FA1-H2O)+H]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - (fa1_exactmass - nl_water) + h_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA1-H2O)+H]+_Q"] = (
                    lipid_dct["[M-(FA1-H2O)+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA1-H2O)+H]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[M-(FA2-H2O)+H]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - (fa2_exactmass - nl_water) + h_exactmass),
                    ms2_ppm * -1,
                )
                lipid_dct["[M-(FA2-H2O)+H]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - (fa2_exactmass - nl_water) + h_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA2-H2O)+H]+_Q"] = (
                    lipid_dct["[M-(FA2-H2O)+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA2-H2O)+H]+_MZ_HIGH"].astype(str)
                )

                lipid_dct["[M-(FA3-H2O)+H]+_MZ_LOW"] = ppm_window_para(
                    (m_exactmass - (fa3_exactmass - nl_water) + h_exactmass),
                    ms2_ppm * -1,
                )
                lipid_dct["[M-(FA3-H2O)+H]+_MZ_HIGH"] = ppm_window_para(
                    (m_exactmass - (fa3_exactmass - nl_water) + h_exactmass), ms2_ppm
                )
                lipid_dct["[M-(FA3-H2O)+H]+_Q"] = (
                    lipid_dct["[M-(FA3-H2O)+H]+_MZ_LOW"].astype(str)
                    + " <= mz <= "
                    + lipid_dct["[M-(FA3-H2O)+H]+_MZ_HIGH"].astype(str)
                )
            # Fragments names when can occur 2 neutral losses of FA. 1 FA with the water and other without
            mg_str = "MG"
            lipid_dct["[MG(FA1)-H2O+H]+_ABBR"] = "[%s(%s)-H2O+H]+" % (mg_str, fa1_abbr)
            lipid_dct["[MG(FA2)-H2O+H]+_ABBR"] = "[%s(%s)-H2O+H]+" % (mg_str, fa2_abbr)
            lipid_dct["[MG(FA3)-H2O+H]+_ABBR"] = "[%s(%s)-H2O+H]+" % (mg_str, fa3_abbr)

            lipid_dct["[MG(FA1)-H2O+H]+_MZ"] = round(
                fa1_exactmass + gly_mg_base_exactmass, 6
            )
            lipid_dct["[MG(FA2)-H2O+H]+_MZ"] = round(
                fa2_exactmass + gly_mg_base_exactmass, 6
            )
            lipid_dct["[MG(FA3)-H2O+H]+_MZ"] = round(
                fa3_exactmass + gly_mg_base_exactmass, 6
            )

            # TODO (georgia.angelidou@uni-leipzig.de): add the fragments of Na [sodium]
            # Fragments when there only 1 neutral loss from the TG.

        elif m_class in ["DG"]:
            mg_str = "MG"
            lipid_dct["[MG(FA1)-H2O+H]+_ABBR"] = "[%s(%s)-H2O+H]+" % (mg_str, fa1_abbr)
            lipid_dct["[MG(FA2)-H2O+H]+_ABBR"] = "[%s(%s)-H2O+H]+" % (mg_str, fa2_abbr)

            lipid_dct["[MG(FA1)-H2O+H]+_MZ"] = round(
                fa1_exactmass + gly_mg_base_exactmass, 6
            )
            lipid_dct["[MG(FA2)-H2O+H]+_MZ"] = round(
                fa2_exactmass + gly_mg_base_exactmass, 6
            )

            lipid_dct["[MG(FA1)-H2O+H]+_MZ_LOW"] = ppm_window_para(
                (fa1_exactmass + gly_mg_base_exactmass), ms2_ppm * -1
            )
            lipid_dct["[MG(FA1)-H2O+H]+_MZ_HIGH"] = ppm_window_para(
                (fa1_exactmass + gly_mg_base_exactmass), ms2_ppm
            )
            lipid_dct["[MG(FA1)-H2O+H]+_Q"] = (
                lipid_dct["[MG(FA1)-H2O+H]+_MZ_LOW"].astype(str)
                + " <= mz <= "
                + lipid_dct["[MG(FA1)-H2O+H]+_MZ_HIGH"].astype(str)
            )

            lipid_dct["[MG(FA2)-H2O+H]+_MZ_LOW"] = ppm_window_para(
                (fa2_exactmass + gly_mg_base_exactmass), ms2_ppm * -1
            )
            lipid_dct["[MG(FA2)-H2O+H]+_MZ_HIGH"] = ppm_window_para(
                (fa2_exactmass + gly_mg_base_exactmass), ms2_ppm
            )
            lipid_dct["[MG(FA2)-H2O+H]+_Q"] = (
                lipid_dct["[MG(FA2)-H2O+H]+_MZ_LOW"].astype(str)
                + " <= mz <= "
                + lipid_dct["[MG(FA2)-H2O+H]+_MZ_HIGH"].astype(str)
            )

        else:
            # TODO (georgia.angelidou@uni-leipzig.de: Info for sphingomyelins
            pass

        return lipid_dct

    def compose_lipid(self, param_dct: dict, ms2_ppm: int = 100) -> pd.DataFrame:
        """
        Main function to generate the combination of lipids.
        Control all functions above.
        Args:
            param_dct (dict): parameter load from GUI or config file
            ms2_ppm (int): The MS/MS tolerance for the calculation of MS/MS query window

        Returns:

        """

        lipid_class = param_dct["lipid_class"]
        lipid_charge = param_dct["charge_mode"]
        if param_dct["exact_position"] == "TRUE":
            position_set = True
        else:
            position_set = False

        usr_fa_df = pd.read_excel(param_dct["fa_whitelist"])
        usr_fa_df = usr_fa_df.fillna(value="F")

        tmp_columns = usr_fa_df.columns.tolist()

        usr_fa_df.columns = usr_fa_df.columns.str.upper()

        # Todo: check SM section below:
        if lipid_class in ["PL", "PA", "PC", "PE", "PG", "PI", "PS", "SM"]:
            if lipid_class in tmp_columns:
                pass
            elif "PL" in tmp_columns:
                pass
            else:
                return False
        elif lipid_class in ["LPL", "LPA", "LPC", "LPE", "LPG", "LPI", "LPS"]:
            if lipid_class in tmp_columns:
                pass
            elif "LPL" in tmp_columns or "PL" in tmp_columns:
                pass
            else:
                return False

        elif lipid_class in ["TG", "DG"]:
            if lipid_class in tmp_columns:
                pass
            else:
                return False
        else:
            return False
        print("[INFO] --> FA white list loaded ...")
        lipid_comb_dct = self.gen_all_comb(
            lipid_class, usr_fa_df, position=position_set
        )

        lipid_info_dct = {}

        abbr_parser = NameParserFA()
        elem_calc = ElemCalc()
        for _lipid in list(lipid_comb_dct.keys()):
            _lipid_dct = lipid_comb_dct[_lipid]

            _fa1_abbr = _lipid_dct["FA1"]
            if _fa1_abbr:
                _fa1_info_dct = abbr_parser.get_fa_info(_fa1_abbr)
                for _fa1_k in list(_fa1_info_dct.keys()):
                    _lipid_dct["FA1_" + _fa1_k] = _fa1_info_dct[_fa1_k]
            else:
                if "FA2" in list(_lipid_dct.keys()):
                    if _lipid_dct["FA2"]:
                        _lipid_dct["FA1"], _lipid_dct["FA2"] = (
                            _lipid_dct["FA2"],
                            _lipid_dct["FA1"],
                        )
                        _fa1_abbr = _lipid_dct["FA1"]
                        _fa1_info_dct = abbr_parser.get_fa_info(_fa1_abbr)
                        for _fa1_k in list(_fa1_info_dct.keys()):
                            _lipid_dct["FA1_" + _fa1_k] = _fa1_info_dct[_fa1_k]
                else:
                    _fa1_info_dct = {}

            if "FA2" in list(_lipid_dct.keys()):
                _fa2_abbr = _lipid_dct["FA2"]
                if _fa2_abbr:
                    _fa2_info_dct = abbr_parser.get_fa_info(_fa2_abbr)
                    for _fa2_k in list(_fa2_info_dct.keys()):
                        _lipid_dct["FA2_" + _fa2_k] = _fa2_info_dct[_fa2_k]
            else:
                _fa2_abbr = ""
                _fa2_info_dct = {}

            # Todo: add and check SM, Cer below
            if lipid_class in ["PA", "PC", "PE", "PG", "PI", "PS", "DG"]:
                _lipid_dct["M_DB"] = _fa1_info_dct["DB"] + _fa2_info_dct["DB"]
                # TODO(georgia.angelidou@uni-leipzig.de): not important (just keep in mind for future correction)
                # consideration the case if the user choose the fa2 position for the different link types
                # or in that they may mistype something
                if _fa1_info_dct["LINK"] in ["FA", "A"]:
                    lipid_bulk_str = "{pl}({c}:{db})".format(
                        pl=lipid_class,
                        c=_fa1_info_dct["C"] + _fa2_info_dct["C"],
                        db=lipid_comb_dct[_lipid]["M_DB"],
                    )
                else:
                    lipid_bulk_str = "{pl}({lk}{c}:{db})".format(
                        pl=lipid_class,
                        lk=_fa1_info_dct["LINK"],
                        c=_fa1_info_dct["C"] + _fa2_info_dct["C"],
                        db=lipid_comb_dct[_lipid]["M_DB"],
                    )
            elif lipid_class in ["LPA", "LPC", "LPE", "LPG", "LPI", "LPS"]:
                _lipid_dct["M_DB"] = _fa1_info_dct["DB"]
                # TODO(georgia.angelidou@uni-leipzig.de): not important (just keep in mind for future correction)
                # consideration the case if the user choose the fa2 position for the different link types
                # or in that they may mistype something
                if _fa1_info_dct["LINK"] in ["FA", "A"]:
                    lipid_bulk_str = "{pl}({c}:{db})".format(
                        pl=lipid_class,
                        c=_fa1_info_dct["C"],
                        db=lipid_comb_dct[_lipid]["M_DB"],
                    )
                else:
                    lipid_bulk_str = "{pl}({lk}{c}:{db})".format(
                        pl=lipid_class,
                        lk=_fa1_info_dct["LINK"],
                        c=_fa1_info_dct["C"],
                        db=lipid_comb_dct[_lipid]["M_DB"],
                    )
            elif lipid_class in ["TG"]:
                _fa3_abbr = _lipid_dct["FA3"]
                _fa3_info_dct = abbr_parser.get_fa_info(_fa3_abbr)

                for _fa3_k in _fa3_info_dct.keys():
                    _lipid_dct["FA3_" + _fa3_k] = _fa3_info_dct[_fa3_k]

                _lipid_dct["M_DB"] = (
                    _fa3_info_dct["DB"] + _fa2_info_dct["DB"] + _fa1_info_dct["DB"]
                )
                # Note: For TG in the current default not consider the different lipids with other type of bond
                # If stay like this need to be mention in somewhere for the user
                if _fa1_info_dct["LINK"] in ["FA", "A"]:
                    lipid_bulk_str = "{tg}({c}:{db})".format(
                        tg=lipid_class,
                        c=(
                            _fa1_info_dct["C"] + _fa2_info_dct["C"] + _fa3_info_dct["C"]
                        ),
                        db=lipid_comb_dct[_lipid]["M_DB"],
                    )
                else:
                    lipid_bulk_str = "{tg}({lk}{c}:{db})".format(
                        tg=lipid_class,
                        lk=_fa1_info_dct["LINK"],
                        c=(
                            _fa1_info_dct["C"] + _fa2_info_dct["C"] + _fa3_info_dct["C"]
                        ),
                        db=lipid_comb_dct[_lipid]["M_DB"],
                    )
            # Todo: check SM section below:
            elif lipid_class in ["SM"]:
                lipid_bulk_str = "{sm}({c}:{db})".format(
                    sm=lipid_class,
                    c=_fa1_info_dct["C"] + _fa2_info_dct["C"],
                    db=lipid_comb_dct[_lipid]["M_DB"],
                )
            else:
                lipid_bulk_str = ""

            _lipid_dct["BULK_ABBR"] = lipid_bulk_str

            _lipid_formula, _lipid_elem_dct = elem_calc.get_formula(lipid_bulk_str)

            _lipid_dct["FORMULA"] = _lipid_formula
            _lipid_dct["EXACTMASS"] = elem_calc.get_exactmass(_lipid_elem_dct)
            for _elem_k in list(_lipid_elem_dct.keys()):
                _lipid_dct["M_" + _elem_k] = _lipid_elem_dct[_elem_k]

            # charged
            _chg_lipid_formula, _chg_lipid_elem_dct = elem_calc.get_formula(
                lipid_bulk_str, charge=lipid_charge
            )
            _lipid_dct[lipid_charge + "_FORMULA"] = _chg_lipid_formula
            _lipid_dct[lipid_charge + "_MZ"] = elem_calc.get_exactmass(
                _chg_lipid_elem_dct
            )

            # fragments

            _lipid_dct = self.calc_fragments(
                _lipid_dct, charge=lipid_charge, ms2_ppm=ms2_ppm
            )

            lipid_info_dct[_lipid] = _lipid_dct
            del _lipid_dct

        lipid_master_df = pd.DataFrame(lipid_comb_dct).T
        lipid_master_df.reset_index(drop=True, inplace=True)

        return lipid_master_df
