# LipidHunter#

LipidHunter 2 has significant improvement from the original LipidHunter.
The major improvements are listed below:

* Special Feature for TG & DG
    + ID based both on FA neutral losses & fragments
    + Accurate ID of isomeric species
    + Correction for all identified FA
    + Correction for the fragment intensity of multiply identified FAs
    + Check for cross-contamination signals e.g. PL fragments
* Other Major Feature Updates
    + Multiprocessing
    + Batch mode
    + 10 times faster processing speed
    + Command line mode
    + KNIME workflow integration
    + Multiple vendor support
    + Improved output style
    + Simplified configuration
    + View run parameters in report
    
![crossplatform_screenshot.png](doc/img/Hunter2_GUI.png)

This repository contains the source code of LipidHunter.

LipidHunter Windows .exe executable version can be found in release page:

https://github.com/SysMedOs/lipidhunter/releases


** Please read the following instructions before you start to run LipidHunter. **

## Instructions ##

### How to install LipidHunter from source code ###
* Download the LipidHunter as zip file for your system

    + Download LipidHunter source Code as .zip. Please notice the date and version of LipidHunter source code.
    + Professional users can use `git` to clone the whole repository, please make sure that you switched to the correct branch.

* Rename the downloaded file to `LipidHunter.zip`
* Unzip `LipidHunter.zip` file to any folder.
* Downloaded LipidHunter test spectra files: [LipidHunter_Test_mzML_File](https://bitbucket.org/SysMedOs/lipidhunter_exe/downloads/LipidHunter_TestFile.zip)

* Python environment

    + LipidHunter 2 is developed under python 3.6, the current version can still run on python 2.7 (not recommended).
    + The best way is to use virtual environment such as `conda`
    + The requirements is listed in [requirements.txt](requirements.txt)

* [LipidHunter user guide](doc/LipidHunter_UserGuide.pdf)


* Errors/bugs
    
    In case you experienced any problems with running LipidHunter, please report an issue in the [issue tracker](https://github.com/SysMedOs/lipidhunter/issues) or contact us.

### License ###

+ LipidHunter is Dual-licensed
    * For academic and non-commercial use: `GPLv2 License`: 
    
        [The GNU General Public License version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

    * For commercial use: please contact the develop team by email.

+ Please cite our publication in an appropriate form. 

### Further questions? ###

* Report any issues here: [https://github.com/SysMedOs/lipidhunter/issues](https://github.com/SysMedOs/lipidhunter/issues)


### Fundings ###
We acknowledge all projects that supports the development of LipidHunter:

+ BMBF - Federal Ministry of Education and Research Germany:

    https://www.bmbf.de/en/

+ e:Med Systems Medicine Network:

    http://www.sys-med.de/en/

+ SysMedOS Project : 

    https://home.uni-leipzig.de/fedorova/sysmedos/
