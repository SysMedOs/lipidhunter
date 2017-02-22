# LipidHunter Instructions #

![LipidHunter_logo_128.png](https://bitbucket.org/repo/oGzkj4/images/583628216-LipidHunter_logo_128.png)

This repository contains the source code of LipidHunter.

There is another repository about the Windows .exe executable version of LipidHunter at 
https://bitbucket.org/SysMedOs/lipidhunterdist

Here we will provide general information how to download the Windows executable version.
 
The LipidHunter Windows executable version is provided for `64bit` version of Windows 7, 8, 8.1 and 10 only.

For more experienced users, or if you have another operating system (`Linux`, `macOS`), there are also instructions about downloading and running the source code.

LipidHunter requires minimum 2.0 GHz CPU and 4 GB RAM to run, for multiple processes of LipidHunter, 8 GB or 16 GB of RAM are recommended.

Though we tried our best to keep the Windows executable distribution up to date. We have no guarantee that the Windows executable distribution corresponding to the latest version of LipidHunter source code. Please check out the original LipidHunter repository for latest features and bug fix.

If you have any problems while using LipidHunter, please report your issue here: [https://bitbucket.org/SysMedOs/lipidhunter/issues](https://bitbucket.org/SysMedOs/lipidhunter/issues)

** Please read the following instructions before you start to run LipidHunter. **

### Instructions ###

* [How to install Windows executable version of LipidHunter](#markdown-header-how-to-install-windows-executable-version-of-lipidhunter)
* [How to install the Source Code](#markdown-header-how-to-install-lipidhunter-from-source-code)
* [License](#markdown-header-license)
* [A step by step tutorial](https://bitbucket.org/SysMedOs/lipidhunter/wiki/Home)
* [Q&A](#markdown-header-further-questions)
* [Fundings](#markdown-header-fundings)



### How to install Windows executable version of LipidHunter ###
* Download the LipidHunter for your system:

    + [LipidHunter for Windows 10 64bit ( `.zip` file around 320 MB)](https://bitbucket.org/SysMedOs/lipidhunterdist/downloads/LipidHunter_Win10_64bit.zip)
    + [LipidHunter for Windows 7,8 and 8.1 64bit ( `.zip` file around 230 MB)](https://bitbucket.org/SysMedOs/lipidhunterdist/downloads/LipidHunter_Win7-8_64bit.zip)
    
* Rename the downloaded file  to `LipidHunter.zip`
    
* Unzip `LipidHunter.zip` file to any folder. We recommend to use `7-zip` to unzip this file.

    + 7-Zip is open source software. [7-Zip home page](http://www.7-zip.org)
            
* **Optional:** Resolve `.dll` issues. LipidHunter is a free open source software. However, some system files may required to execute the LipidHunter.exe. If you receive .dll error from LipidHunter, try to copy following system files to the LipidHunter folder and try again.
    + Copy following `.dll` and `.DRV` files from your system folder: `C:\WINDOWS\system32\` to the LipidHunter folder:

        + `.DRV` files:
        
                WINSPOOL.DRV - C:\WINDOWS\system32\WINSPOOL.DRV
            
        + `.dll` files:
    
                OLEAUT32.dll - C:\WINDOWS\system32\OLEAUT32.dll
                USER32.dll - C:\WINDOWS\system32\USER32.dll
                IMM32.dll - C:\WINDOWS\system32\IMM32.dll
                SHELL32.dll - C:\WINDOWS\system32\SHELL32.dll
                ole32.dll - C:\WINDOWS\system32\ole32.dll
                ODBC32.dll - C:\WINDOWS\system32\ODBC32.dll
                COMCTL32.dll - C:\WINDOWS\system32\COMCTL32.dll
                ADVAPI32.dll - C:\WINDOWS\system32\ADVAPI32.dll
                SHLWAPI.dll - C:\WINDOWS\system32\SHLWAPI.dll
                WS2_32.dll - C:\WINDOWS\system32\WS2_32.dll
                GDI32.dll - C:\WINDOWS\system32\GDI32.dll
                WINMM.dll - C:\WINDOWS\system32\WINMM.dll
                VERSION.dll - C:\WINDOWS\system32\VERSION.dll
                KERNEL32.dll - C:\WINDOWS\system32\KERNEL32.dll
                COMDLG32.dll - C:\WINDOWS\system32\COMDLG32.dll
                OPENGL32.dll - C:\WINDOWS\system32\OPENGL32.dll
            
        + For windows 7, 8, 8.1 additional `.dll` files may required:
        
                msvcrt.dll - C:\WINDOWS\system32\msvcrt.dll
                ntdll.dll - C:\WINDOWS\system32\ntdll.dll

* Start LipidHunter by `LipidHunter.exe`

* Run the test files
        
        Files will be provided later

* Run your data

* Problems
    
    In case, of problems running the program after following the above steps, contact with us.

### How to install LipidHunter from source code ###
* Download the LipidHunter as zip file for your system

    + [LipidHunter source Code](https://bitbucket.org/SysMedOs/lipidhunter/wiki/Install%20the%20Source%20code)

* Rename the downloaded file to `LipidHunter.zip`
    
* Unzip `LipidHunter.zip` file to any folder.

* Follow the link below for more information to run the program.

    + [Lipid Hunter source code User Guide](https://bitbucket.org/SysMedOs/lipidhunter/wiki/Install%20the%20Source%20code)

### License ###

+ LipidHunter is Dual-licensed
    * For academic and non-commercial use: `GPLv2 License` Please read more information by the following link: 
    
        [The GNU General Public License version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

    * For commercial use: please contact the develop team by email.

+ Please cite our publication in an appropriate form. 

### Further questions? ###

* Read our [wiki](https://bitbucket.org/SysMedOs/lipidhunter/wiki/Home)
* Report your issue here: [https://bitbucket.org/SysMedOs/lipidhunter/issues](https://bitbucket.org/SysMedOs/lipidhunter/issues)


### Fundings ###
We acknowlege all projects that supports the development of LipidHunter:

+ BMBF - Federal Ministry of Education and Research Germany:

    https://www.bmbf.de/en/

+ e:Med Systems Medicine Network:

    http://www.sys-med.de/en/

+ SysMedOS Project : 

    https://home.uni-leipzig.de/fedorova/sysmedos/