set PATH=%PATH%;%HOMEPATH%\Miniconda2;%HOMEPATH%\Miniconda2\Scripts;C:\Miniconda2\Scripts;C:\Miniconda2
set APPNAME=teigen
conda create -y -c mjirik -c SimpleITK -c menpo --no-default-packages -n %APPNAME% pywget wget pip %APPNAME%

call activate %APPNAME%
rem application acutalization if the installer is running for second time (and prev line fails because of existance of teigen environment)
conda install -y -c mjirik %APPNAME%
rem :windows specific
rem conda install -c jmargeta scikit-fmm

rem python -m wget https://raw.githubusercontent.com/mjirik/lisa/master/requirements_pip.txt
rem pip install -r requirements_pip.txt
rem del requirements_pip.txt
