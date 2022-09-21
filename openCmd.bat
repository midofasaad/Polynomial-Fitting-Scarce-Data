@echo off
color 0A
title Filter Evaluation 
if EXIST .\env38 GOTO CHECKDBPW
echo creating virtual environment: env38
python -m venv env38
echo setting intranet pip repository (codecraft)
pip config set global.index-url https://artifactory.cc.bmwgroup.net/artifactory/api/pypi/external-pypi-org/simple
echo activating virtual environment and updating pip in env38
cmd /C ".\env38\Scripts\activate & python -m pip install --upgrade pip"
echo installing required modules from requirements.txt into env38
cmd /C ".\env38\Scripts\activate & python -m pip install -r requirements.txt" 
:CHECKDBPW
echo activating virtual environment: env38
cmd /K ".\env38\Scripts\activate & echo start test server with python app.py"
exit