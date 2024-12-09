# CVKLab-internal-modules
Assorted small functions for our data analysis

image_analysis_functions.py
  - Assorted reading, writing, editing image functions
  - Analysis for TIF, IMS, and OMG files

flim_analysis_functions.py
  - Functions for Fluorescence Lifetime Imaging Analysis
  - Needs phasor analysis

MDAnalysis_functions.py
  - Functions to use the MDAnalysis module for analyzing molecular dynamics trajectories
  - Analyzes and displays GRO, XTC files
  - Occasionally requires Jupyter or JupyterLab

standard_fccs_functions.py
  - functions for reading, writing, and analyzing camera-based fluorescence cross-correlation spectroscopy data
    
pendant_droplet_functions.py
  - functions for modeling and analyzing images of pendant droplet tensiometry


# To automatically load in your local Spyder, edit this code:
```python
### edit these lines for each module or computer
module = 'image_analysis_functions'  # edit this for each module
nickname = 'iaf' # edit this for each module
temp_save_fold = r'C:\Users\cvkelly\Downloads/'  # edit this for each computer

### this code should be constant
import requests, sys
github_url = 'https://raw.githubusercontent.com/CVKellyWSU/CVKLab-internal-modules/refs/heads/main/'+module+'.py'
temp_py_file = github_url[(1+github_url.rfind('/')):]
temp_py_text = requests.get(github_url)
temp_py_text = temp_py_text.text
f = open(temp_save_fold+temp_py_file,'w')
f.write(temp_py_text)
f.close()
sys.path.append(temp_save_fold)
exec('import '+module+' as '+nickname)  

