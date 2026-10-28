"""Executa o notebook existente com kernel Python e caches graváveis."""
import os
from pathlib import Path
import nbformat
from nbclient import NotebookClient
ROOT=Path(__file__).resolve().parents[1]
for key,sub in [('MPLCONFIGDIR','matplotlib'),('XDG_CACHE_HOME',''),('IPYTHONDIR','ipython'),('JUPYTER_RUNTIME_DIR','jupyter')]:
    target=ROOT/'logs/cache'/sub;target.mkdir(parents=True,exist_ok=True);os.environ.setdefault(key,str(target))
p=ROOT/'notebooks/LPA_EDA.ipynb'
nb=nbformat.read(p,as_version=4)
NotebookClient(nb,timeout=120,resources={'metadata':{'path':str(ROOT)}},kernel_name='python3').execute()
nbformat.write(nb,p)
print('Notebook executado:',sum(c.cell_type=='code' for c in nb.cells),'células de código sem erro')
