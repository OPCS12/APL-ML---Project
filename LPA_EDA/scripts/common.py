"""Persistência e gráficos determinísticos; não modifica os dados originais."""
from pathlib import Path
import hashlib, json
import os
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'logs/cache/matplotlib'))
os.environ.setdefault('XDG_CACHE_HOME',str(Path(__file__).resolve().parents[1]/'logs/cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
SEED = 172057

def table(df, name):
    df.to_csv(ROOT/'tabelas'/f'{name}.csv', index=False)

def figure(name):
    plt.tight_layout()
    plt.savefig(ROOT/'figuras'/f'{name}.png', dpi=160, bbox_inches='tight')
    plt.close('all')

def checksum(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, default=str))
