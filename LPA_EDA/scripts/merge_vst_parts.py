"""Une partes TSV gzip por linhas, preservando bytes e ordem dos registros.

Não converte valores numéricos nem remove genes repetidos. Somente cabeçalhos
repetidos, previamente verificados como idênticos, deixam de ser copiados.
"""
import argparse
from pathlib import Path
import gzip, hashlib, json, os, shutil
from datetime import datetime, timezone
from collections import Counter


def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()


def merge(parts, output, audit_path):
    if output.exists(): raise FileExistsError(f'Não sobrescrever arquivo existente: {output}')
    temporary=output.with_suffix(output.suffix+'.partial')
    if temporary.exists(): raise FileExistsError(temporary)
    expected=hashlib.sha256();id_hash=hashlib.sha256();headers=[];records=[]
    ids=Counter();n=0;cols=None;previous_ended=True
    try:
        with open(temporary,'xb') as compressed:
            with gzip.GzipFile(filename='',mode='wb',fileobj=compressed,mtime=0) as dest:
                for part in parts:
                    count=0;first=last=None;part_body=hashlib.sha256()
                    with gzip.open(part,'rb') as src:
                        header=src.readline()
                        if not header.endswith(b'\n'): raise ValueError(f'Cabeçalho incompleto: {part}')
                        headers.append(header)
                        if cols is None:
                            cols=header.rstrip(b'\r\n').split(b'\t')
                            if len(set(cols))!=len(cols): raise ValueError('Cabeçalhos duplicados')
                            dest.write(header);expected.update(header)
                        elif header!=headers[0]: raise ValueError(f'Cabeçalho diferente: {part}')
                        if not previous_ended: raise ValueError('Fronteira sem fim de linha; união suspensa')
                        for line in src:
                            fields=line.rstrip(b'\r\n').split(b'\t')
                            if len(fields)!=len(cols): raise ValueError(f'{part}, linha {count+2}: largura inconsistente')
                            gene=fields[0].decode('utf-8')
                            if not gene.startswith('ENSG'): raise ValueError(f'ID não esperado: {gene}')
                            ids[gene]+=1;id_hash.update(fields[0]+b'\n')
                            count+=1;n+=1
                            if first is None:first=gene
                            last=gene;previous_ended=line.endswith(b'\n')
                            dest.write(line);expected.update(line);part_body.update(line)
                    records.append({'file':part.name,'compressed_bytes':part.stat().st_size,'compressed_sha256':sha256(part),'data_rows':count,'first_gene':first,'last_gene':last,'data_bytes_sha256':part_body.hexdigest(),'header_identical':True,'header_copied':len(records)==0})
        actual=hashlib.sha256();observed_ids=hashlib.sha256();recount=0
        with gzip.open(temporary,'rb') as src:
            hdr=src.readline();actual.update(hdr)
            if hdr!=headers[0]: raise AssertionError('Cabeçalho de saída divergente')
            for line in src:
                actual.update(line);observed_ids.update(line.split(b'\t',1)[0]+b'\n');recount+=1
        if actual.digest()!=expected.digest() or observed_ids.digest()!=id_hash.digest() or recount!=n:
            raise AssertionError('Saída não corresponde aos bytes concatenados das entradas')
        temporary.rename(output);output.chmod(0o444)
    except Exception:
        if temporary.exists():temporary.unlink()
        raise
    report={'processed_utc':datetime.now(timezone.utc).isoformat(),'method':'partes na ordem dos argumentos; único cabeçalho; bytes das linhas de dados preservados','parts':records,'output':str(output),'genes_rows':n,'unique_gene_ids':len(ids),'duplicated_gene_ids':{k:v for k,v in ids.items() if v>1},'sample_columns':len(cols)-1,'first_sample':cols[1].decode(),'last_sample':cols[-1].decode(),'removed_repeated_headers':len(parts)-1,'uncompressed_sha256':actual.hexdigest(),'ordered_gene_ids_sha256':id_hash.hexdigest(),'compressed_sha256':sha256(output),'byte_preservation_verified':True,'gene_order_verified':True,'gzip_crc_verified':True,'original_unsplit_file_equivalence':'não verificada; arquivo integral anterior não estava acessível'}
    audit_path.write_text(json.dumps(report,ensure_ascii=False,indent=2))
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--parts',type=Path,nargs='+',required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--audit',type=Path,required=True)
    a=p.parse_args();print(json.dumps(merge(a.parts,a.output,a.audit),ensure_ascii=False,indent=2))
