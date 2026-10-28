"""Leitura estrita: recusa cabeçalhos ambíguos em vez de inventar amostras."""
import csv, gzip, re, json
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from common import ROOT, SEED, table, figure, checksum

def read_matrix(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt',encoding='utf-8-sig',newline='') as f:
        sample=f.read(65536)
    # Testa somente separadores explícitos; espaços não são divididos implicitamente.
    candidates={sep:len(next(csv.reader(sample.splitlines(),delimiter=sep))) for sep in ['\t',',',';']}
    sep=max(candidates,key=candidates.get)
    if candidates[sep]<2: raise ValueError('Separador não identificável sem ambiguidade')
    widths=Counter();bad=[]
    with op(path,'rt',encoding='utf-8-sig',newline='') as f:
        reader=csv.reader(f,delimiter=sep);header=next(reader)
        for line,row in enumerate(reader,2):
            widths[len(row)]+=1
            if len(row)!=len(header): bad.append({'line':line,'fields':len(row),'expected':len(header)})
    audit={'file':str(path),'sha256':checksum(path),'delimiter':repr(sep),'header_fields':len(header),'row_widths':dict(widths),'inconsistent_rows':len(bad),'duplicate_headers':len(header)-len(set(header))}
    table(pd.DataFrame(bad,columns=['line','fields','expected']),'vst_linhas_inconsistentes')
    (ROOT/'logs/vst_estrutura.json').write_text(json.dumps(audit,indent=2))
    if bad or len(set(header))!=len(header) or any(not h.strip() for h in header): raise ValueError('Cabeçalho vazio/duplicado ou largura inconsistente; consultar logs/vst_estrutura.json')
    # Detecta o cabeçalho conhecido e qualquer célula com dois IDs APL.
    if any(len(re.findall(r'APL\d+',h))>1 for h in header): raise ValueError('Cabeçalho contém múltiplos IDs APL na mesma célula; não foi dividido')
    raw=pd.read_csv(path,sep=sep,dtype=str,keep_default_na=False,compression='infer')
    table(raw.head(5),'vst_primeiras_linhas');table(raw.tail(5),'vst_ultimas_linhas')
    ids=raw.iloc[:,0]
    gene_rows=ids.str.fullmatch(r'ENSG\d+(?:\.\d+)?').all()
    gene_cols=all(re.fullmatch(r'ENSG\d+(?:\.\d+)?',h) for h in header[1:])
    if gene_rows==gene_cols: raise ValueError('Orientação não identificada inequivocamente')
    matrix=raw.iloc[:,1:].apply(pd.to_numeric,errors='coerce');matrix.index=ids
    invalid=matrix.isna() | ~np.isfinite(matrix)
    if invalid.any().any():
        ri,ci=np.where(invalid.to_numpy());table(pd.DataFrame({'row_id':matrix.index[ri],'column':matrix.columns[ci],'raw_value':[raw.iloc[r,c+1] for r,c in zip(ri,ci)]}),'vst_valores_invalidos')
        raise ValueError(f'{invalid.to_numpy().sum()} valores ausentes/não numéricos/infinitos; nenhuma imputação aplicada')
    if gene_cols: matrix=matrix.T
    if matrix.index.duplicated().any() or matrix.columns.duplicated().any(): raise ValueError('IDs de genes/amostras duplicados: nenhuma agregação automática')
    if not all(re.fullmatch(r'APL\d+',str(c)) for c in matrix.columns): raise ValueError('IDs de amostra fora do padrão confirmado APLn; revisão necessária')
    audit.update(genes=len(matrix),samples=matrix.shape[1],orientation='genes_rows' if gene_rows else 'genes_columns',missing=0)
    return matrix,audit

def run(path,clinical,s5,second=None):
    m,audit=read_matrix(path)
    if second:
        b,_=read_matrix(second)
        same_ids=m.index.equals(b.index) and m.columns.equals(b.columns)
        audit['second_exact_ids_order_and_values']=bool(same_ids and np.array_equal(m.to_numpy(),b.to_numpy()))
        audit['second_same_id_sets']=bool(set(m.index)==set(b.index) and set(m.columns)==set(b.columns))
        if audit['second_same_id_sets']: audit['second_equal_after_explicit_reorder']=bool(np.array_equal(m.to_numpy(),b.loc[m.index,m.columns].to_numpy()))
    else: audit['second_matrix']='não fornecida; equivalência não verificada'
    # 17 dígitos significativos permitem round-trip dos floats; fontes permanecem imutáveis.
    m.to_csv(ROOT/'dados_processados/expressao_vst.tsv.gz',sep='\t',float_format='%.17g')
    restored=pd.read_csv(ROOT/'dados_processados/expressao_vst.tsv.gz',sep='\t',index_col=0,float_precision='round_trip')
    if not np.array_equal(m.to_numpy(),restored.to_numpy()): raise ValueError('Falha no round-trip da matriz exportada')
    ids=pd.DataFrame({'patient_id':sorted(set(m.columns)|set(clinical.patient_id)|set(s5.iloc[:,1]))})
    ids['original_vst']=ids.patient_id.where(ids.patient_id.isin(m.columns))
    ids=ids.merge(clinical[['patient_id','original_s22']],how='left')
    ids['original_s5']=ids.patient_id.where(ids.patient_id.isin(s5.iloc[:,1]))
    ids['status']=np.select([ids.original_vst.notna()&ids.original_s22.notna(),ids.original_vst.notna()],['VST–S22 confirmada','VST sem clínica'],default='sem expressão observada')
    ids['rule']='igualdade exata após retirada de sufixo documentado de S22';ids['evidence']='nota † S22';ids['ambiguity']=''
    table(ids,'correspondencia_ids')
    joined=m.T.join(clinical.set_index('patient_id'),how='inner')
    joined.to_csv(ROOT/'dados_processados/clinica_transcriptomica_confirmada.tsv.gz',sep='\t',float_format='%.17g')
    for label,subset in [('antes',clinical),('apos',clinical[clinical.patient_id.isin(m.columns)])]:
        for c in ['sanz','APL9']: table(subset[c].value_counts(dropna=False).rename_axis(c).reset_index(name='n'),f'grupos_{c}_{label}_integracao')
    v=m.var(axis=1);means=m.mean(axis=1);floor=float(m.min().min())
    stats=pd.DataFrame({'ensembl_id':m.index,'mean':means,'median':m.median(axis=1),'variance':v,'std':m.std(axis=1),'near_floor_fraction':np.isclose(m.to_numpy(),floor,atol=1e-4,rtol=0).mean(axis=1)})
    table(stats,'genes_estatisticas');table(stats.sort_values('variance',ascending=False).head(100),'genes_maior_variabilidade');table(stats.sort_values('near_floor_fraction',ascending=False).head(100),'genes_proximos_piso')
    ss=m.describe().T;ss['IQR']=ss['75%']-ss['25%'];table(ss.reset_index(names='patient_id'),'amostras_estatisticas')
    gene_dup=m.duplicated(keep=False);sample_dup=m.T.duplicated(keep=False)
    table(pd.DataFrame({'ensembl_id':m.index,'numerically_duplicated':gene_dup,'zero_variance':v.eq(0),'near_zero_variance':v.lt(1e-8)}),'genes_sinais_auditoria')
    table(pd.DataFrame({'patient_id':m.columns,'numerically_duplicated':sample_dup}),'amostras_sinais_auditoria')
    audit.update(zero_variance=int(v.eq(0).sum()),near_zero_variance=int(v.lt(1e-8).sum()),duplicate_gene_rows=int(gene_dup.sum()),duplicate_sample_columns=int(sample_dup.sum()),floor=floor,near_floor_fraction=float(stats.near_floor_fraction.mean()),candidates_with_expression=int(clinical[clinical.early_death_candidate].patient_id.isin(m.columns).sum()))
    vals=m.to_numpy().ravel();rng=np.random.default_rng(SEED)
    sampled=rng.choice(vals,min(len(vals),100000),replace=False)
    plt.hist(vals,bins=100);plt.xlabel('VST');figure('vst_global_histograma')
    sns.kdeplot(sampled);plt.xlabel('VST; KDE amostra fixa de até 100 mil células');figure('vst_global_densidade')
    for c in m.columns[:4]: sns.histplot(m[c],bins=60,label=c,element='step',fill=False)
    plt.legend();figure('vst_amostras_selecionadas')
    plt.figure(figsize=(18,5));plt.boxplot(m.to_numpy(),showfliers=False);plt.xlabel('Amostra na ordem original');plt.ylabel('VST');figure('vst_boxplot_amostras')
    plt.scatter(means,v,s=2,alpha=.25);plt.xlabel('Média VST');plt.ylabel('Variância');figure('media_variancia')
    annpath=ROOT/'anotacoes/biomart.csv'
    annotation=None
    if annpath.exists():
        annotation=pd.read_csv(annpath,dtype={'ensembl_gene_id':str})
        table(annotation,'anotacao_completa')
        # Um-para-muitos é auditado antes de qualquer seleção de subconjunto.
        table(annotation[annotation.ensembl_gene_id.duplicated(keep=False)],'anotacao_um_para_muitos')
        groups=annotation.groupby('ensembl_gene_id').gene_biotype.agg(lambda x:';'.join(sorted(set(x.dropna()))))
        counts=groups.value_counts().rename_axis('biotype').reset_index(name='n');counts['fraction']=counts.n/counts.n.sum();table(counts,'biotipos')
        counts.set_index('biotype').n.plot.barh(figsize=(9,8));figure('biotipos')
    subsets={'global':m.index}
    if annotation is not None:
        for name in ['protein_coding','lncRNA']:
            eligible=set(annotation.loc[annotation.gene_biotype.eq(name),'ensembl_gene_id'])
            subsets[name]=[g for g in m.index if re.sub(r'\.\d+$','',g) in eligible]
    clinical_index=clinical.set_index('patient_id').reindex(m.columns)
    for subset,genes in subsets.items():
        usable=v.reindex(genes);usable=usable[usable>1e-8].nlargest(1000).index
        k=min(10,len(usable),m.shape[1]-1)
        if k<3: continue
        # PCA centrada, sem dividir pelo desvio: VST já estabiliza a variância.
        model=PCA(n_components=k,svd_solver='full');scores=model.fit_transform(m.loc[usable].T)
        pcs=pd.DataFrame(scores,index=m.columns,columns=[f'PC{i+1}' for i in range(k)])
        table(pcs.reset_index(names='patient_id'),f'pca_{subset}_scores')
        table(pd.DataFrame({'component':pcs.columns,'variance_ratio':model.explained_variance_ratio_}),f'pca_{subset}_variance')
        table(pd.DataFrame(model.components_.T,index=usable,columns=pcs.columns).reset_index(names='ensembl_id'),f'pca_{subset}_loadings')
        plt.plot(range(1,k+1),model.explained_variance_ratio_,'o-');plt.ylabel('Fração da variância');figure(f'pca_{subset}_scree')
        clinical_index['FLT3_any']=clinical_index[['FLT3_ITD','FLT3_TKD']].max(axis=1)
        for c in ['sanz','APL9','NRAS','FLT3_any','sex','early_death_candidate']:
            for pc in ['PC2','PC3']:
                plot=pcs.join(clinical_index[[c]]);plot[c]=plot[c].astype('string').fillna('Sem metadado')
                sns.scatterplot(data=plot,x='PC1',y=pc,hue=c,s=30);plt.title(subset+(' — rótulo inferido' if c=='early_death_candidate' else ''));figure(f'pca_{subset}_{c}_{pc}')
        robust=(pcs-pcs.median()).abs().div((pcs-pcs.median()).abs().median().replace(0,np.nan))
        table(pd.DataFrame({'patient_id':pcs.index,'max_pc_mad_ratio':robust.max(axis=1),'flag_gt6':robust.max(axis=1)>6,'action':'sinalização somente; não excluir'}),f'outliers_{subset}')
        distance=squareform(pdist(m.loc[usable].T,metric='euclidean'))
        tree=linkage(m.loc[usable].T,method='ward')
        sns.clustermap(pd.DataFrame(distance,index=m.columns,columns=m.columns),row_linkage=tree,col_linkage=tree,xticklabels=False,yticklabels=False);figure(f'distancias_{subset}')
        top=usable[:40];z=m.loc[top].sub(m.loc[top].mean(axis=1),axis=0).div(m.loc[top].std(axis=1),axis=0)
        sns.clustermap(z,col_linkage=tree,xticklabels=False,yticklabels=True,figsize=(12,12));figure(f'heatmap_{subset}_40genes')
        tests=[]
        for pc in pcs.columns[:3]:
            ok=clinical_index.age.notna();rho,p=spearmanr(pcs.loc[ok,pc],clinical_index.loc[ok,'age']);tests.append({'component':pc,'variable':'age','rho':rho,'p':p})
        t=pd.DataFrame(tests);t['q_BH']=multipletests(t.p,method='fdr_bh')[1];table(t,f'pca_{subset}_clinica_exploratoria')
    # Assinaturas publicadas: Ensembl diretamente de S16; não seleção independente.
    sig=pd.read_csv(ROOT/'tabelas/assinatura_389_publicada.csv');sigids=sig.iloc[:,0].astype(str)
    common=m.index.intersection(sigids)
    table(pd.DataFrame({'ensembl_id':sigids,'present':sigids.isin(m.index)}),'assinatura_389_disponibilidade')
    if len(common)>=2:
        selected=m.loc[common].var(axis=1).nlargest(30).index
        sns.heatmap(m.loc[selected].T.corr(),vmin=-1,vmax=1,cmap='vlag');figure('coexpressao_30genes_publicados')
    apl9=sig[sig['Gene Symbol'].isin(['GATA1','CLCN4','CACNA2D2','OPTN','KRT1','LTF','S100A12','LEF1','LGALS1'])]
    table(apl9,'APL9_ids_fonte_S16')
    for _,row in apl9.iterrows():
        gene=row.iloc[0]
        if gene in m.index:
            plot=clinical_index.copy();plot['VST']=m.loc[gene].values
            for c in ['sanz','APL9','early_death_candidate']:
                sns.boxplot(data=plot,x=c,y='VST');plt.title(str(row['Gene Symbol'])+'; descrição exploratória');figure(f'APL9_{row["Gene Symbol"]}_{c}')
    table(pd.DataFrame([audit]),'auditoria_expressao')
    return audit
