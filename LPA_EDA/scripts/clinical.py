"""Auditoria clínica: conserva IDs originais e só normaliza sufixos documentados."""
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests
from common import ROOT, table, figure

COLS=['original_s22','sex','age','sanz','consolidation','wbc_gt10','platelets_ge40','FLT3_ITD','FLT3_TKD','WT1','NRAS','ARID1A','RREB1','KRAS','EBF3','PML_RARA','APL9']

def run():
    path=ROOT/'dados_originais/suplemento_lin_2021.xlsx'
    x=pd.ExcelFile(path)
    table(pd.DataFrame([{'sheet':s,'rows_including_title':pd.read_excel(x,s,header=None).shape[0], 'columns':pd.read_excel(x,s,header=None).shape[1]} for s in x.sheet_names]),'excel_sheets')
    raw=pd.read_excel(x,'Table S22',header=None)
    table(raw.head(8),'s22_primeiras_linhas');table(raw.tail(8),'s22_ultimas_linhas')
    # A primeira linha é título, a segunda cabeçalho; seleção explícita de pacientes.
    body=raw.iloc[2:].copy()
    keep=body.iloc[:,0].astype(str).str.fullmatch(r'APL\d+(?:AC|A|C)?')
    removed=body.loc[~keep].copy();removed['reason']='linha vazia ou nota; não é paciente'
    table(removed,'s22_linhas_nao_pacientes')
    d=body.loc[keep].copy();d.columns=COLS
    d.insert(1,'patient_id',d.original_s22.str.replace(r'(AC|A|C)$','',regex=True))
    if d.patient_id.duplicated().any(): raise ValueError('IDs S22 normalizados duplicados: revisão necessária')
    for c in ['age',*COLS[5:15]]: d[c]=pd.to_numeric(d[c],errors='raise')
    # Não é desfecho confirmado; False significa apenas consolidação registrada.
    d['early_death_candidate']=d.consolidation.isna()
    d['sanz_binary']=d.sanz.map({'High-risk':'HR','Low-risk':'não HR','Intermediate-risk':'não HR'})
    s5=pd.read_excel(x,'Table S5',header=1)
    table(s5,'s5_preservada')
    s5ids=s5.iloc[:,1].astype(str)
    if s5ids.duplicated().any(): raise ValueError('IDs S5 duplicados')
    all_ids=sorted(set(s5ids)|set(d.patient_id),key=lambda s:int(s[3:]))
    mapping=pd.DataFrame({'patient_id':all_ids}).merge(d[['patient_id','original_s22']],how='left').merge(pd.DataFrame({'patient_id':s5ids,'original_s5':s5ids}),how='left')
    mapping['original_vst']=pd.NA
    mapping['status']=np.where(mapping.original_s5.notna() & mapping.original_s22.notna(),'S5–S22 confirmada; VST pendente','presença em somente uma tabela; VST pendente')
    mapping['rule']='remoção apenas do sufixo A/C/AC de S22; igualdade exata do ID restante'
    mapping['evidence']='nota † da Tabela S22'
    mapping['ambiguity']='arquivo VST indisponível para auditoria'
    table(mapping,'correspondencia_ids')
    table(d,'clinica_completa');d.to_csv(ROOT/'dados_processados/clinica.csv',index=False)
    completeness=pd.DataFrame({'variable':COLS[1:],'observed':[d[c].notna().sum() for c in COLS[1:]],'missing':[d[c].isna().sum() for c in COLS[1:]],'missing_fraction':[d[c].isna().mean() for c in COLS[1:]]})
    table(completeness,'completude_clinica')
    table(pd.DataFrame({'patient_id':d.patient_id,'missing':d[COLS[1:]].isna().sum(axis=1)}),'ausencia_por_paciente')
    patterns=d[COLS[1:]].isna().astype(int).astype(str).agg(''.join,axis=1).value_counts().rename_axis('pattern_1_missing').reset_index(name='n')
    table(patterns,'padroes_ausencia')
    numeric=d.age.describe(percentiles=[.25,.5,.75]).to_frame().T
    numeric['missing']=d.age.isna().sum();numeric['IQR']=numeric['75%']-numeric['25%'];table(numeric.reset_index(names='variable'),'idade_descritiva')
    rows=[]
    categorical=[c for c in COLS[1:] if c!='age']+['sanz_binary','early_death_candidate']
    for c in categorical:
        for level,n in d[c].fillna('Ausente').value_counts(dropna=False).items():
            rows.append({'variable':c,'level':level,'n':n,'percent_total':100*n/len(d),'denominator':len(d)})
    table(pd.DataFrame(rows),'tabela1_categoricas')
    for group in ['sanz','APL9','early_death_candidate']:
        table(d.groupby(group,dropna=False).age.describe().reset_index(),f'idade_por_{group}')
        for c in categorical:
            if c!=group: table(pd.crosstab(d[c].fillna('Ausente'),d[group].fillna('Ausente')).reset_index(),f'contingencia_{c}_{group}')
    candidates=d[d.early_death_candidate]
    table(candidates,'mortalidade_candidatos_nao_confirmados')
    table(candidates[['patient_id','original_s22','sanz','APL9']].merge(pd.DataFrame({'patient_id':s5ids,'s5_present':True}),how='left'),'candidatos_correspondencia_s5')
    # Discordâncias são sinalizadas, não corrigidas: ≥40 na planilha versus >40 na definição.
    issues=[]
    for _,r in d.iterrows():
        expected='High-risk' if r.wbc_gt10==1 else ('Low-risk' if r.platelets_ge40==1 else 'Intermediate-risk')
        if r.sanz!=expected: issues.append({'patient_id':r.patient_id,'issue':'Sanz difere da reconstrução categórica; investigar limiar de 40','observed':r.sanz,'reconstruction':expected})
        if not 0<=r.age<=120: issues.append({'patient_id':r.patient_id,'issue':'idade implausível'})
        for c in COLS[5:15]:
            if pd.notna(r[c]) and r[c] not in [0,1]: issues.append({'patient_id':r.patient_id,'issue':f'código inesperado {c}: {r[c]}'})
    table(pd.DataFrame(issues,columns=['patient_id','issue','observed','reconstruction']),'consistencia_clinica')
    # Testes exploratórios com tamanhos de efeito; não constituem prognóstico.
    tests=[]
    for c in ['sanz','FLT3_ITD','FLT3_TKD','WT1','NRAS','ARID1A','RREB1','KRAS','EBF3','sex']:
        ct=pd.crosstab(d[c],d.APL9)
        chi,p,_,exp=chi2_contingency(ct)
        method='chi-quadrado';effect=(chi/(ct.to_numpy().sum()*min(ct.shape[0]-1,ct.shape[1]-1)))**.5
        if ct.shape==(2,2): _,p=fisher_exact(ct);method='Fisher bilateral'
        tests.append({'variable':c,'n_complete':int(ct.to_numpy().sum()),'test':method,'p':p,'cramers_v':effect,'min_expected':float(exp.min())})
    tests=pd.DataFrame(tests);tests['q_BH']=multipletests(tests.p,method='fdr_bh')[1];table(tests,'associacoes_APL9_exploratorias')
    # Figuras com denominadores explícitos; dados ausentes aparecem como categoria.
    sns.heatmap(d[COLS[1:]].isna(),cbar=False);plt.ylabel('Pacientes S22');figure('matriz_ausencia')
    completeness.set_index('variable').missing_fraction.plot.barh();plt.xlabel('Proporção ausente');figure('completude')
    d.age.plot.hist(bins=20);plt.xlabel('Idade (anos)');figure('idade_histograma')
    sns.boxplot(data=d,x='sanz',y='age');figure('idade_sanz')
    for c in ['sanz','APL9','sex','consolidation']:
        d[c].fillna('Ausente').value_counts().plot.bar();plt.ylabel(f'Pacientes (N={len(d)})');figure('distribuicao_'+c)
    ct=pd.crosstab(d.sanz,d.APL9);ct.div(ct.sum(axis=1),axis=0).plot.bar(stacked=True);plt.ylabel('Proporção entre APL9 observados');figure('sanz_APL9')
    pd.crosstab(d.sanz,d.early_death_candidate).plot.bar(stacked=True);plt.ylabel('Pacientes; candidato ≠ morte confirmada');figure('candidatos_sanz')
    sns.heatmap(pd.crosstab(d.APL9,d.early_death_candidate),annot=True,fmt='d');figure('candidatos_APL9')
    table(pd.read_excel(x,'Table S16',header=1),'assinatura_389_publicada')
    for c in ['PML_RARA','APL9']:
        table(pd.crosstab(d.sanz,d[c].isna()).reset_index(),f'ausencia_{c}_sanz')
    summary={'s22_patients':len(d),'s5_samples':len(s5ids),'s5_s22_matches':len(set(s5ids)&set(d.patient_id)), 'clinical_without_s5':sorted(set(d.patient_id)-set(s5ids)), 'candidate_count':len(candidates),'candidate_APL9':candidates.APL9.value_counts().to_dict(),'candidates_in_s5':int(candidates.patient_id.isin(s5ids).sum()),'sanz':d.sanz.value_counts().to_dict(),'APL9':d.APL9.value_counts().to_dict(),'APL9_missing':int(d.APL9.isna().sum()),'sanz_categorical_discrepancies':len(issues),'age_median':float(d.age.median())}
    return d,s5,summary
