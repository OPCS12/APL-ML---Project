"""Executa módulos disponíveis e registra explicitamente etapas impedidas."""
import argparse, json, logging
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
from common import ROOT, checksum, table, write_json
from clinical import run as clinical
from expression import run as expression, read_matrix

def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix',type=Path);p.add_argument('--second-matrix',type=Path);a=p.parse_args()
    logging.basicConfig(filename=ROOT/'logs/execucao.log',level=logging.INFO,format='%(asctime)s %(levelname)s %(message)s')
    stamp=datetime.now(timezone.utc).isoformat()
    provenance=[]
    for f in (ROOT/'dados_originais').iterdir():
        if f.is_file(): provenance.append({'file':f.name,'bytes':f.stat().st_size,'sha256':checksum(f),'processed_utc':stamp})
    table(pd.DataFrame(provenance),'procedencia_checksums')
    d,s5,summary=clinical()
    summary['processed_utc']=stamp;summary['expression_status']='BLOQUEADA: anexo acima de 32 MiB; cópia GEO bloqueada pela política de rede'
    if a.matrix:
        if not a.matrix.is_file(): raise FileNotFoundError(a.matrix)
        try:
            summary['expression']=expression(a.matrix,d,s5,a.second_matrix)
            m,_=read_matrix(a.matrix)
            (ROOT/'dados_processados/gene_ids.txt').write_text('\n'.join(m.index)+'\n')
            summary['expression_status']='executada'
        except Exception as e:
            logging.exception('Auditoria da matriz interrompeu análises de expressão')
            summary['expression_status']=f'FALHOU: {e}'
            write_json(ROOT/'relatorios/resumo.json',summary)
            raise
    write_json(ROOT/'relatorios/resumo.json',summary)
    audit=[{'component':'Excel S22','rows':len(d),'status':'auditado','decision':'346 pacientes preservados; notas e linha vazia registradas separadamente'}, {'component':'Excel S5','rows':len(s5),'status':'auditado','decision':'preservado; não equivale à conferência da matriz'}, {'component':'VST','rows':None,'status':summary['expression_status'],'decision':'nenhuma dimensão assumida'}, {'component':'biomaRt','rows':None,'status':'não executado','decision':'R ausente; matriz e acesso Ensembl necessários'}, {'component':'segunda versão VST','rows':None,'status':'não fornecida','decision':'transposição não comparada'}]
    table(pd.DataFrame(audit),'auditoria_integridade')
    report=f'''# EDA de LPA — resultados executados e limites

Processamento UTC: {stamp}. Originais preservados com SHA-256. Sem modelos preditivos.

## Resultados clínicos observados
A S22 contém **{len(d)} pacientes**, enquanto S5 documenta **{len(s5)} amostras de RNA-seq**.
São **{summary['s5_s22_matches']} correspondências S5–S22**; isso não confirma presença na VST.
Sanz: {summary['sanz']}. Mediana de idade: {summary['age_median']} anos.
APL9: {summary['APL9']}, com {summary['APL9_missing']} ausentes.
A nota de rodapé S22 menciona n=321 com sobrevida, mas os rótulos APL9 não nulos devem ser contados diretamente, como feito aqui; não se deve converter a nota em filtro.

## Mortalidade candidata
Há **{summary['candidate_count']} registros sem consolidação**, com distribuição {summary['candidate_APL9']}.
Destes, {summary['candidates_in_s5']} constam em S5. A presença individual na expressão permanece desconhecida.
O artigo descreve 22 mortes na indução e duas perdas de seguimento na coorte de 348; a S22 tem 346 registros.
A concordância agregada não valida um desfecho individual. Ausência de consolidação também pode refletir não realização, recusa, perda de seguimento, omissão de registro ou outras circunstâncias não documentadas.
`early_death_candidate=False` significa consolidação registrada, não comprovação de sobrevivência.
É necessário solicitar aos autores IDs e datas dos eventos, momento da coleta e motivos da ausência antes de modelar mortalidade.

## Consistência e ausências
Foram sinalizados {summary['sanz_categorical_discrepancies']} registros com inconsistências nos testes clínicos especificados.
O limiar de plaquetas difere: ≥40 no suplemento versus >40 na definição de Sanz; as categorias não permitem resolver o valor exato na fronteira.
Zeros mutacionais são categoria observada de ausência da mutação, não dados faltantes. WBC e plaquetas são categorias, não contagens.
Não houve imputação nem exclusão de pacientes. Ver tabelas de completude, padrões de ausência e auditoria individual.
Os testes APL9 são exploratórios, com Fisher em 2×2 e qui-quadrado em demais tabelas, Cramér V e FDR BH dentro dessa família. Células esperadas pequenas limitam a aproximação qui-quadrado. Não demonstram prognóstico independente.

## Expressão e anotação
{summary['expression_status']}.
Genes, amostras reais, valores ausentes, duplicatas numéricas, piso VST, variâncias, outliers, PCA e composição de biotipos **não foram medidos nesta execução**.
A contagem de S5 não substitui a auditoria da matriz. A segunda versão VST não foi fornecida; não é possível afirmar equivalência ou perda de informação.
A anotação biomaRt não foi executada: R não está instalado e a instalação do sistema requer uma capacidade administrativa indisponível. Não foram atribuídos símbolos/biotipos fictícios.
O script R exige release explícito, preserva não mapeados, um-para-muitos e cache. O release usado no estudo não está identificado no PDF principal; requer métodos suplementares. Um release contemporâneo não comprova compatibilidade.
A S16 fornece IDs Ensembl dos 389 genes publicados, exportados como referência, sem tratá-los como seleção independente do Sanz.
Sem FASTQ, contagens brutas e metadados de lote, métricas de alinhamento individuais da matriz e efeitos de lote não podem ser auditados plenamente; S5 contém métricas agregadas úteis, preservadas separadamente.

## Preparação futura e interpretação
Os dados clínicos permitem descrição e auditoria; não permitem nesta execução concluir sobre heterogeneidade transcriptômica ou valor incremental.
Comparar futuramente clínico/Sanz, clínico+APL9 e clínico+assinatura parcimoniosa apenas depois de validar os eventos. Não usar consolidação como preditor ao diagnóstico.
Com cerca de 22 eventos candidatos e milhares de genes, o risco de sobreajuste é alto. Seleção por variância/genes, imputação, escala, redução dimensional e ajuste de hiperparâmetros deverão ser reaprendidos em cada partição de treinamento, com validação interna aninhada ou bootstrap adequado e validação externa.
Relatar ROC-AUC, PR-AUC com prevalência de referência, calibração e intervalos de incerteza na etapa posterior, sem estimá-los agora. Reproduzir Sanz não equivale a predizer morte precoce; associação não demonstra valor incremental.
Nenhum gene foi selecionado definitivamente com base em candidatos de mortalidade da coorte inteira. APL9 e os 389 genes têm histórico de desenvolvimento nesta coorte, limitando a independência de comparações futuras.

## Referências e lacunas documentais
Ver fundamentacao.md. O PDF fornecido é Lin et al. (2021), DOI 10.1158/1078-0432.CCR-20-4375. Jin et al. (2024) não foi fornecido nem identificado bibliograficamente de forma inequívoca; nenhuma afirmação específica foi atribuída a esse trabalho. O guia de Biaggi foi indicado pelo usuário, mas seu conteúdo integral não foi consultado; a organização didática não implica reprodução ou validação desse texto.
'''
    if summary.get('expression_status')=='executada':
        start=report.index('## Expressão e anotação')
        end=report.index('## Preparação futura e interpretação')
        report=report[:start]+'## Expressão e anotação\n\nAuditoria executada: '+json.dumps(summary['expression'],ensure_ascii=False,indent=2)+'\n\nConsulte auditoria_expressao.csv e as tabelas/figuras geradas. Resultados exigem revisão científica; nenhum modelo preditivo foi treinado.\n\n'+report[end:]
    (ROOT/'relatorios/EDA.md').write_text(report)
    flow=f'''flowchart TD\n A["S22: {len(d)} pacientes"] --> C["S5–S22: {summary['s5_s22_matches']} IDs confirmados"]\n B["S5: {len(s5)} registros RNA-seq"] --> C\n A --> D["Sem registro S5: {len(d)-summary['s5_s22_matches']}"]\n C --> E["VST: acesso pendente; elegibilidade não confirmada"]\n A --> F["{summary['candidate_count']} candidatos inferidos; não mortes confirmadas"]\n'''
    (ROOT/'relatorios/fluxo_integracao.mmd').write_text(flow)
    from reporting import render
    render(summary)
    print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__': main()
