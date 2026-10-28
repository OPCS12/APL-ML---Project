"""Gera um notebook que apresenta saídas reais e implementações equivalentes."""
import nbformat as nbf
import json
from common import ROOT
state=json.loads((ROOT/'relatorios/resumo.json').read_text())
complete=state.get('expression_status')=='executada'
nb=nbf.v4.new_notebook()
nb.metadata={'kernelspec':{'display_name':'Python (LPA EDA)','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.12'}}
cells=[]
md=lambda s:cells.append(nbf.v4.new_markdown_cell(s))
code=lambda s:cells.append(nbf.v4.new_code_cell(s))
md('# EDA transcriptômica em LPA\n\nDados reais, auditoria antes de filtragem. **Execução parcial**: clínica e expressão global executadas; biomaRt pendente. Sem modelos preditivos. Consulte `relatorios/pendencias.md`.')
md((ROOT/'relatorios/fundamentacao.md').read_text())
code("from pathlib import Path\nimport sys, json, pandas as pd\nfrom IPython.display import display, Markdown, Image\nROOT=Path.cwd()\nif ROOT.name=='notebooks': ROOT=ROOT.parent\nassert (ROOT/'scripts/run.py').exists(), 'Abra a partir de LPA_EDA ou notebooks'\nsys.path.insert(0,str(ROOT/'scripts'))\nsummary=json.loads((ROOT/'relatorios/resumo.json').read_text())\ndisplay(pd.Series(summary,dtype=object).to_frame('resultado'))")
steps=[
('1 — Estrutura e procedência','Integridade antecede inferência.','Quais arquivos e registros estão efetivamente disponíveis?','SHA-256, inspeção de abas, linhas e cabeçalhos; VST com parsing estrito.','procedencia_checksums','excel_sheets','Arquivos auditados; união de quatro partes com 36.119 genes × 323 amostras, sem alterar os bytes das linhas de dados.','Sem identidade confiável das amostras, não é válido dividir pacientes em treino/teste.'),
('2 — Anotação biomaRt','Símbolos e biotipos dependem de release e montagem.','Quantos IDs são mapeados sem perda?','R/biomaRt, left join, registro de assembly/release e cache por conjunto de IDs.','auditoria_integridade',None,'Não executada: R ausente; matriz já disponível. Nenhum símbolo ou biotipo foi inventado.','Congelar anotação e documentar equivalências antes da modelagem.'),
('3 — Integração de IDs','Correspondência exige evidência documental.','S5 e S22 representam os mesmos pacientes?','Remover somente sufixos A/C/AC documentados na nota †, conservando os IDs originais.','correspondencia_ids',None,'323 correspondências VST–S5–S22 confirmadas; IDs preservados.','Evitar amostras do mesmo paciente em partições distintas.'),
('4 — Características basais','Estatísticas descritivas delimitam a população observada.','Como a coorte se distribui por idade, sexo, Sanz e APL9?','Resumo numérico de idade e frequências com denominador total; WBC e plaquetas mantidos categóricos.','idade_descritiva','tabela1_categoricas','346 pacientes; idade mediana 38 anos. Os grupos não demonstram associação prognóstica independente.','Descrever diferenças de distribuição e não usar consolidação como preditor basal.'),
('5 — Ausências e consistência','Ausência não equivale a teste mutacional negativo.','Quais variáveis têm disponibilidade limitada?','Completude por paciente/variável, padrões e comparação descritiva, sem imputação.','completude_clinica','consistencia_clinica','23 ausências em APL9; o limiar plaquetário não permite recuperar a contagem original.','Imputação futura estimada apenas no treinamento; não presumir mecanismo de ausência.'),
('6 — Distribuição VST','VST não é contagem; piso não significa ausência biológica.','Há valores inválidos, duplicatas, genes invariantes ou amostras atípicas?','Histograma, KDE com semente, estatísticas por gene/amostra e tolerância absoluta 1e-4 para o piso.','auditoria_expressao','amostras_estatisticas','Executada: 36.119 genes × 323 amostras; sem ausentes, perfis duplicados ou genes invariantes. Parâmetros de tolerância e variância próxima de zero (1e-8) são escolhas computacionais, não limites biológicos.','Não excluir outliers automaticamente; avaliar sensibilidade depois de esclarecer a origem.'),
('7 — PCA e agrupamentos','Variância e distâncias descrevem heterogeneidade, não predição.','Existe estrutura global associada a atributos clínicos?','Até 1.000 genes de maior variância positiva; PCA centrada, Ward/euclidiana; heatmap de 40 genes.','pca_global_variance','pca_global_scores','PCA global executada. PCA de protein-coding/lncRNA depende de anotação; ausência de lote limita atribuição biológica.','Refazer seleção/PCA em cada treino quando usados em predição.'),
('8 — Relações clínicas e assinaturas','Associação descritiva deve considerar circularidade e multiplicidade.','Sanz e mutações apresentam estrutura em relação a GI/GII?','Fisher 2×2, qui-quadrado nos demais casos, Cramér V e BH na família declarada; S16 como fonte publicada.','associacoes_APL9_exploratorias','assinatura_389_publicada','Ver resultados calculados abaixo; células esparsas limitam testes assintóticos. Os 389 genes já foram selecionados por grupos clínicos.','Não interpretar reutilização desta coorte como validação independente de APL9.'),
('9 — Mortalidade candidata','Totais concordantes não validam eventos individuais.','Os registros sem consolidação reproduzem os agregados?','Indicador inferido, lista individual e conferência S5–VST, sem converter em óbito confirmado.','candidatos_correspondencia_vst','contingencia_APL9_early_death_candidate','22 candidatos (5 GI e 17 GII), todos na S5 e VST. APL9 é observado, mortalidade é inferida, evento individual confirmado está indisponível.','Não treinar modelos de morte precoce antes de confirmar rótulos e temporalidade.'),
('10 — Planejamento posterior','Poucos eventos e muitos genes favorecem otimismo e instabilidade.','Como comparar valor incremental com avaliação honesta?','Planejar clínico, clínico+APL9 e clínico+assinatura parcimoniosa com validação interna e externa.','auditoria_integridade',None,'Nenhum modelo treinado; nenhuma métrica preditiva calculada.','Pré-processamento dentro do treinamento; ROC-AUC, PR-AUC, calibração e incerteza somente em etapa posterior.')]
for title,theory,obj,method,t1,t2,interp,ml in steps:
 md(f'## {title}\n\n**A. Fundamentação:** {theory}\n\n**B. Objetivo:** {obj}\n\n**C. Método:** {method}\n\n**D. Código executável:** os módulos integrais estão nas células finais e em `scripts/`; as células abaixo apresentam os resultados já produzidos.\n\n**E. Resultados reais:**')
 code("display(pd.read_csv(ROOT/'tabelas/%s.csv'))"%t1+("\ndisplay(pd.read_csv(ROOT/'tabelas/%s.csv'))"%t2 if t2 else ''))
 md(f'**F. Interpretação:** {interp}\n\n**G. Implicações para ML:** {ml}')
md('## Figuras executadas\n\nGráficos clínicos e transcriptômicos calculados com os dados recebidos.')
code("for name in ['idade_histograma','idade_sanz','matriz_ausencia','completude','sanz_APL9','candidatos_APL9','vst_global_histograma','pca_global_sanz_PC2','pca_global_APL9_PC2','heatmap_global_40genes']:\n    display(Image(filename=str(ROOT/'figuras'/f'{name}.png')))" )
md((ROOT/'relatorios/auditoria_uniao_partes.md').read_text())
md('## Relatório acadêmico com resultados reais')
code("display(Markdown((ROOT/'relatorios/EDA.md').read_text()))")
md('## Código modular completo\n\nAs próximas células definem os módulos equivalentes. Não executam chamadas de rede ou análise de expressão automaticamente. Para reproduzir: `python scripts/run.py`; para VST real, acrescentar `--matrix CAMINHO`.')
for f in ['common.py','clinical.py','expression.py']:
 md('### '+f);code((ROOT/'scripts'/f).read_text().replace('Path(__file__).resolve().parents[1]','ROOT'))
md('### Script R biomaRt\n\n```r\n'+(ROOT/'scripts/annotate_biomart.R').read_text()+'\n```')
md('### Orquestração Python\n\n```python\n'+(ROOT/'scripts/run.py').read_text()+'\n```')
md((ROOT/'relatorios/pendencias.md').read_text())
nb.cells=cells
nbf.write(nb,ROOT/'notebooks/LPA_EDA.ipynb')
