# LPA_EDA

Projeto em Python para auditoria e EDA descritiva de LPA; R/biomaRt para anotação. Não treina modelos preditivos. Comece por [relatorios/EDA.md](relatorios/EDA.md), que distingue resultados executados e etapas pendentes. A fundamentação precede o código no notebook.

## Situação desta execução

As quatro partes VST foram recebidas, preservadas e unidas na sequência 01–04, com um único cabeçalho e linhas de dados idênticas byte a byte às recebidas. A matriz conjunta contém **36.119 genes × 323 amostras**; IDs e ordem não foram alterados. Consulte [auditoria da união](relatorios/auditoria_uniao_partes.md). O script `scripts/merge_vst_parts.py` reproduz a união e recusa sobrescrever arquivos existentes.

Clínica, expressão global, PCA, agrupamentos e visualizações dos nove genes APL9 foram executados com os dados reais. Há 323 correspondências VST–S22, incluindo os 22 candidatos inferidos (5 GI e 17 GII), sem confirmação individual de óbito. APL320 e APL321 são colunas separadas nestas partes. Os originais e a união foram preservados somente leitura.

R/biomaRt ainda está indisponível; composição de biotipos e PCA por protein-coding/lncRNA permanecem pendentes. Não há segunda versão independente VST anexada. Não se comprovou equivalência com o arquivo integral inacessível enviado anteriormente. Ver `relatorios/pendencias.md`.

## Instalação local

Use Python 3.12. Na pasta LPA_EDA:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/run.py
```

Na máquina preparada, use `/workspace/.venvs/lpa-eda/bin/python` (não depende de ativação). Configure caches graváveis se necessário:

```bash
export MPLCONFIGDIR=/workspace/.cache/lpa-eda/matplotlib
export XDG_CACHE_HOME=/workspace/.cache/lpa-eda
mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
```

As versões em requirements.txt correspondem ao ambiente validado; não existe lockfile anterior. `logs/versions.json` registra as versões importadas. Originais devem ser mantidos fora do controle público de versão quando houver restrições de compartilhamento.

## Dados e auditoria estrutural

Os quatro arquivos originais estão em `dados_originais/GSE172057_VST_parte_0N.txt.gz`; a matriz unificada com representação textual preservada está em `dados_processados/GSE172057_VST_unificada.txt.gz`. SHA-256 e limites dos blocos constam em `relatorios/auditoria_uniao_partes.json`. A união foi validada por segunda leitura integral, hash de conteúdo e sequência de IDs. Nenhuma linha gênica foi removida ou reordenada; somente três cabeçalhos redundantes foram omitidos após confirmação de identidade.

Para reproduzir a análise a partir da matriz disponível:

```bash
python scripts/run.py --matrix dados_processados/GSE172057_VST_unificada.txt.gz
```

A execução sem `--matrix` refaz somente a clínica e assinala que a expressão não foi executada nessa rodada. Para comparar duas fontes independentes reais, use também `--second-matrix CAMINHO`; não confunda uma transposição construída pelo próprio código com validação independente.

O leitor recusa campos inconsistentes, cabeçalhos vazios/duplicados, IDs ambíguos, dados ausentes/não finitos e duplicatas de IDs. A matriz recebida passou por essas verificações sem correções. As tabelas e figuras derivadas podem ser atualizadas por novas execuções; os originais e a matriz unificada não são sobrescritos. O processamento conjunto não comprova que o arquivo integral anterior continha exatamente os mesmos bytes: isso exige acesso à fonte anterior ou seu checksum descompactado.

## Anotação R

Instale R e o pacote seguindo as instruções oficiais do Bioconductor:

```r
install.packages('BiocManager', repos='https://cloud.r-project.org')
BiocManager::install('biomaRt', ask=FALSE, update=FALSE)
```

Depois da auditoria bem-sucedida da matriz, obtenha dos métodos suplementares o release compatível. O exemplo abaixo usa release 103 apenas como anotação contemporânea exploratória, sem alegar ser a original:

```bash
Rscript scripts/annotate_biomart.R dados_processados/gene_ids.txt 103
python scripts/run.py --matrix dados_processados/GSE172057_VST_unificada.txt.gz
```

Não há chamada automática ao serviço sem release explícito. Cache por IDs/release, metadados de assembly e `sessionInfo()` ficam em anotacoes/. Campos ausentes e relações um-para-muitos permanecem na tabela. O serviço deve estar acessível por HTTPS com TLS verificado. Se falhar, continue a EDA global sem inventar biotipos. PCA de protein-coding/lncRNA depende de anotação obtida; o método requer ao menos quatro amostras e três genes variáveis.

## Notebook e resultados

Abra `notebooks/LPA_EDA.ipynb` com Jupyter em ambiente com as dependências instaladas. O notebook chama os mesmos módulos, mostra tabelas e imagens reais, e assinala as etapas não executadas. `scripts/run.py` refaz os resultados clínicos e o relatório; figuras e tabelas geradas podem ser sobrescritas por uma nova execução, nunca os originais. Consulte o timestamp UTC em resumo.json. Revise os logs se uma execução falhar: arquivos de resultados de execuções antigas não demonstram sucesso atual.

## Limites e próximos passos

Antes de inferências prognósticas, confirmar com os autores os eventos individuais, os dois perdidos de seguimento, IDs da VST e release de anotação. Consolidação é pós-diagnóstico e não deve ser preditor. Os 22 candidatos são rótulos inferidos; 0 de mutação é observado negativo. WBC/plaquetas são categorias. APL9 e os 389 genes não são assinaturas desenvolvidas independentemente desta coorte.

Os requisitos finais ainda não integralmente atendidos estão em `relatorios/pendencias.md`. Não publique o ambiente como totalmente configurado até resolver os requisitos de expressão e R e validar a execução correspondente.
