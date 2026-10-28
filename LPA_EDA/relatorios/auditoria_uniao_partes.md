# Auditoria da união das partes VST

As partes foram lidas como arquivos gzip independentes e unidas na sequência 01 → 02 → 03 → 04, sem ordenação dos genes. O cabeçalho das quatro partes é idêntico byte a byte: foi escrito uma vez, e as outras três ocorrências foram registradas e omitidas. `ID` continua sendo a coluna de identificadores dos genes.

| Parte | Genes | Primeiro ID | Último ID |
| --- | ---: | --- | --- |
| GSE172057_VST_parte_01.txt.gz | 10,000 | ENSG00000000003 | ENSG00000159917 |
| GSE172057_VST_parte_02.txt.gz | 10,000 | ENSG00000159921 | ENSG00000219395 |
| GSE172057_VST_parte_03.txt.gz | 10,000 | ENSG00000219409 | ENSG00000255362 |
| GSE172057_VST_parte_04.txt.gz | 6,119 | ENSG00000255363 | ENSG00000273493 |

Resultado: **36,119 genes × 323 amostras**, além da coluna `ID`. Há 36,119 identificadores de genes únicos; nenhum ID de gene foi removido, agregado ou reordenado. As colunas de amostras começam em APL2 e terminam em APL379, preservando todos os intervalos e a ordem recebida. APL320 e APL321 são colunas separadas nestas partes. Isso descreve os arquivos recebidos agora, sem estabelecer o formato do arquivo integral inacessível enviado anteriormente.

Todos os registros têm 324 campos separados por tabulação. O CRC gzip das quatro partes e da saída foi verificado pela leitura completa. As linhas de dados foram copiadas como bytes, sem conversão ou arredondamento de números. Uma segunda leitura confirmou que o conteúdo descompactado é idêntico ao cabeçalho único seguido dos quatro blocos originais, e que a sequência de IDs foi preservada.

- SHA-256 do conteúdo descompactado: `0032c4f8beab4dead030512beb1ca7966c9902bc0e88ac54e25dc448db489627`.
- SHA-256 da sequência de IDs: `679f2d145d27813aa4e9efa484614a537c33aed7848cddf8d7c1570030011726`.
- SHA-256 do gzip de saída: `1328d23813e337239c675e9c5a01fd67bc69a54179145e06368ca191ff6bb86e`.

Saída: `dados_processados/GSE172057_VST_unificada.txt.gz`. Os quatro arquivos recebidos foram preservados em `dados_originais/` com permissão somente leitura. A saída também foi marcada somente leitura. Os checksums individuais, limites de cada bloco e decisões constam em `auditoria_uniao_partes.json`.

Não é possível comprovar equivalência com o arquivo integral anterior sem acesso a esse arquivo ou ao seu checksum descompactado. A ausência de genes ou amostras em relação ao estudo original precisa de evidência externa; a união garante preservação integral destas quatro partes.

Para reproduzir a união em outro caminho (o script recusa sobrescrever saída existente):

```bash
python scripts/merge_vst_parts.py --parts dados_originais/GSE172057_VST_parte_01.txt.gz dados_originais/GSE172057_VST_parte_02.txt.gz dados_originais/GSE172057_VST_parte_03.txt.gz dados_originais/GSE172057_VST_parte_04.txt.gz --output dados_processados/uniao_reproduzida.txt.gz --audit relatorios/auditoria_reproduzida.json
```
