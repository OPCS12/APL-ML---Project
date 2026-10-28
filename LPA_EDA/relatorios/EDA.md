# EDA de LPA — resultados observados

Processamento UTC: 2026-10-27T23:43:52.957935+00:00. Originais preservados; nenhuma imputação, modelo preditivo ou exclusão de pacientes por qualidade. Consulte `auditoria_uniao_partes.md` para a reconstrução exata da matriz.

## Coorte clínica
A Tabela S22 contém **346 pacientes**, e a S5 contém **323 registros RNA-seq**, todos com correspondência na S22 após remoção apenas dos sufixos A/C/AC documentados na nota †. IDs originais foram conservados. A S22 inclui 23 pacientes sem registro S5; não se presume presença na VST.
Sanz na S22: {'Intermediate-risk': 172, 'High-risk': 109, 'Low-risk': 65}. Mediana de idade: 38.0 anos. APL9: {'GII': 163, 'GI': 160}, com 23 ausências.
A nota § menciona n=321 com sobrevida, mas há 323 rótulos APL9 observados na tabela: não se deve converter a nota em filtro nem presumir eventos/datas disponíveis.

As contagens WBC e plaquetas só estão disponíveis como categorias. O limiar de plaquetas do suplemento é ≥40, enquanto a definição de Sanz do artigo usa >40. Os grupos originais foram preservados; a reconstrução categórica sinalizou 0 inconsistências nos testes implementados.
Zero mutacional indica resultado negativo observado, não ausência de dado. Não houve imputação. Tabelas de completude, padrões de ausência e estatísticas por grupo estão em `tabelas/`.

## Mortalidade precoce candidata
Há **22 pacientes sem registro de consolidação**, classificados exclusivamente como `early_death_candidate`: {'GII': 17, 'GI': 5}. Todos os 22 constam na S5.
O artigo relata 22 mortes durante indução e dois perdidos de seguimento entre 348 pacientes. A S22 tem 346 registros; a concordância agregada não confirma identidade dos óbitos. Ausência de consolidação também pode refletir omissão, não realização, recusa ou outras circunstâncias não documentadas. False significa somente consolidação registrada.

**22 dos 22 candidatos possuem ID inequívoco na VST**. Isso confirma a correspondência das amostras, não o desfecho.

## Auditoria estrutural e integração
A matriz conjunta tem **36,119 genes × 323 amostras**. Orientação: genes nas linhas. `ID` é a coluna dos genes; a primeira amostra é APL2. Nas quatro partes recebidas, APL320 e APL321 são cabeçalhos separados, sem ambiguidade estrutural. Não houve renomeação de amostras nem correção de valores.
A verificação numérica identificou **0 valores ausentes/não numéricos/infinitos**; IDs duplicados foram bloqueados pelo leitor e não foram encontrados. Há 323 correspondências VST–S22 e 0 amostras VST sem clínica. A tabela individual de correspondência inclui também todos os pacientes sem expressão, sem exclusão silenciosa.
O SHA-256 das linhas descompactadas e o hash da sequência de genes confirmam a preservação byte a byte na união 01–04. Os três cabeçalhos redundantes foram omitidos após verificação de identidade. As partes e o arquivo unificado permanecem somente leitura.
O arquivo exportado `expressao_vst.tsv.gz` teve round-trip numérico verificado. A matriz principal unificada conserva inclusive a representação textual dos valores. A equivalência com o anexo integral anterior e com uma segunda versão transposta não foi verificada: essas fontes continuam indisponíveis.

## Distribuição e variabilidade VST
O mínimo observado é **3.44844073036933**; 24.29% das células estão a até 1e-4 desse mínimo. Isso descreve o piso numérico da transformação, sem significar ausência de transcrição ou zero bruto.
Há 0 genes de variância nula e 0 abaixo do limiar computacional 1e-8 (incluindo os de variância nula). Foram sinalizadas 0 linhas gênicas pertencentes a conjuntos de perfis numericamente idênticos e 0 colunas de amostras com perfil idêntico. Essas contagens incluem todos os membros dos conjuntos, não somente repetições adicionais. Nenhum foi removido.
As medianas de expressão por amostra variam de 4.2661 a 5.7660. 244 genes têm pelo menos 90% das amostras próximas ao piso. As tabelas apresentam médias, medianas, IQR, variância, desvio e genes mais variáveis, sem inferência prognóstica.

## PCA, agrupamentos e assinaturas
A PCA global usa até 1.000 genes com maior variância acima de 1e-8, centrados e sem padronização por desvio-padrão. A fração de variância explicada por PC1, PC2 e PC3 é, respectivamente, 20.85%, 13.03% e 6.92%. Isso mede a variância no subconjunto selecionado, não em todos os genes.
Foram produzidos escores, loadings, scree plot e projeções separadas por Sanz, GI/GII, NRAS, FLT3, sexo e candidatos inferidos. Ward/euclidiana e heatmaps usam o mesmo subconjunto; expressão padronizada é mostrada para os 40 genes mais variáveis. O critério exploratório de distância à mediana superior a 6 MAD nas componentes sinalizou 2 amostras; ele não é diagnóstico de erro e nenhuma foi excluída.
Resumos das componentes por grupos estão disponíveis. Diferenças entre médias ou agrupamentos visuais não demonstram separação clínica validada, causalidade ou valor incremental. Sem informação de lote, não é possível decompor variação biológica e técnica adequadamente.
Dos 389 IDs da assinatura publicada S16, 389 foram encontrados. 9 genes APL9 foram identificados pela correspondência documentada Ensembl–símbolo da S16, sem substituir a anotação global solicitada. Suas distribuições por grupos e coexpressão foram visualizadas; não se recalculou o escore publicado usando coeficientes sobre uma escala sem compatibilidade estabelecida. Os 389 genes foram previamente selecionados por grupos clínicos e não constituem seleção independente de Sanz.

## Anotação e limitações
R/biomaRt ainda não foi executado: R está ausente e a máquina não oferece instalação administrativa. Não foram inventados símbolos nem biotipos. A composição global de protein-coding, lncRNAs e demais classes, a tabela completa de anotação e as PCAs por biotipo permanecem pendentes. O script R conserva não mapeados, um-para-muitos e cache; é necessário confirmar o release e a montagem do estudo nos métodos suplementares. Usar um release contemporâneo não comprova compatibilidade.
Sem FASTQ e contagens brutas, a EDA VST não permite auditar qualidade das bases, cobertura ou alinhamento diretamente, nem usar DESeq2 como se fossem contagens. As métricas S5 foram preservadas como evidência documental independente.
Os testes clínicos são exploratórios: Fisher em 2×2, qui-quadrado nos demais, tamanho de efeito Cramér V e FDR BH na família declarada. Células esparsas limitam aproximações assintóticas. PCA–idade usa Spearman e BH em sua família. Nenhum desses resultados valida prognóstico.

## Implicações para modelagem posterior
É indispensável validar os eventos individuais com os autores, incluindo datas e momento da coleta. Consolidação ocorre após o diagnóstico e não deve ser preditor basal. Cerca de 22 eventos potenciais e milhares de genes tornam a seleção adaptativa muito instável: considerar assinaturas parcimoniosas, com seleção, imputação, escala e redução dimensional reaprendidas no treinamento de cada partição. Comparar futuramente clínico/Sanz, clínico+APL9 e clínico+assinatura somente com validação interna apropriada e externa. APL9 foi desenvolvida nesta coorte, o que limita a independência de uma avaliação interna posterior.
Relatar ROC-AUC, PR-AUC com referência de prevalência, calibração e incerteza na fase posterior. Reproduzir Sanz não equivale a predizer morte precoce. Não houve treino ou avaliação de modelos nem seleção definitiva de genes pelo desfecho nesta EDA.

## Referências
Ver fundamentacao.md. Lin et al. (2021), DOI 10.1158/1078-0432.CCR-20-4375, foi inspecionado no PDF fornecido. Jin et al. (2024) não foi identificado de forma inequívoca: requer título/DOI/PDF. O guia Biaggi foi indicado pelo usuário, mas seu conteúdo integral não foi consultado. Não se atribuíram resultados não verificados a essas fontes.
