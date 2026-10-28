# Registro de decisões

1. PDF e Excel originais copiados e tornados somente leitura; SHA-256 registrado; nenhuma célula corrigida.
2. S22 lida com título e cabeçalho separados; todas as linhas fora do padrão de paciente registradas individualmente, incluindo notas e linha vazia.
3. Remoção apenas de sufixos A/C/AC, com evidência na nota †; IDs originais preservados. Duplicatas normalizadas interrompem integração.
4. S5 usada como evidência documental de RNA-seq; não como substituição da matriz.
5. Nenhuma amostra, gene ou paciente excluído por qualidade; nenhum valor imputado.
6. Consolidação ausente produz exclusivamente `early_death_candidate`, sem transformar em morte confirmada.
7. Grupos Sanz publicados preservados mesmo quando há divergência da reconstrução por categorias.
8. VST não transferida: limite do executor de 32 MiB. Consulta pública via NCBI negada pelo proxy. Sem conclusão sobre suas dimensões ou integridade.
9. Segunda VST não fornecida. Não se compara uma transposição construída pelo próprio código como se fosse uma segunda fonte independente.
10. R ausente; sudo não disponível e apt-get sem permissão. biomaRt preparado, não executado; release original não confirmado.
11. Sem treino de modelos, métricas preditivas ou seleção definitiva de genes. Testes clínicos exploratórios possuem família BH explícita.
12. Repositório sem commits/referências remotas no teste realizado. Novos arquivos pertencem ao projeto solicitado pelo texto anexado, separado das restrições de preservação de arquivos preexistentes do onboarding.


## Atualização após recebimento das partes

13. Quatro arquivos gzip independentes recebidos e lidos integralmente com CRC válido. Cada parte repete o mesmo cabeçalho, idêntico byte a byte.
14. União na ordem 01–04, sem ordenar os genes: 10.000 + 10.000 + 10.000 + 6.119 = 36.119 registros. Três cabeçalhos redundantes omitidos; nenhuma linha gênica omitida, alterada, agregada ou repetida pela união.
15. Saída lida novamente e comparada por SHA-256 das linhas descompactadas e hash da sequência de IDs; verificação aprovada. Originais e unificada somente leitura.
16. Cabeçalho atual com ID + 323 amostras, começando em APL2. APL320 e APL321 separados nas partes recebidas; isso não permite inferir qual era o cabeçalho da fonte anterior inacessível.
17. Leitura numérica estrita, exportação com round-trip, integração VST–S22, estatísticas VST, PCA global, agrupamentos e assinatura S16 executados. Sem exclusões por qualidade nem treino de modelos. Todos os 22 candidatos possuem correspondência VST, mas continuam sem confirmação de óbito.
18. R/biomaRt e PCA por biotipo permanecem pendentes; resultados antigos de bloqueio de transferência em logs descrevem o momento anterior e não invalidam o acesso agora obtido por partes.
