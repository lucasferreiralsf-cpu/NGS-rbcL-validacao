# Validaçao NGS rbcL 

Este repositório documenta a análise de dados NGS Illumina MiSeq para validação do marcador **rbcL** em plantas.

## Objetivos
- Realizar controle de qualidade (QC) das leituras.
- Aplicar trimming para remover adaptadores e bases de baixa qualidade.
- Alinhar sequências contra banco de referência (GenBank/BOLD).
- Identificar espécies vegetais com base no marcador rbcL.

## Estrutura
- `data/` → arquivos FASTQ originais (não versionados).
- `scripts/` → scripts em R para cada etapa da análise.
- `results/` → relatórios e gráficos gerados.
- `docs/` → documentação e protocolos.

## Ferramentas
- R + pacotes [ShortRead](ca://s?q=Como_usar_pacote_ShortRead_no_R), [Rqc](ca://s?q=Como_usar_pacote_Rqc_no_R), [Biostrings](ca://s?q=Como_usar_pacote_Biostrings_no_R).
- FastQC, Trimmomatic, BLAST.

## Fluxo inicial
1. Baixar dados do BaseSpace.
2. QC com FastQC/Rqc.
3. Trimming com Trimmomatic/Cutadapt.
4. Alinhamento/identificação com BLAST ou BOLD.


