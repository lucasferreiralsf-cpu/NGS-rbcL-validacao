# Pipeline de análise NGS rbcL

Este documento descreve o fluxo de análise aplicado aos dados FASTQ obtidos via Illumina MiSeq para validação do marcador **rbcL** em plantas.

---

## 1. Controle de Qualidade (QC)
- Ferramentas: [FastQC](ca://s?q=Como_usar_FastQC), [Rqc](ca://s?q=Como_usar_pacote_Rqc_no_R).
- Objetivo: avaliar qualidade das leituras, identificar adaptadores, checar distribuição de bases.

---

## 2. Trimming
- Ferramentas: [Trimmomatic](ca://s?q=Como_usar_Trimmomatic), [Cutadapt](ca://s?q=Como_usar_Cutadapt).
- Objetivo: remover adaptadores e bases de baixa qualidade.

---

## 3. Alinhamento
- Ferramentas: [BLAST](ca://s?q=Como_usar_BLAST), [BOLD](ca://s?q=Como_usar_BOLD).
- Objetivo: alinhar sequências contra banco de referência para identificar espécies.

---

## 4. Identificação
- Base: marcador **rbcL**.
- Objetivo: determinar espécies vegetais presentes nas amostras.

---

## Estrutura de saída
- `results/qc/` → relatórios de qualidade.
- `results/trimming/` → arquivos pós-trimming.
- `results/alignment/` → resultados de alinhamento.
- `results/identification/` → tabelas e gráficos finais.
