################################################################################
# ANÁLISE DE DADOS NGS DO MARCADOR rbcL UTILIZANDO DADA2
#
# Objetivo:
# Processar dados de sequenciamento Illumina MiSeq (2x300 bp) para obtenção
# de ASVs do marcador rbcL e validação taxonômica por BLAST.
#
# Amostras analisadas:
# - Controle Branco (BCO)
# - Algodão
# - Milho
# - Soja
#
# Observações da análise:
# - Os primers amplificaram corretamente o gene rbcL.
# - Leituras R1 e R2 apresentaram alinhamento ao rbcL com ~99,7% de identidade.
# - A baixa eficiência de merge sugere sobreposição limitada entre R1 e R2.
# - Os dados indicam que uma abordagem single-end (R1) pode ser mais robusta.
################################################################################


################################################################################
# ETAPA 1 - INSTALAÇÃO DOS PACOTES (EXECUTAR APENAS UMA VEZ)
################################################################################

# Instala o gerenciador do Bioconductor
install.packages("BiocManager")

# Instala o DADA2
BiocManager::install("dada2")

# Pacotes auxiliares
install.packages(c(
  "ggplot2",
  "tidyverse"
))

# Leitura de arquivos FASTQ
BiocManager::install("ShortRead")

# Exportação de FASTA
install.packages("seqinr")


################################################################################
# ETAPA 2 - CARREGAR PACOTES
################################################################################

library(dada2)
library(ShortRead)
library(ggplot2)
library(seqinr)


################################################################################
# ETAPA 3 - DEFINIR O DIRETÓRIO DOS DADOS
################################################################################

# Ajustar o caminho para a pasta contendo os FASTQ

path <- "C:/Users/lucas.ferreira/OneDrive - Intecso Soluções e Inovações em Agronegócios LTDA/Área de Trabalho/NGS_rbcL"

# Verificar se a pasta existe

dir.exists(path)

# Visualizar os arquivos disponíveis

list.files(path)


################################################################################
# ETAPA 4 - IDENTIFICAR OS ARQUIVOS FORWARD (R1) E REVERSE (R2)
################################################################################

# Leituras Forward

fnFs <- sort(
  list.files(
    path,
    pattern = "_R1_001.fastq.gz$",
    full.names = TRUE
  )
)

# Leituras Reverse

fnRs <- sort(
  list.files(
    path,
    pattern = "_R2_001.fastq.gz$",
    full.names = TRUE
  )
)

# Conferir

fnFs
fnRs


################################################################################
# ETAPA 5 - AVALIAÇÃO DA QUALIDADE DAS LEITURAS
################################################################################

# Forward

plotQualityProfile(fnFs[1])

# Reverse

plotQualityProfile(fnRs[1])

# Resultado observado:
# - Qualidade excelente (Q35–Q38)
# - Leituras de aproximadamente 301 bp
# - Sem degradação significativa no final dos reads


################################################################################
# ETAPA 6 - QUANTIDADE DE READS POR AMOSTRA
################################################################################

sapply(fnFs, function(x) countFastq(x)$records)

# Resultado observado:
#
# BCO       = 5.297
# Algodão   = 167.470
# Milho     = 114.221
# Soja      = 116.910


################################################################################
# ETAPA 7 - AVALIAR O COMPRIMENTO DAS LEITURAS
################################################################################

# Ler uma amostra

fq <- readFastq(fnFs[1])

# Distribuição dos comprimentos

table(width(sread(fq)))

# Resumo para todas as amostras

for(i in 1:length(fnFs)) {
  
  cat("\n\n")
  cat(basename(fnFs[i]), "\n")
  
  print(
    summary(
      width(
        sread(
          readFastq(fnFs[i])
        )
      )
    )
  )
  
}

# Objetivo:
# Verificar presença de reads truncados e tamanho real das leituras.


################################################################################
# ETAPA 8 - CRIAR A PASTA DE ARQUIVOS FILTRADOS
################################################################################

filt_path <- file.path(path, "filtered")

dir.create(
  filt_path,
  showWarnings = FALSE
)


################################################################################
# ETAPA 9 - DEFINIR NOMES DOS FASTQ FILTRADOS
################################################################################

# Forward

filtFs <- file.path(
  filt_path,
  basename(fnFs)
)

# Reverse

filtRs <- file.path(
  filt_path,
  basename(fnRs)
)


################################################################################
# ETAPA 10 - FILTRAGEM DOS READS
################################################################################

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  maxN = 0,
  maxEE = c(2,2),
  truncQ = 2,
  minLen = 200,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

# Avaliar resultado

out

# Resultado observado:
#
# BCO      = 58 reads
# Algodão  = 137.983 reads
# Milho    = 98.915 reads
# Soja     = 112.979 reads


################################################################################
# ETAPA 11 - APRENDIZADO DO MODELO DE ERRO
################################################################################

# O DADA2 aprende as taxas de erro do sequenciamento para
# separar variantes verdadeiras de erros técnicos.

errF <- learnErrors(
  filtFs,
  multithread = TRUE
)

errR <- learnErrors(
  filtRs,
  multithread = TRUE
)

# Visualização

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Resultado observado:
# - Modelos convergiram corretamente.
# - Padrão típico esperado para dados Illumina.


################################################################################
# ETAPA 12 - INFERÊNCIA DAS ASVs
################################################################################

dadaFs <- dada(
  filtFs,
  err = errF,
  multithread = TRUE
)

dadaRs <- dada(
  filtRs,
  err = errR,
  multithread = TRUE
)

# Inspecionar

dadaFs
dadaRs

# Resultado observado:
#
# BCO      = 1 ASV
# Algodão  = 1421 ASVs
# Milho    = 1536 ASVs
# Soja     = 268 ASVs


################################################################################
# ETAPA 13 - INSPECIONAR AS PRIMEIRAS SEQUÊNCIAS
################################################################################

fq <- readFastq(fnFs[2])

head(as.character(sread(fq)), 10)

# Objetivo:
# Verificar presença dos primers e ausência dos overhangs Illumina.

# Resultado observado:
#
# - Overhangs Illumina removidos
# - Índices removidos
# - Primer rbcL ainda presente


################################################################################
# ETAPA 14 - MERGE DOS READS FORWARD E REVERSE
################################################################################

mergers <- mergePairs(
  dadaFs,
  filtFs,
  dadaRs,
  filtRs,
  verbose = TRUE
)

# Verificar processamento

length(mergers)

# Inspecionar a primeira amostra

head(mergers[[1]])

# Número de merges por amostra

sapply(mergers, nrow)


################################################################################
# ETAPA 15 - CRIAR A TABELA DE ASVs
################################################################################

seqtab <- makeSequenceTable(mergers)

# Dimensões

dim(seqtab)

# Comprimento das ASVs

table(nchar(colnames(seqtab)))

# Visualizar as sequências

colnames(seqtab)

nchar(colnames(seqtab))

# Resultado observado:
#
# 4 amostras x 1320 ASVs
#
# Comprimento variando de:
# ~200 bp até ~589 bp


################################################################################
# ETAPA 16 - ACOMPANHAR A RETENÇÃO DE READS
################################################################################

getN <- function(x) sum(getUniques(x))

track <- cbind(
  input = out[,1],
  filtered = out[,2],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN)
)

track

# Objetivo:
# Acompanhar quantos reads sobrevivem em cada etapa.


################################################################################
# ETAPA 17 - AVALIAR DISTRIBUIÇÃO DOS COMPRIMENTOS DAS ASVs
################################################################################

asv.lengths <- nchar(colnames(seqtab))

summary(asv.lengths)

hist(
  asv.lengths,
  breaks = 100,
  main = "Comprimento das ASVs merged"
)

# Resultado observado:
# Forte heterogeneidade de tamanhos.


################################################################################
# ETAPA 18 - INSPECIONAR AS ASVs MAIS ABUNDANTES
################################################################################

asv.sum <- colSums(seqtab)

head(
  sort(
    asv.sum,
    decreasing = TRUE
  ),
  20
)

top20 <- order(
  colSums(seqtab),
  decreasing = TRUE
)[1:20]

seqtab[, top20]

# Objetivo:
# Verificar quais amostras contribuem para as ASVs dominantes.


################################################################################
# ETAPA 19 - REMOÇÃO DE QUIMERAS
################################################################################

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE
)

dim(seqtab.nochim)

sum(seqtab.nochim) / sum(seqtab)

table(
  nchar(colnames(seqtab.nochim))
)

# Resultado observado:
#
# ASVs:
# 1320 -> 1306
#
# Retenção:
# 98,05%
#
# Pouca evidência de quimeras.


################################################################################
# ETAPA 20 - EXPORTAR AS 5 ASVs MAIS ABUNDANTES
################################################################################

asv.sum <- colSums(seqtab.nochim)

top5_idx <- order(
  asv.sum,
  decreasing = TRUE
)[1:5]

top5_seqs <- colnames(seqtab.nochim)[top5_idx]

top5_abund <- asv.sum[top5_idx]

top5_abund


################################################################################
# ETAPA 21 - CRIAR NOMES DAS ASVs
################################################################################

asv_names <- paste0(
  "ASV_",
  1:5,
  "_Abund_",
  top5_abund
)

asv_names


################################################################################
# ETAPA 22 - EXPORTAR FASTA
################################################################################

write.fasta(
  sequences = lapply(top5_seqs, s2c),
  names = asv_names,
  file.out = "Top5_ASVs_rbcL.fasta"
)

# Verificar

file.exists("Top5_ASVs_rbcL.fasta")


################################################################################
# ETAPA 23 - CONFERIR O FASTA
################################################################################

cat(
  readLines("Top5_ASVs_rbcL.fasta"),
  sep = "\n"
)


################################################################################
# ETAPA 24 - VALIDAR LEITURA FORWARD (R1) POR BLAST
################################################################################

fq <- readFastq(fnFs[2])

seq <- as.character(
  sread(fq)[1]
)

cat(seq)

# Resultado observado:
#
# rbcL cloroplastidial
# Cobertura = 100%
# Identidade ≈ 99,67%
# Gossypium spp.
#
# Conclusão:
# Amplificação correta do marcador rbcL.


################################################################################
# ETAPA 25 - VALIDAR LEITURA REVERSE (R2) POR BLAST
################################################################################

fqR <- readFastq(fnRs[2])

seqR <- as.character(
  sread(fqR)[1]
)

cat(seqR)

# Resultado observado:
#
# rbcL cloroplastidial
# Cobertura = 100%
# Identidade ≈ 99,3–99,7%
#
# Conclusão:
# O primer reverse também amplificou corretamente o gene rbcL.


################################################################################
# CONCLUSÕES DA ANÁLISE
################################################################################

# ✅ Sequenciamento com excelente qualidade.
#
# ✅ Primers amplificaram corretamente o marcador rbcL.
#
# ✅ R1 e R2 apresentaram alinhamento ao gene rbcL com ~99,7% de identidade.
#
# ✅ Controle branco praticamente limpo após filtragem.
#
# ✅ Baixa taxa de quimeras (~2%).
#
# ⚠️ A etapa de merge apresentou eficiência limitada.
#
# ⚠️ O amplicon encontra-se próximo do limite operacional do sequenciamento
#    paired-end 2x300 bp.
#
# ⚠️ Recomenda-se avaliar também uma pipeline single-end utilizando apenas
#    os reads Forward (R1), que apresentaram excelente qualidade e BLAST
#    consistente para rbcL.
################################################################################
