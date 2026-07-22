## 1. Instalação dos pacotes (executar apenas uma vez)
# Instala o gerenciador do Bioconductor
install.packages("BiocManager")

# Instala o DADA2
BiocManager::install("dada2")

# Instala pacotes auxiliares
install.packages(c(
  "ggplot2",
  "tidyverse"
))

# Instala o pacote para leitura de FASTQ
BiocManager::install("ShortRead")

## 2. Carregar os pacotes
library(dada2)
library(ShortRead)
library(ggplot2)

## 3. Definir o diretório de trabalho
# Pasta contendo os arquivos FASTQ

path <- "C:/Users/lucas.ferreira/OneDrive - Intecso Soluções e Inovações em Agronegócios LTDA/Área de Trabalho/NGS_rbcL" #personalizar o caminho dos arquivos nessa etapa

# Verificar se a pasta existe
dir.exists(path)

## 4. Verificar os arquivos disponíveis
list.files(path)

## 5. Separar leituras Forward (R1) e Reverse (R2)
# Leituras forward
fnFs <- sort(
  list.files(
    path,
    pattern = "_R1_001.fastq.gz$",
    full.names = TRUE
  )
)

# Leituras reverse
fnRs <- sort(
  list.files(
    path,
    pattern = "_R2_001.fastq.gz$",
    full.names = TRUE
  )
)

# Visualizar os arquivos encontrados:
fnFs
fnRs

## 6. Avaliação da qualidade das leituras
# Primeira amostra (Forward)

plotQualityProfile(fnFs[1])

# Primeira amostra (Reverse)

plotQualityProfile(fnRs[1])

## 7. Quantidade de reads por amostra

sapply(fnFs, function(x) countFastq(x)$records)

## 8. Avaliar comprimento das sequências
# Ler uma amostra:

fq <- readFastq(fnFs[1])

# Contar os comprimentos observados:

table(width(sread(fq)))

## 9. Resumo dos comprimentos para todas as amostras

for(i in 1:length(fnFs)){
  
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

## 10. Criar pasta para os arquivos filtrados

filt_path <- file.path(path, "filtered")

dir.create(
  filt_path,
  showWarnings = FALSE
)

## 11. Criar nomes dos arquivos filtrados
# Forward filtrados

filtFs <- file.path(
  filt_path,
  basename(fnFs)
)

# Reverse filtrados

filtRs <- file.path(
  filt_path,
  basename(fnRs)
)


## 12. Filtragem de reads
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

# Avaliar a filtragem

out

## 13. Aprendizado do modelo de erro (DADA2)

errF <- learnErrors(
  filtFs,
  multithread = TRUE
)

errR <- learnErrors(
  filtRs,
  multithread = TRUE
)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
