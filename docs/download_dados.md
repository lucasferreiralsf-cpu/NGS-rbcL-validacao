# Download de dados via BaseSpace CLI

Este documento descreve o processo utilizado para instalar e configurar o **BaseSpace CLI** no macOS e realizar o download dos arquivos FASTQ do run **INTECSO_18JUNHO2026**.

---

## 1. Instalação do CLI Illumina

mkdir -p $HOME/bin
curl -o $HOME/bin/bs https://launch.basespace.illumina.com/CLI/latest/amd64-osx/bs
chmod u+x $HOME/bin/bs
echo 'export PATH=$HOME/bin:$PATH' >> ~/.zshrc
source ~/.zshrc
bs --version

#Saida esperada
BaseSpaceCLI 1.7.0 -- built on 2026-01-20 at 17:10

## 2. Autenticação
bs auth

## 3. Listar projetos disponiveis
bs list projects
identificar o id do projeto de interesse 
ex.: | INTECSO_18JUNHO2026_2026-06-18T17_17_51_6e9e5db | 507459959 | 930000653 |
O PROJECT_ID do run INTECSO_18JUNHO2026 é 507459959.

## 4. Download dos dados 
bs download project --id 507459959 --output ~/NGS-rbcL-validacao/data/
Os arquivos .fastq.gz de cada biosample (ALFACE, ARROZ, BATATA_DOCE etc.) são baixados diretamente para a pasta data/.

## 5. Conferencia de arquivos
ls -lh ~/NGS-rbcL-validacao/data/
→ Lista os arquivos FASTQ baixados.

## 6. Controle de versão
Os arquivos .fastq.gz não são versionados no GitHub.
O .gitignore contém regras para ignorar arquivos grandes de NGS:

*.fastq
*.fastq.gz
*.bam
*.sam
*.sra
*.fq

