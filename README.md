# FastaGen

```markdown
# FastaGen - Ferramenta de Análise Genômica

## Descrição

FastaGen é uma ferramenta leve e amigável para análise de sequências genômicas. Ela compara um conjunto de sequências de consulta (FASTA/FASTQ) com um genoma de referência (FASTA) para:

- Calcular o conteúdo GC (%)
- Detectar variantes de nucleotídeo único (SNVs)
- Predizer o impacto na proteína (silencioso ou alterado)
- Gerar relatórios resumidos (Excel) e arquivos no formato VCF (Variant Call Format)

A ferramenta está disponível em duas versões:
- **Interface de linha de comando (CLI)** – adequada para processamento em lote e integração em pipelines.
- **Interface gráfica (GUI)** – construída com Tkinter, oferecendo layout moderno com logs em tempo real, barra de progresso e abas de resultados.

## Funcionalidades

- Detecção automática do formato (FASTA/FASTQ)
- Alinhamento pareado global usando o PairwiseAligner do Biopython
- Detecção de substituições de bases (ignorando indels para simplificar)
- Classificação das amostras com base no GC (normal/suspeito com limites ajustáveis)
- Tradução das sequências de DNA para proteína e comparação com a referência
- Exportação dos resultados para Excel (tabela resumo e tabela detalhada de mutações)
- Exportação das variantes para o formato VCF v4.2
- Interface gráfica com multithreading (a interface não congela durante a análise)

## Requisitos

- Python 3.7 ou superior
- Pacotes Python necessários:
  - biopython
  - pandas
  - openpyxl

Opcional (para a versão GUI):
- Tkinter (geralmente incluso no Python no Windows; no Linux pode ser necessário `python3-tk`)

## Instalação

1. Clone ou baixe o repositório.

2. Instale as dependências usando pip:

   ```bash
   pip install biopython pandas openpyxl
   ```

3. (Opcional) Para a versão GUI no Linux:
   ```bash
   sudo apt-get install python3-tk
   ```

## Uso

### Interface de Linha de Comando (CLI)

Execute o script `fastagen.py` no terminal.

```bash
python fastagen.py -i amostra.fasta -r referencia.fasta [opcoes]
```

#### Argumentos obrigatórios

| Argumento | Descrição |
|-----------|-----------|
| `-i, --input` | Caminho para as sequências de consulta (FASTA/FASTQ) |
| `-r, --reference` | Caminho para o genoma de referência (FASTA) |

#### Argumentos opcionais

| Argumento | Descrição | Padrão |
|-----------|-----------|--------|
| `-o, --output` | Prefixo dos arquivos de saída | `resultado` |
| `-f, --format` | Formato de entrada: `fasta`, `fastq` ou `auto` | `auto` |
| `--gc-min` | Percentual GC mínimo para classificação "Normal" | 40 |
| `--gc-max` | Percentual GC máximo para classificação "Normal" | 60 |

#### Exemplo

```bash
python fastagen.py -i minhas_sequencias.fastq -r referencia.fasta -o minha_analise --gc-min 35 --gc-max 65
```

### Interface Gráfica (GUI)

Execute o script `fastagen_ui.py`:

```bash
python fastagen_ui.py
```

A janela principal contém:

- **Painel esquerdo** – seleção de arquivos, ajuste de parâmetros e botão de execução.
- **Painel direito** – três abas:
  - Resultados – tabela resumo de todas as sequências processadas.
  - Log – mensagens de processamento em tempo real.
  - Mutações detalhadas – lista de todas as substituições detectadas.

Após a análise, clique nos botões "Salvar Excel" ou "Salvar VCF" para exportar os resultados.

## Entrada/Saída

### Arquivos de entrada

- **Arquivo de consulta**: formato FASTA ou FASTQ. Cada registro representa uma sequência amostra.
- **Arquivo de referência**: formato FASTA. A primeira sequência do arquivo é usada como referência.

### Arquivos de saída

- **Relatório Excel** (`<prefixo>_relatorio.xlsx`):
  - Aba "Resumo": uma linha por sequência de consulta com as colunas:
    - ID, Tamanho (pb), GC (%), Mutações, Impacto na Proteína (Silencioso/Alterado/Indefinido), Status (Normal/Suspeito)
  - Aba "Mutacoes": lista detalhada das mutações (somente se encontradas) com colunas:
    - ID, Posição (base 0), Ref, Alt

- **Arquivo VCF** (`<prefixo>_variantes.vcf`): VCF versão 4.2 com uma entrada por substituição. As posições são convertidas para coordenadas base 1.

## Como Funciona

1. O genoma de referência é carregado.
2. Para cada sequência de consulta:
   - O conteúdo GC é calculado.
   - A sequência de consulta é alinhada globalmente à referência usando um esquema de pontuação simples (match=1, mismatch=0, gap=0).
   - Substituições são identificadas (posições onde as bases diferem e nenhuma é um gap).
   - As sequências alinhadas (com gaps removidos) são traduzidas para proteína e comparadas para determinar o impacto na proteína.
   - Com base no percentual GC e nos limites definidos pelo usuário, a sequência é classificada como "Normal" ou "Suspeito".
3. Uma tabela resumo e uma tabela de mutações são geradas e podem ser exportadas.

## Limitações

- Alinhamento global sem penalidades para gaps: útil para comparar sequências de comprimento semelhante (ex.: genomas virais, amplicons). Não recomendado para sequências muito divergentes ou variantes estruturais grandes.
- Apenas substituições são relatadas; inserções e deleções são ignoradas na versão atual.
- A saída VCF usa um nome de cromossomo fixo "chr1". Para referências com múltiplos cromossomos, a interpretação da posição deve ser feita com cuidado.

## Licença

Este projeto está sob a licença mit.

## Autor

Desenvolvido como parte de um pipeline MVP de bioinformática.
```
by k.
