#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FastaGen CLI - Análise genômica em linha de comando
Uso: python fastagen_cli.py -i amostra.fasta -r referencia.fasta [opcoes]
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import pandas as pd

# ============================================================
# FUNÇÕES DE ANÁLISE
# ============================================================
def calcular_gc(seq):
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100 if len(seq) > 0 else 0.0

def alinhar_sequencias(seq_query, seq_ref):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0
    aligner.extend_gap_score = 0
    alns = aligner.align(seq_query, seq_ref)
    try:
        primeiro = next(alns)
        return str(primeiro[0]), str(primeiro[1])
    except StopIteration:
        return None, None

def detectar_mutacoes(seq_aln_query, seq_aln_ref):
    mutacoes = []
    for i, (q, r) in enumerate(zip(seq_aln_query, seq_aln_ref)):
        if q != r and q != '-' and r != '-':
            mutacoes.append({"Posição": i, "Ref": r, "Alt": q})
    return mutacoes

def impacto_proteico(seq_aln_query, seq_aln_ref):
    try:
        prot_query = Seq(seq_aln_query.replace("-", "")).translate()
        prot_ref   = Seq(seq_aln_ref.replace("-", "")).translate()
        return "Alterado" if prot_query != prot_ref else "Silencioso"
    except Exception:
        return "Indefinido"

def salvar_excel(df_resumo, df_mutacoes, arquivo):
    with pd.ExcelWriter(arquivo, engine='openpyxl') as writer:
        df_resumo.to_excel(writer, index=False, sheet_name='Resumo')
        if not df_mutacoes.empty:
            df_mutacoes.to_excel(writer, index=False, sheet_name='Mutacoes')
    print(f" Excel salvo: {arquivo}")

def salvar_vcf(mutacoes, arquivo):
    with open(arquivo, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for m in mutacoes:
            pos = m["Posição"] + 1
            f.write(f"chr1\t{pos}\t.\t{m['Ref']}\t{m['Alt']}\t.\tPASS\t.\n")
    print(f" VCF salvo: {arquivo}")

def processar_analise(amostra_path, ref_path, formato, gc_min, gc_max, prefixo):
    # Carregar referência
    print(f"Carregando referência: {ref_path}")
    ref_record = next(SeqIO.parse(ref_path, "fasta"))
    ref_seq = str(ref_record.seq)

    # Carregar amostra
    records = list(SeqIO.parse(amostra_path, formato))
    total = len(records)
    print(f"Processando {total} sequências de {amostra_path}...")

    resultados = []
    todas_mutacoes = []

    for i, record in enumerate(records, 1):
        seq = str(record.seq)
        gc = calcular_gc(seq)

        # Alinhamento
        aln_query, aln_ref = alinhar_sequencias(seq, ref_seq)
        if aln_query is None:
            print(f"  [{i}/{total}] Falha no alinhamento de {record.id}. Ignorado.")
            continue

        mutacoes = detectar_mutacoes(aln_query, aln_ref)
        impacto = impacto_proteico(aln_query, aln_ref)

        if gc < gc_min or gc > gc_max:
            status = "Suspeito"
        else:
            status = "Normal"

        # Armazena mutações
        for m in mutacoes:
            m["ID"] = record.id
            todas_mutacoes.append(m)

        resultados.append({
            "ID": record.id,
            "Tamanho (bp)": len(seq),
            "GC (%)": round(gc, 2),
            "Mutações": len(mutacoes),
            "Impacto Proteico": impacto,
            "Status": status
        })

        print(f"  [{i}/{total}] {record.id}: GC={gc:.1f}%, mutações={len(mutacoes)}, status={status}")

    # Criar DataFrames
    df_resumo = pd.DataFrame(resultados)
    df_mutacoes = pd.DataFrame(todas_mutacoes)

    # Salvar arquivos
    salvar_excel(df_resumo, df_mutacoes, f"{prefixo}_relatorio.xlsx")
    if not df_mutacoes.empty:
        salvar_vcf(todas_mutacoes, f"{prefixo}_variantes.vcf")
    else:
        print(" Nenhuma mutação encontrada. VCF não gerado.")

    # Estatísticas finais
    print("\n" + "="*50)
    print("RESUMO FINAL")
    print("="*50)
    print(f"Sequências processadas : {len(df_resumo)}")
    print(f"Mutações totais        : {len(df_mutacoes)}")
    print(f"GC médio               : {df_resumo['GC (%)'].mean():.2f}%")
    print(f"Sequências suspeitas   : {df_resumo[df_resumo['Status']=='Suspeito'].shape[0]}")
    print(f"Impacto proteico alterado: {df_resumo[df_resumo['Impacto Proteico']=='Alterado'].shape[0]}")
    print("="*50)

def main():
    parser = argparse.ArgumentParser(
        description="FastaGen CLI - Análise genômica (GC, mutações, impacto proteico)",
        epilog="Exemplo: python fastagen_cli.py -i amostra.fastq -r referencia.fasta -o resultado"
    )
    parser.add_argument("-i", "--input", required=True, help="Arquivo de amostra (FASTA/FASTQ)")
    parser.add_argument("-r", "--reference", required=True, help="Arquivo de referência (FASTA)")
    parser.add_argument("-o", "--output", default="resultado", help="Prefixo dos arquivos de saída (padrão: resultado)")
    parser.add_argument("-f", "--format", choices=["fasta", "fastq"], default=None,
                        help="Formato do arquivo de entrada (detectado automaticamente se não informado)")
    parser.add_argument("--gc-min", type=int, default=40, help="GC mínimo para classificar como Normal (padrão: 40)")
    parser.add_argument("--gc-max", type=int, default=60, help="GC máximo para classificar como Normal (padrão: 60)")

    args = parser.parse_args()

    # Verificar existência dos arquivos
    if not Path(args.input).exists():
        sys.exit(f"Erro: arquivo de entrada '{args.input}' não encontrado.")
    if not Path(args.reference).exists():
        sys.exit(f"Erro: arquivo de referência '{args.reference}' não encontrado.")

    # Detectar formato se não informado
    formato = args.format
    if formato is None:
        if args.input.lower().endswith(('.fastq', '.fq')):
            formato = "fastq"
        elif args.input.lower().endswith(('.fasta', '.fa', '.fna')):
            formato = "fasta"
        else:
            sys.exit("Erro: não foi possível detectar o formato. Use -f fasta ou fastq.")

    # Validar limites de GC
    if args.gc_min < 0 or args.gc_max > 100 or args.gc_min >= args.gc_max:
        sys.exit("Erro: limites de GC inválidos. Use --gc-min e --gc-max com 0-100 e gc_min < gc_max.")

    processar_analise(args.input, args.reference, formato, args.gc_min, args.gc_max, args.output)

if __name__ == "__main__":
    main()