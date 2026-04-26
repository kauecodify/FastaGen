#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FastaGen
Análise genômica: GC, mutações, impacto proteico, exportação Excel/VCF
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from threading import Thread
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import os
import sys

# ============================================================
# FUNÇÕES DE ANÁLISE (Núcleo)
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

def processar_analise(amostra_path, ref_path, formato, gc_min, gc_max, callback_progresso, callback_log):
    """Executa a análise e retorna (df_resumo, df_mutacoes)."""
    try:
        ref_record = next(SeqIO.parse(ref_path, "fasta"))
        ref_seq = str(ref_record.seq)

        records = list(SeqIO.parse(amostra_path, formato))
        total = len(records)
        if total == 0:
            raise Exception("Nenhuma sequência encontrada no arquivo de amostra.")

        resultados = []
        todas_mutacoes = []

        for i, record in enumerate(records, 1):
            seq = str(record.seq)
            gc = calcular_gc(seq)

            aln_query, aln_ref = alinhar_sequencias(seq, ref_seq)
            if aln_query is None:
                callback_log(f"⚠️ Falha no alinhamento de {record.id}. Ignorado.")
                continue

            mutacoes = detectar_mutacoes(aln_query, aln_ref)
            impacto = impacto_proteico(aln_query, aln_ref)

            if gc < gc_min or gc > gc_max:
                status = "Suspeito"
            else:
                status = "Normal"

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

            callback_progresso(i, total, record.id)

        df_resumo = pd.DataFrame(resultados)
        df_mutacoes = pd.DataFrame(todas_mutacoes)
        return df_resumo, df_mutacoes
    except Exception as e:
        raise e

# ============================================================
# CLASSE DA INTERFACE TKINTER (MODERNA)
# ============================================================
class FastaGenApp:
    def __init__(self, root):
        self.root = root
        self.root.title("FastaGen MVP - Análise Genômica")
        self.root.geometry("1200x750")
        self.root.minsize(1000, 600)
        self.root.configure(bg='#f5f5f5')

        # Variáveis de controle
        self.amostra_path = tk.StringVar()
        self.referencia_path = tk.StringVar()
        self.formato = tk.StringVar(value="auto")
        self.gc_min = tk.IntVar(value=40)
        self.gc_max = tk.IntVar(value=60)
        self.prefixo = tk.StringVar(value="resultado")

        # Configurar estilos ttk
        self.configurar_estilos()

        # Criar layout principal
        self.criar_layout()

        # Armazenar resultados
        self.df_resumo = None
        self.df_mutacoes = None

    def configurar_estilos(self):
        style = ttk.Style()
        style.theme_use('clam')

        # Cores principais: roxo (#6A0DAD, #8B5CF6, #D8B4FE)
        style.configure('TFrame', background='#f5f5f5')
        style.configure('TLabel', background='#f5f5f5', font=('Segoe UI', 10))
        style.configure('TLabelframe', background='#f5f5f5', foreground='#6A0DAD', font=('Segoe UI', 10, 'bold'))
        style.configure('TLabelframe.Label', background='#f5f5f5', foreground='#6A0DAD', font=('Segoe UI', 10, 'bold'))

        # Botões
        style.configure('Roxo.TButton', background='#6A0DAD', foreground='white', font=('Segoe UI', 10, 'bold'), borderwidth=0, focuscolor='none')
        style.map('Roxo.TButton', background=[('active', '#8B5CF6'), ('disabled', '#cccccc')])

        # Entry
        style.configure('TEntry', fieldbackground='white', borderwidth=1, relief='solid')
        style.configure('TCombobox', fieldbackground='white')

        # Progressbar
        style.configure('Roxo.Horizontal.TProgressbar', background='#6A0DAD', troughcolor='#e0e0e0', thickness=12)

        # Treeview
        style.configure('Treeview', background='white', foreground='#333', fieldbackground='white', rowheight=25, font=('Segoe UI', 9))
        style.map('Treeview', background=[('selected', '#D8B4FE')])
        style.configure('Treeview.Heading', background='#6A0DAD', foreground='white', font=('Segoe UI', 10, 'bold'), relief='flat')
        style.map('Treeview.Heading', background=[('active', '#8B5CF6')])

    def criar_layout(self):
        # Frame principal (container)
        main_container = ttk.Frame(self.root)
        main_container.pack(fill=tk.BOTH, expand=True, padx=15, pady=15)

        # ===== CABEÇALHO (logo e título) =====
        header = ttk.Frame(main_container)
        header.pack(fill=tk.X, pady=(0, 15))

        titulo = tk.Label(header, text="FastaGen", font=('Segoe UI', 28, 'bold'), fg='#6A0DAD', bg='#f5f5f5')
        titulo.pack(side=tk.LEFT)

        subtitulo = tk.Label(header, text="Análise Genômica", font=('Segoe UI', 12), fg='#666', bg='#f5f5f5')
        subtitulo.pack(side=tk.LEFT, padx=(10, 0), pady=(10, 0))

        # Linha divisória
        separator = ttk.Separator(main_container, orient='horizontal')
        separator.pack(fill=tk.X, pady=5)

        # ===== CORPO PRINCIPAL (2 colunas) =====
        body = ttk.Frame(main_container)
        body.pack(fill=tk.BOTH, expand=True, pady=10)

        # Coluna esquerda (configurações)
        left_panel = ttk.Frame(body, width=350)
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 15))
        left_panel.pack_propagate(False)

        # Coluna direita (resultados)
        right_panel = ttk.Frame(body)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # ----- Painel esquerdo: Arquivos -----
        file_frame = ttk.LabelFrame(left_panel, text="📂 Arquivos de entrada", padding=10)
        file_frame.pack(fill=tk.X, pady=(0, 15))

        # Amostra
        ttk.Label(file_frame, text="Amostra (FASTA/FASTQ):").grid(row=0, column=0, sticky=tk.W, pady=5)
        entry_amostra = ttk.Entry(file_frame, textvariable=self.amostra_path, width=30)
        entry_amostra.grid(row=0, column=1, padx=5, sticky=tk.EW)
        btn_amostra = ttk.Button(file_frame, text="📂", command=self.selecionar_amostra, width=3)
        btn_amostra.grid(row=0, column=2, padx=2)

        # Referência
        ttk.Label(file_frame, text="Referência (FASTA):").grid(row=1, column=0, sticky=tk.W, pady=5)
        entry_ref = ttk.Entry(file_frame, textvariable=self.referencia_path, width=30)
        entry_ref.grid(row=1, column=1, padx=5, sticky=tk.EW)
        btn_ref = ttk.Button(file_frame, text="📂", command=self.selecionar_referencia, width=3)
        btn_ref.grid(row=1, column=2, padx=2)

        file_frame.columnconfigure(1, weight=1)

        # ----- Painel esquerdo: Opções -----
        opts_frame = ttk.LabelFrame(left_panel, text="⚙️ Parâmetros", padding=10)
        opts_frame.pack(fill=tk.X, pady=(0, 15))

        # Formato
        ttk.Label(opts_frame, text="Formato:").grid(row=0, column=0, sticky=tk.W, pady=5)
        formato_combo = ttk.Combobox(opts_frame, textvariable=self.formato, values=["auto", "fasta", "fastq"], state="readonly", width=10)
        formato_combo.grid(row=0, column=1, sticky=tk.W, padx=5)

        # GC min/max
        ttk.Label(opts_frame, text="GC min (Normal):").grid(row=1, column=0, sticky=tk.W, pady=5)
        spin_min = ttk.Spinbox(opts_frame, from_=0, to=100, textvariable=self.gc_min, width=8)
        spin_min.grid(row=1, column=1, sticky=tk.W, padx=5)

        ttk.Label(opts_frame, text="GC max (Normal):").grid(row=2, column=0, sticky=tk.W, pady=5)
        spin_max = ttk.Spinbox(opts_frame, from_=0, to=100, textvariable=self.gc_max, width=8)
        spin_max.grid(row=2, column=1, sticky=tk.W, padx=5)

        # Prefixo
        ttk.Label(opts_frame, text="Prefixo saída:").grid(row=3, column=0, sticky=tk.W, pady=5)
        entry_prefixo = ttk.Entry(opts_frame, textvariable=self.prefixo, width=15)
        entry_prefixo.grid(row=3, column=1, sticky=tk.W, padx=5)

        # Botão executar
        self.btn_executar = ttk.Button(left_panel, text="EXECUTAR ANÁLISE", style='Roxo.TButton', command=self.iniciar_analise)
        self.btn_executar.pack(pady=10, fill=tk.X)

        # Barra de progresso e status
        self.progresso = ttk.Progressbar(left_panel, style='Roxo.Horizontal.TProgressbar', mode='determinate')
        self.progresso.pack(fill=tk.X, pady=5)

        self.lbl_status = ttk.Label(left_panel, text="Pronto para análise", font=('Segoe UI', 9, 'italic'))
        self.lbl_status.pack()

        # ----- Painel direito: Notebook com abas -----
        notebook = ttk.Notebook(right_panel)
        notebook.pack(fill=tk.BOTH, expand=True)

        # Aba 1: Resultados (tabela)
        tab_results = ttk.Frame(notebook)
        notebook.add(tab_results, text="📊 Resultados")

        # Treeview com scroll
        tree_frame = ttk.Frame(tab_results)
        tree_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        self.tree = ttk.Treeview(tree_frame, columns=("ID", "Tamanho", "GC", "Mutações", "Impacto", "Status"), show="headings", height=18)
        self.tree.heading("ID", text="ID")
        self.tree.heading("Tamanho", text="Tamanho (bp)")
        self.tree.heading("GC", text="GC (%)")
        self.tree.heading("Mutações", text="Mutações")
        self.tree.heading("Impacto", text="Impacto")
        self.tree.heading("Status", text="Status")
        self.tree.column("ID", width=250)
        self.tree.column("Tamanho", width=100)
        self.tree.column("GC", width=80)
        self.tree.column("Mutações", width=80)
        self.tree.column("Impacto", width=120)
        self.tree.column("Status", width=100)

        scroll_y = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL, command=self.tree.yview)
        scroll_x = ttk.Scrollbar(tree_frame, orient=tk.HORIZONTAL, command=self.tree.xview)
        self.tree.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
        self.tree.grid(row=0, column=0, sticky='nsew')
        scroll_y.grid(row=0, column=1, sticky='ns')
        scroll_x.grid(row=1, column=0, sticky='ew')
        tree_frame.grid_rowconfigure(0, weight=1)
        tree_frame.grid_columnconfigure(0, weight=1)

        # Aba 2: Log da análise
        tab_log = ttk.Frame(notebook)
        notebook.add(tab_log, text="📝 Log")

        self.log_text = scrolledtext.ScrolledText(tab_log, wrap=tk.WORD, font=('Consolas', 9), bg='#fafafa', fg='#333')
        self.log_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Aba 3: Relatório de mutações (opcional)
        tab_muts = ttk.Frame(notebook)
        notebook.add(tab_muts, text="🧬 Mutações detalhadas")

        self.tree_muts = ttk.Treeview(tab_muts, columns=("ID", "Posição", "Ref", "Alt"), show="headings")
        self.tree_muts.heading("ID", text="ID")
        self.tree_muts.heading("Posição", text="Posição (0‑based)")
        self.tree_muts.heading("Ref", text="Referência")
        self.tree_muts.heading("Alt", text="Amostra")
        self.tree_muts.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Botões de exportação
        export_frame = ttk.Frame(right_panel)
        export_frame.pack(fill=tk.X, pady=(10, 0))

        self.btn_excel = ttk.Button(export_frame, text="📎 Salvar Excel", command=self.salvar_excel, state=tk.DISABLED, style='Roxo.TButton')
        self.btn_excel.pack(side=tk.LEFT, padx=5)

        self.btn_vcf = ttk.Button(export_frame, text="📎 Salvar VCF", command=self.salvar_vcf, state=tk.DISABLED, style='Roxo.TButton')
        self.btn_vcf.pack(side=tk.LEFT, padx=5)

    # ========== MÉTODOS AUXILIARES ==========
    def selecionar_amostra(self):
        path = filedialog.askopenfilename(title="Arquivo de amostra", filetypes=[("Sequências", "*.fasta *.fa *.fastq"), ("Todos", "*.*")])
        if path:
            self.amostra_path.set(path)

    def selecionar_referencia(self):
        path = filedialog.askopenfilename(title="Arquivo de referência", filetypes=[("FASTA", "*.fasta *.fa"), ("Todos", "*.*")])
        if path:
            self.referencia_path.set(path)

    def adicionar_log(self, msg):
        self.log_text.insert(tk.END, msg + "\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    def atualizar_progresso(self, atual, total, seq_id):
        percent = int((atual / total) * 100)
        self.progresso["value"] = percent
        self.lbl_status["text"] = f"Processando: {atual}/{total} - {seq_id}"
        self.root.update_idletasks()

    def iniciar_analise(self):
        if not self.amostra_path.get() or not self.referencia_path.get():
            messagebox.showerror("Erro", "Selecione os arquivos de amostra e referência.")
            return

        if self.gc_min.get() >= self.gc_max.get():
            messagebox.showerror("Erro", "GC mínimo deve ser menor que GC máximo.")
            return

        self.btn_executar.config(state=tk.DISABLED)
        self.btn_excel.config(state=tk.DISABLED)
        self.btn_vcf.config(state=tk.DISABLED)
        self.progresso["value"] = 0
        self.log_text.delete(1.0, tk.END)
        self.tree.delete(*self.tree.get_children())
        self.tree_muts.delete(*self.tree_muts.get_children())

        formato = self.formato.get()
        if formato == "auto":
            if self.amostra_path.get().lower().endswith(".fastq"):
                formato = "fastq"
            else:
                formato = "fasta"

        def thread_analise():
            try:
                df_resumo, df_mutacoes = processar_analise(
                    self.amostra_path.get(),
                    self.referencia_path.get(),
                    formato,
                    self.gc_min.get(),
                    self.gc_max.get(),
                    self.atualizar_progresso,
                    self.adicionar_log
                )
                self.df_resumo = df_resumo
                self.df_mutacoes = df_mutacoes

                # Popula tabela resumo
                for _, row in df_resumo.iterrows():
                    self.tree.insert("", tk.END, values=(
                        row["ID"], row["Tamanho (bp)"], row["GC (%)"],
                        row["Mutações"], row["Impacto Proteico"], row["Status"]
                    ))

                # Popula tabela de mutações
                if not df_mutacoes.empty:
                    for _, row in df_mutacoes.iterrows():
                        self.tree_muts.insert("", tk.END, values=(row["ID"], row["Posição"], row["Ref"], row["Alt"]))

                self.adicionar_log(f"\n✅ Análise concluída! {len(df_resumo)} sequências processadas.")
                self.lbl_status["text"] = "Análise concluída."
                self.btn_excel.config(state=tk.NORMAL)
                if not df_mutacoes.empty:
                    self.btn_vcf.config(state=tk.NORMAL)
                else:
                    self.adicionar_log("ℹ️ Nenhuma mutação encontrada. VCF não disponível.")
            except Exception as e:
                self.adicionar_log(f"❌ ERRO: {str(e)}")
                messagebox.showerror("Erro na análise", str(e))
            finally:
                self.btn_executar.config(state=tk.NORMAL)

        Thread(target=thread_analise, daemon=True).start()

    def salvar_excel(self):
        if self.df_resumo is None:
            return
        path = filedialog.asksaveasfilename(defaultextension=".xlsx", filetypes=[("Excel", "*.xlsx")],
                                             initialfile=f"{self.prefixo.get()}_relatorio.xlsx")
        if path:
            with pd.ExcelWriter(path, engine='openpyxl') as writer:
                self.df_resumo.to_excel(writer, index=False, sheet_name='Resumo')
                if self.df_mutacoes is not None and not self.df_mutacoes.empty:
                    self.df_mutacoes.to_excel(writer, index=False, sheet_name='Mutacoes')
            self.adicionar_log(f"✅ Excel salvo: {path}")

    def salvar_vcf(self):
        if self.df_mutacoes is None or self.df_mutacoes.empty:
            return
        path = filedialog.asksaveasfilename(defaultextension=".vcf", filetypes=[("VCF", "*.vcf")],
                                             initialfile=f"{self.prefixo.get()}_variantes.vcf")
        if path:
            with open(path, 'w') as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                for _, row in self.df_mutacoes.iterrows():
                    pos = row["Posição"] + 1
                    f.write(f"chr1\t{pos}\t.\t{row['Ref']}\t{row['Alt']}\t.\tPASS\t.\n")
            self.adicionar_log(f"✅ VCF salvo: {path}")

# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    root = tk.Tk()
    app = FastaGenApp(root)
    root.mainloop()