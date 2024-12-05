import csv
import os
import time
from collections import defaultdict

# codon to amino acid mapping
codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# function to convert codon to amino acid
def codon_to_amino_acid(codon):
    return codon_table.get(codon.upper(), None)

# MaxHeap implementation
class MaxHeap:
    def __init__(self):
        self.heap = []

    def insert(self, key, value):
        self.heap.append((value, key))
        self._heapify_up(len(self.heap) - 1)

    def extract_max(self):
        if not self.heap:
            return None
        max_item = self.heap[0]
        last_item = self.heap.pop()
        if self.heap:
            self.heap[0] = last_item
            self._heapify_down(0)
        return max_item

    def _heapify_up(self, index):
        parent = (index - 1) // 2
        if parent >= 0 and self.heap[parent][0] < self.heap[index][0]:
            self.heap[parent], self.heap[index] = self.heap[index], self.heap[parent]
            self._heapify_up(parent)

    def _heapify_down(self, index):
        largest = index
        size = len(self.heap)
        left = 2 * index + 1
        right = 2 * index + 2
        if left < size and self.heap[left][0] > self.heap[largest][0]:
            largest = left
        if right < size and self.heap[right][0] > self.heap[largest][0]:
            largest = right
        if largest != index:
            self.heap[largest], self.heap[index] = self.heap[index], self.heap[largest]
            self._heapify_down(largest)

    def is_empty(self):
        return len(self.heap) == 0

class Transcript:
    def __init__(self, gene_name):
        self.gene_name = gene_name
        self.amino_acid_codons = {}
        self.total_amino_acid_counts = {}
        self.heaps = {}

    def add_sequence(self, sequence):
        codons = [sequence[i:i+3] for i in range(0, len(sequence) - 2, 3)]
        for codon in codons:
            amino_acid = codon_to_amino_acid(codon)
            if amino_acid is None:
                continue
            if amino_acid not in self.amino_acid_codons:
                self.amino_acid_codons[amino_acid] = {}
                self.total_amino_acid_counts[amino_acid] = 0
            self.amino_acid_codons[amino_acid][codon] = self.amino_acid_codons[amino_acid].get(codon, 0) + 1
            self.total_amino_acid_counts[amino_acid] += 1

    def calculate_usage_rates(self):
        for amino_acid in self.amino_acid_codons:
            total_count = self.total_amino_acid_counts[amino_acid]
            heap = MaxHeap()
            for codon in self.amino_acid_codons[amino_acid]:
                count = self.amino_acid_codons[amino_acid][codon]
                usage_rate = count / total_count
                heap.insert(codon, usage_rate)
            self.heaps[amino_acid] = heap

    def get_optimal_codons(self):
        optimal_codons = {}
        for amino_acid in self.heaps:
            heap = self.heaps[amino_acid]
            if not heap.is_empty():
                usage_rate, codon = heap.extract_max()
                optimal_codons[amino_acid] = (codon, usage_rate)
        return optimal_codons

def process_file(filename):
    start_time = time.time()
    transcripts = {}
    transcript_counts = defaultdict(int)

    with open(filename, 'r', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        for row in reader:
            gene_name = row['gene_name']
            sequence = row['sequence']
            transcript_counts[gene_name] += 1
            if gene_name not in transcripts:
                transcripts[gene_name] = Transcript(gene_name)
            transcripts[gene_name].add_sequence(sequence)

    for transcript in transcripts.values():
        transcript.calculate_usage_rates()

    output_data = []
    for gene_name, transcript in transcripts.items():
        optimal_codons = transcript.get_optimal_codons()
        for amino_acid, (codon, usage_rate) in optimal_codons.items():
            output_data.append({
                'gene_name': gene_name,
                'amino_acid': amino_acid,
                'optimal_codon': codon,
                'usage_rate': f"{usage_rate:.4f}"
            })

    elapsed_time = time.time() - start_time
    return output_data, elapsed_time

def main():
    csv_files = [
        "csvs/P42_Brain_Ribo_rep1.csv",
        "csvs/P42_Brain_Ribo_rep2.csv",
        "csvs/P42_Heart_Ribo_rep1.csv",
        "csvs/P42_Heart_Ribo_rep2.csv",
        "csvs/P42_Kidney_Ribo_rep1.csv",
        "csvs/P42_Kidney_Ribo_rep2.csv",
        "csvs/P42_Liver_Ribo_rep1.csv",
        "csvs/P42_Lung_Ribo_rep1.csv",
        "csvs/P42_Lung_Ribo_rep2.csv",
        "csvs/P42_Retina_Ribo_rep2.csv",
    ]

    output_folder = "output_csvs"
    os.makedirs(output_folder, exist_ok=True)

    for filename in csv_files:
        process_file(filename, output_folder)

    print(f"Optimal codons for all files have been saved to the '{output_folder}' folder.")

if __name__ == "__main__":
    main()