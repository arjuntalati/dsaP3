import csv
import pandas as pd
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
        # value is the priority (usage rate), key is the codon
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
        self.amino_acid_codons = {}  # amino_acid -> {codon: count}
        self.total_amino_acid_counts = {}  # amino_acid -> total count
        self.heaps = {}  # amino_acid -> MaxHeap

    def add_sequence(self, sequence):
        # process the sequence into codons
        codons = [sequence[i:i+3] for i in range(0, len(sequence) - 2, 3)]
        for codon in codons:
            amino_acid = codon_to_amino_acid(codon)
            if amino_acid is None:
                continue  # skip invalid codons
            # initialize  the dictoinaries
            if amino_acid not in self.amino_acid_codons:
                self.amino_acid_codons[amino_acid] = {}
                self.total_amino_acid_counts[amino_acid] = 0
            # count codon usage
            self.amino_acid_codons[amino_acid][codon] = self.amino_acid_codons[amino_acid].get(codon, 0) + 1
            self.total_amino_acid_counts[amino_acid] += 1

    def calculate_usage_rates(self):
        # for each amino acid, calculate usage rates and build the max heap
        for amino_acid in self.amino_acid_codons:
            total_count = self.total_amino_acid_counts[amino_acid]
            heap = MaxHeap()
            for codon in self.amino_acid_codons[amino_acid]:
                count = self.amino_acid_codons[amino_acid][codon]
                usage_rate = count / total_count
                # insert into max heap
                heap.insert(codon, usage_rate)
            self.heaps[amino_acid] = heap

    def get_optimal_codons(self):
        # return the codon with the highest usage rate for each amino acid
        optimal_codons = {}
        for amino_acid in self.heaps:
            heap = self.heaps[amino_acid]
            if not heap.is_empty():
                usage_rate, codon = heap.extract_max()
                optimal_codons[amino_acid] = (codon, usage_rate)
        return optimal_codons
    
    def scale_usage_rates(self, scale_factor):
        for amino_acid in self.heaps:
            heap = self.heaps[amino_acid]
            adjusted_heap = MaxHeap()
            temp_list = []
            total_scaled_usage = 0
            # extract all items, adjust usage rates, calculate total
            while not heap.is_empty():
                usage_rate, codon = heap.extract_max()
                adjusted_rate = usage_rate * scale_factor
                total_scaled_usage += adjusted_rate
                temp_list.append((adjusted_rate, codon))
            # normalize the adjusted rates
            normalized_temp_list = []
            for adjusted_rate, codon in temp_list:
                normalized_rate = adjusted_rate / total_scaled_usage if total_scaled_usage != 0 else 0
                normalized_temp_list.append((normalized_rate, codon))
            # reinsert normalized items
            for normalized_rate, codon in normalized_temp_list:
                adjusted_heap.insert(codon, normalized_rate)
            self.heaps[amino_acid] = adjusted_heap



def main():
    transcripts = {}  # gene_name -> Transcript instance
    transcript_counts = {}  # gene_name -> occurrence count

    # list of CSV files
    csv_files = [
        "P42_Liver_Ribo_rep2.1_no_headers_gene_only_processed.csv",
        "P42_Retina_Ribo_rep1.1_no_headers_gene_only_processed.csv",
        # and the rest of the files...
    ]

    for filename in csv_files:
        with open(filename, 'r', encoding='utf-8') as file:
            reader = csv.DictReader(file, delimiter=',')
            for row in reader:
                gene_name = row['gene_name']
                sequence = row['sequence']
                # update transcript counts
                transcript_counts[gene_name] = transcript_counts.get(gene_name, 0) + 1
                # Initialize Transcript instance 
                if gene_name not in transcripts:
                    transcripts[gene_name] = Transcript(gene_name)
                # add sequence to the transcript
                transcripts[gene_name].add_sequence(sequence)

    # for each transcript, calculate usage rates
    for gene_name in transcripts:
        transcripts[gene_name].calculate_usage_rates()

    # scale usage rates by transcript frequency
    for gene_name in transcripts:
        frequency = transcript_counts[gene_name]
        transcripts[gene_name].scale_usage_rates(frequency)

    # store optimal codons for all proteins
    all_optimal_codons = {}

    # get optimal codons for all proteins
    for gene_name in transcripts:
        transcript = transcripts[gene_name]
        optimal_codons = transcript.get_optimal_codons()
        all_optimal_codons[gene_name] = optimal_codons

    # save all optimal codons
    output_filename = "C:\\Users\\adiaz\\Downloads\\optimal_codons.csv"
    with open(output_filename, 'w', newline='') as csvfile:
        fieldnames = ['gene_name', 'amino_acid', 'optimal_codon', 'usage_rate']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for gene_name, codons in all_optimal_codons.items():
            for amino_acid, (codon, usage_rate) in codons.items():
                writer.writerow({
                    'gene_name': gene_name,
                    'amino_acid': amino_acid,
                    'optimal_codon': codon,
                    'usage_rate': f"{usage_rate:.4f}"
                })

    print(f"Optimal codons for all proteins have been saved to '{output_filename}'.")

if __name__ == "__main__":
    main()