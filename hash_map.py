import pandas as pd

class HashMap:
    def __init__(self, size=1000):
        self.size = size
        self.map = [None] * size

    def _hash(self, key):
        """Hash function to calculate the index."""
        return hash(key) % self.size

    def insert(self, key, value):
        """Insert or update a key-value pair."""
        index = self._hash(key)

        if self.map[index] is None:
            self.map[index] = [[key, value]]
        else:
            for pair in self.map[index]:
                if pair[0] == key:
                    pair[1] = value
                    return
            self.map[index].append([key, value])

    def get(self, key):
        """Retrieve a value by key."""
        index = self._hash(key)
        if self.map[index] is not None:
            for pair in self.map[index]:
                if pair[0] == key:
                    return pair[1]
        return None

    def keys(self):
        """Get all keys."""
        keys = []
        for bucket in self.map:
            if bucket is not None:
                for pair in bucket:
                    keys.append(pair[0])
        return keys

    def items(self):
        """Get all key-value pairs."""
        items = []
        for bucket in self.map:
            if bucket is not None:
                for pair in bucket:
                    items.append((pair[0], pair[1]))
        return items

    def __repr__(self):
        """String representation of the hash map."""
        result = "{"
        for bucket in self.map:
            if bucket is not None:
                for pair in bucket:
                    result += f"{repr(pair[0])}: {repr(pair[1])}, "
        result += "}"
        return result

class CodonHashMap:
    def __init__(self):
        self.transcripts = HashMap()

    def update_codon(self, transcript_id, amino_acid, codon):
        # Check if transcript exists
        transcript_data = self.transcripts.get(transcript_id)
        if not transcript_data:
            transcript_data = HashMap()
            self.transcripts.insert(transcript_id, transcript_data)

        # Check if amino acid exists in the transcript
        amino_acid_data = transcript_data.get(amino_acid)
        if not amino_acid_data:
            amino_acid_data = HashMap()
            transcript_data.insert(amino_acid, amino_acid_data)

        # Update codon count
        current_count = amino_acid_data.get(codon) or 0
        amino_acid_data.insert(codon, current_count + 1)

    def get_transcript(self, transcript_id):
        return self.transcripts.get(transcript_id)

    def __repr__(self):
        return repr(self.transcripts)

# Dictionary of codons to amino acids
CODON_TABLE = {
    "TTT": "F", "TTC": "F",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C",
    "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "TAA": "*", "TAG": "*", "TGA": "*",
}

# Parsing CSV file
def parse_csv(file_path):
    """Parse the CSV file and return a DataFrame."""
    data = pd.read_csv(file_path, delimiter=',', quotechar='"')
    return data

# Processing codon data
def process_gene_data(data, codon_map, gene_counts):
    total_rows = len(data)
    for idx, (_, row) in enumerate(data.iterrows(), 1):
        transcript_id = str(row['gene_name']).strip()
        sequence = str(row['sequence']).strip().upper()

        # Update gene count
        current_count = gene_counts.get(transcript_id) or 0
        gene_counts.insert(transcript_id, current_count + 1)

        # Split sequence into codons and update counts
        codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]
        for codon in codons:
            if codon in CODON_TABLE:
                amino_acid = CODON_TABLE[codon]
                codon_map.update_codon(transcript_id, amino_acid, codon)

    return gene_counts

# Normalize codon usage
def normalize_codon_usage(codon_map):
    normalized = {}

    for transcript_id, transcript_data in codon_map.transcripts.items():
        normalized[transcript_id] = {}

        for amino_acid, amino_acid_data in transcript_data.items():
            total_codons = sum(
                count for _, count in amino_acid_data.items()
            )

            if total_codons == 0:  # avoid division by zero
                continue

            normalized[transcript_id][amino_acid] = {
                codon: count / total_codons
                for codon, count in amino_acid_data.items()
            }

    return normalized

# Aggregate genome-wide optimality
def aggregate_optimality(codon_map, gene_counts):
    genome_wide = {}

    for transcript_id, transcript_data in codon_map.transcripts.items():
        scale_factor = gene_counts.get(transcript_id) or 0
        if scale_factor == 0:  # skip transcripts with no scale factor
            continue

        for amino_acid, amino_acid_data in transcript_data.items():
            if amino_acid not in genome_wide:
                genome_wide[amino_acid] = {}

            for codon, count in amino_acid_data.items():
                current_codon_count = genome_wide[amino_acid].get(codon, 0)
                genome_wide[amino_acid][codon] = current_codon_count + count * scale_factor

    # Find most optimal codons
    optimal_codon = {}
    for amino_acid, codon_counts in genome_wide.items():
        optimal_codon[amino_acid] = max(
            codon_counts.items(),
            key=lambda x: x[1] if x else 0
        )[0]

    return optimal_codon