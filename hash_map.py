class HashMap:
    def __init__(self, size=100):
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

    def __repr__(self):
        """String representation of the hash map."""
        return str(self.map)



class CodonHashMap:
    def __init__(self):
        self.transcripts = HashMap()
    
    def update_codon(self, transcript_id, amino_acid, codon):
        # check if transcript exists
        transcript_data = self.transcripts.get(transcript_id)
        if not transcript_data:
            self.transcripts.insert(transcript_id, HashMap())
            transcript_data = self.transcripts.get(transcript_id)
        
        # check if amino acid exists in the transcript
        amino_acid_data = transcript_data.get(amino_acid)
        if not amino_acid_data:
            transcript_data.insert(amino_acid, HashMap())
            amino_acid_data = transcript_data.get(amino_acid)
        
        # update codon count
        current_count = amino_acid_data.get(codon) or 0
        amino_acid_data.insert(codon, current_count + 1)
    
    def get_transcript(self, transcript_id):
        return self.transcripts.get(transcript_id)
    
    def __repr__(self):
        return repr(self.transcripts)


# mapping odons to amino acids
CODON_TABLE = {
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",  # Alanine
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",  # Arginine
    "AAT": "N", "AAC": "N",                          # Asparagine
    "GAT": "D", "GAC": "D",                          # Aspartic Acid
}


# parsing csv File
import csv

def parse_csv(file_path):
    data = []
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        print(f"CSV Headers: {reader.fieldnames}")  # Debugging: print the headers
        for row in reader:
            data.append(row)
    return data


# processing Codon Data
def process_gene_data(data, codon_map):
    gene_counts = HashMap()
    
    for row in data:
        transcript_id = row['gene_name']
        sequence = row['sequence']
        
        # update gene count
        current_count = gene_counts.get(transcript_id) or 0
        gene_counts.insert(transcript_id, current_count + 1)
        
        # split sequence into codons and update counts
        codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        for codon in codons:
            if codon in CODON_TABLE:
                amino_acid = CODON_TABLE[codon]
                codon_map.update_codon(transcript_id, amino_acid, codon)
    
    return gene_counts


# normalize Codon Usage
def normalize_codon_usage(codon_map):
    normalized = {}
    
    for transcript_id in codon_map.transcripts.map:
        if transcript_id is None:
            continue
        
        transcript_data = transcript_id[1]
        normalized[transcript_id[0]] = {}
        
        for amino_acid in transcript_data.map:
            if amino_acid is None:
                continue
            
            amino_acid_data = amino_acid[1]
            total_codons = sum(amino_acid_data.map[i][1] for i in range(len(amino_acid_data.map)) if amino_acid_data.map[i])
            normalized[transcript_id[0]][amino_acid[0]] = {
                codon[0]: codon[1] / total_codons for codon in amino_acid_data.map if codon
            }
    
    return normalized


# aggregate genome-wide optimality
def aggregate_optimality(codon_map, gene_counts):
    genome_wide = HashMap()
    
    for transcript_id in codon_map.transcripts.map:
        if transcript_id is None:
            continue
        
        transcript_data = transcript_id[1]
        scale_factor = gene_counts.get(transcript_id[0])
        
        for amino_acid in transcript_data.map:
            if amino_acid is None:
                continue
            
            amino_acid_data = amino_acid[1]
            for codon in amino_acid_data.map:
                if codon is None:
                    continue
                
                current_score = genome_wide.get(amino_acid[0]) or HashMap()
                current_codon_count = current_score.get(codon[0]) or 0
                current_score.insert(codon[0], current_codon_count + codon[1] * scale_factor)
                genome_wide.insert(amino_acid[0], current_score)
    
    # find most optimal codons
    optimal_codon = {}
    for amino_acid in genome_wide.map:
        if amino_acid is None:
            continue
        
        codon_counts = genome_wide.get(amino_acid[0])
        optimal_codon[amino_acid[0]] = max(codon_counts.map, key=lambda x: x[1] if x else 0)[0]
    
    return optimal_codon


# main Function
if __name__ == "__main__":
    file_path = "example.csv"  #example
    
    # parse data
    data = parse_csv(file_path)

    
    # initialize Codon Map
    codon_map = CodonHashMap()
    
    # process data and count genes
    gene_counts = process_gene_data(data, codon_map)
    
    # normalize codon usage
    normalized_usage = normalize_codon_usage(codon_map)
    
    # aggregate genome-wide optimality
    genome_optimality = aggregate_optimality(codon_map, gene_counts)
    
    # output results
    print("Normalized Codon Usage:", normalized_usage)
    print("Genome-Wide Optimal Codons:", genome_optimality)
