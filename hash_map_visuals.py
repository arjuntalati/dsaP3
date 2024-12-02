import os
from hash_map import CodonHashMap, parse_csv, process_gene_data
import plotly.graph_objects as go


def process_folder(folder_path):
    codon_maps = {}
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".csv"):
            file_path = os.path.join(folder_path, file_name)
            print(f"Processing file: {file_name}")
            data = parse_csv(file_path)
            codon_map = CodonHashMap()
            process_gene_data(data, codon_map)
            codon_maps[file_name] = codon_map
    return codon_maps


def extract_data_for_plot(codon_map):
    labels = []
    parents = []
    values = []

    # transcripts
    for transcript_entry in codon_map.transcripts.map:
        if transcript_entry is None:
            continue
        transcript_id, transcript_data = transcript_entry[0], transcript_entry[1]
        labels.append(transcript_id)
        parents.append("")
        transcript_total = 0

        # amino acids
        for amino_acid_entry in transcript_data.map:
            if amino_acid_entry is None:
                continue
            amino_acid, amino_acid_data = amino_acid_entry[0], amino_acid_entry[1]
            labels.append(amino_acid)
            parents.append(transcript_id)
            amino_acid_total = 0

            # codons
            for codon_entry in amino_acid_data.map:
                if codon_entry is None:
                    continue
                codon, codon_count = codon_entry[0], codon_entry[1]
                labels.append(codon)
                parents.append(amino_acid)
                values.append(codon_count)
                amino_acid_total += codon_count

            values.append(amino_acid_total)
            transcript_total += amino_acid_total

        values.append(transcript_total)

    return labels, parents, values


def generate_plot(codon_maps):
    fig = go.Figure()

    for file_name, codon_map in codon_maps.items():
        labels, parents, values = extract_data_for_plot(codon_map)
        fig.add_trace(go.Sunburst(
            labels=labels,
            parents=parents,
            values=values,
            branchvalues="total",
            name=file_name
        ))

    fig.update_layout(
        title="Codon Usage Visualization",
        margin=dict(t=50, l=25, r=25, b=25),
    )
    fig.show()


if __name__ == "__main__":
    folder_path = "csvs"
    codon_maps = process_folder(folder_path)
    generate_plot(codon_maps)