import os
import pandas as pd
import circlify
import plotly.graph_objects as go
from flask import Flask, render_template_string, request

app = Flask(__name__)

# Load all data from CSVs
def load_data(csv_folder):
    data = {}
    for file_name in os.listdir(csv_folder):
        if file_name.endswith(".csv"):
            file_path = os.path.join(csv_folder, file_name)
            sample_name = os.path.splitext(file_name)[0]
            data[sample_name] = pd.read_csv(file_path)
    return data

# Generate circle packing visualization
def generate_circle_packing(data, level, parent_name=None):
    """
    Generate circle packing visualization for the current level.
    - data: DataFrame containing the hierarchical data.
    - level: Current hierarchy level ("file", "gene_name", "amino_acid").
    - parent_name: Parent entity name (gene name, amino acid name, etc.).
    """
    next_level_map = {"file": "gene_name", "gene_name": "amino_acid", "amino_acid": "optimal_codon"}
    next_level = next_level_map.get(level)

    # Filter data based on the level
    if level == "file":
        filtered_data = data
        group_col = "gene_name"
    elif level == "gene_name":
        filtered_data = data[data["gene_name"] == parent_name]
        group_col = "amino_acid"
    elif level == "amino_acid":
        filtered_data = data[data["amino_acid"] == parent_name]
        group_col = "optimal_codon"
    else:
        return None, None

    grouped = filtered_data.groupby(group_col).agg({"usage_rate": "sum"}).reset_index()

    # Prepare data for circlify
    circle_data = [{"id": row[group_col], "datum": row["usage_rate"]} for _, row in grouped.iterrows()]
    if not circle_data:
        return None, None

    # Generate circles
    circles = circlify.circlify(circle_data, show_enclosure=True)

    # Create a Plotly figure
    fig = go.Figure()
    children = []

    for circle in circles:
        if circle.level == 0:  # Skip enclosing circle
            continue

        id_name = circle.ex["id"]
        children.append(id_name)

        fig.add_trace(go.Scatter(
            x=[circle.x],
            y=[circle.y],
            mode="markers+text",
            marker=dict(
                size=circle.r * 500,
                color="lightblue",
                opacity=0.6,
                line=dict(color="black", width=1)
            ),
            text=id_name,
            textposition="middle center",
            hoverinfo="text",
            hovertext=f"{id_name}",
        ))

    fig.update_layout(
        title=f"{parent_name or 'Samples'} - {level}",
        showlegend=False,
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False),
        width=800,
        height=800
    )

    return fig.to_html(full_html=False), children

# Flask routes
@app.route("/")
def index():
    csv_folder = "output_csvs"
    samples = os.listdir(csv_folder)
    return render_template_string("""
        <h1>Sample Files</h1>
        <div>
            {% for sample in samples %}
            <div style="display: inline-block; margin: 10px;">
                <a href="/circle/file/{{ sample }}" style="text-decoration: none;">
                    <div style="width: 100px; height: 100px; background-color: lightblue; border-radius: 50%; display: flex; align-items: center; justify-content: center;">
                        {{ sample }}
                    </div>
                </a>
            </div>
            {% endfor %}
        </div>
    """, samples=samples)

@app.route("/circle/<level>/<name>")
def circle(level, name):
    csv_folder = "output_csvs"
    all_data = load_data(csv_folder)

    # Identify the correct file and load data
    if level == "file":
        data = all_data[os.path.splitext(name)[0]]
        parent_name = name
    else:
        # For gene_name or amino_acid, infer the parent file
        sample_name = request.args.get("sample")
        data = all_data[sample_name]
        parent_name = name

    # Generate visualization
    html, children = generate_circle_packing(data, level, parent_name)

    # Map level for next navigation
    next_level_map = {"file": "gene_name", "gene_name": "amino_acid", "amino_acid": "optimal_codon"}
    next_level = next_level_map.get(level)

    return render_template_string("""
        <h1>Circle Packing for {{ name }}</h1>
        <div>{{ html|safe }}</div>
        <br>
        {% if children %}
        <h2>Click on a circle to explore deeper:</h2>
        <div>
            {% for child in children %}
            <a href="/circle/{{ next_level }}/{{ child }}?sample={{ parent_name }}" style="margin: 10px; display: inline-block; text-decoration: none;">
                <div style="width: 100px; height: 100px; background-color: lightgreen; border-radius: 50%; display: flex; align-items: center; justify-content: center;">
                    {{ child }}
                </div>
            </a>
            {% endfor %}
        </div>
        {% endif %}
        <br>
        <a href="/">Back to Sample Selection</a>
    """, name=name, html=html, children=children, next_level=next_level, parent_name=os.path.splitext(name)[0])

if __name__ == "__main__":
    app.run(debug=True)