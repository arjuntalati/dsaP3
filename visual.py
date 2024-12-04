import os
import pandas as pd
import circlify
import plotly.graph_objects as go
from flask import Flask, render_template_string, jsonify, request, redirect, url_for

app = Flask(__name__)

# Load all data from CSVs
def load_data(csv_folder):
    data = {}
    for file_name in os.listdir(csv_folder):
        if file_name.endswith(".csv"):
            file_path = os.path.join(csv_folder, file_name)
            sample_name = os.path.splitext(file_name)[0]
            df = pd.read_csv(file_path)
            df['sample_name'] = sample_name  # Add sample_name to DataFrame
            data[sample_name] = df
    return data

# Generate circle packing data
def generate_circle_packing(data, level, parent_name=None):
    next_level_map = {
        "gene_name": "amino_acid",
        "amino_acid": "optimal_codon"
    }
    group_col_map = {
        "gene_name": "gene_name",
        "amino_acid": "amino_acid",
        "optimal_codon": "optimal_codon"
    }
    group_col = group_col_map.get(level)
    prev_level_map = {v: k for k, v in next_level_map.items()}

    # Adjusted filtering logic
    if level == "gene_name":
        filtered_data = data
    elif parent_name is not None:
        previous_level = prev_level_map.get(level)
        prev_group_col = group_col_map.get(previous_level)
        filtered_data = data[data[prev_group_col] == parent_name]
    else:
        filtered_data = data

    if filtered_data.empty:
        return None, None

    # Use counts or usage_rate based on level
    if level == "optimal_codon":
        grouped = filtered_data.groupby(group_col).agg({"usage_rate": "sum"}).reset_index()
        size_col = 'usage_rate'
    else:
        grouped = filtered_data.groupby(group_col).size().reset_index(name='count')
        size_col = 'count'

    # Prepare data for circlify
    circle_data = [
        {"id": row[group_col], "datum": row[size_col]}
        for _, row in grouped.iterrows()
    ]
    if not circle_data:
        return None, None

    # Generate circles
    circles = circlify.circlify(
        circle_data,
        show_enclosure=False,
        target_enclosure=circlify.Circle(x=0, y=0, r=1)
    )

    # Prepare data for Plotly
    plot_data = []
    for circle in circles:
        if circle.level != 1:  # We only need circles at level 1
            continue
        id_name = circle.ex["id"]
        plot_data.append({
            "x": circle.x,
            "y": circle.y,
            "r": circle.r,
            "id": id_name,
            "datum": circle.ex["datum"],
            "level": level
        })

    return plot_data, grouped[group_col].tolist()

@app.route("/", methods=["GET", "POST"])
def index():
    csv_folder = "output_csvs"
    all_data = load_data(csv_folder)
    samples = list(all_data.keys())

    if request.method == "POST":
        # Get selected samples from form
        selected_samples = request.form.getlist('samples')
        if len(selected_samples) != 2:
            return "Please select exactly two samples.", 400
        # Redirect to compare page with selected samples
        return redirect(url_for('compare', sample1=selected_samples[0], sample2=selected_samples[1]))

    return render_template_string("""
        <h1>Select Two Samples to Compare</h1>
        <form method="post">
            {% for sample in samples %}
                <input type="checkbox" name="samples" value="{{ sample }}"> {{ sample }}<br>
            {% endfor %}
            <br>
            <input type="submit" value="Compare">
        </form>
    """, samples=samples)

@app.route("/compare")
def compare():
    sample1 = request.args.get('sample1')
    sample2 = request.args.get('sample2')

    if not sample1 or not sample2:
        return "Two samples are required for comparison.", 400

    return render_template_string("""
        <h1>Comparing Samples: {{ sample1 }} and {{ sample2 }}</h1>
        <div id="samples-container" style="display: flex; flex-wrap: wrap;">
            <div style="margin: 20px;">
                <h2>{{ sample1 }}</h2>
                <div id="circle-container-{{ sample1 }}"></div>
            </div>
            <div style="margin: 20px;">
                <h2>{{ sample2 }}</h2>
                <div id="circle-container-{{ sample2 }}"></div>
            </div>
        </div>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <script>
            const samples = [ "{{ sample1 }}", "{{ sample2 }}" ];
            const history = {};  // To keep track of navigation history for each sample

            samples.forEach(function(sampleName) {
                const containerId = 'circle-container-' + sampleName;
                history[sampleName] = [];  // Initialize history for each sample
                loadVisualization(sampleName, 'gene_name', null, containerId);
            });

            function plotVisualization(plotData, level, parentName, sampleName, containerId) {
                const shapes = plotData.map(p => ({
                    type: 'circle',
                    xref: 'x',
                    yref: 'y',
                    x0: p.x - p.r,
                    y0: p.y - p.r,
                    x1: p.x + p.r,
                    y1: p.y + p.r,
                    line: {
                        color: 'black',
                    },
                    fillcolor: 'lightblue',
                    opacity: 0.6,
                }));

                const annotations = plotData.map(p => ({
                    x: p.x,
                    y: p.y,
                    text: p.id,
                    showarrow: false,
                    font: {
                        size: 12,
                        color: 'black'
                    },
                }));

                const clickTrace = {
                    x: plotData.map(p => p.x),
                    y: plotData.map(p => p.y),
                    mode: 'markers',
                    marker: {
                        size: 0.1,
                        color: 'rgba(0,0,0,0)',
                    },
                    text: plotData.map(p => `${p.id}<br>Value: ${p.datum}`),
                    customdata: plotData.map(p => [p.id, level]),
                    hoverinfo: 'text',
                };

                const layout = {
                    title: `${parentName || sampleName} - ${level}`,
                    showlegend: false,
                    xaxis: { showgrid: false, zeroline: false, visible: false },
                    yaxis: { showgrid: false, zeroline: false, visible: false },
                    width: 600,
                    height: 600,
                    shapes: shapes,
                    annotations: annotations,
                    hovermode: 'closest'
                };

                const data = [clickTrace];

                Plotly.newPlot(containerId, data, layout);

                // Attach click event
                const container = document.getElementById(containerId);
                container.on('plotly_click', function(event) {
                    const clickedPoint = event.points[0];
                    const clickedName = clickedPoint.customdata[0];
                    const currentLevel = clickedPoint.customdata[1];

                    const nextLevelMap = {"gene_name": "amino_acid", "amino_acid": "optimal_codon"};
                    const nextLevel = nextLevelMap[currentLevel];

                    if (nextLevel) {
                        let newParentName = clickedName;
                        // Save current state to history
                        history[sampleName].push({ level: currentLevel, parentName: parentName });
                        loadVisualization(sampleName, nextLevel, newParentName, containerId);
                    }
                });

                // Add Back button if not at top level
                if (history[sampleName].length > 0) {
                    const backButtonId = 'back-button-' + sampleName;
                    if (!document.getElementById(backButtonId)) {
                        const backButton = document.createElement('button');
                        backButton.id = backButtonId;
                        backButton.innerText = 'Back';
                        backButton.onclick = function() {
                            const previousState = history[sampleName].pop();
                            if (previousState) {
                                loadVisualization(sampleName, previousState.level, previousState.parentName, containerId);
                            }
                        };
                        container.parentElement.insertBefore(backButton, container);
                    }
                } else {
                    // Remove Back button if at top level
                    const backButton = document.getElementById('back-button-' + sampleName);
                    if (backButton) {
                        backButton.remove();
                    }
                }
            }

            function loadVisualization(sampleName, level, parentName, containerId) {
                let url = `/api/visualize?level=${level}&sample=${encodeURIComponent(sampleName)}`;
                if (parentName) {
                    url += `&parent=${encodeURIComponent(parentName)}`;
                }

                fetch(url).then(response => response.json()).then(data => {
                    if (data.error) {
                        console.error(data.error);
                        return;
                    }
                    plotVisualization(data.plot_data, data.level, parentName, sampleName, containerId);
                });
            }
        </script>
    """, sample1=sample1, sample2=sample2)

@app.route("/api/visualize")
def api_visualize():
    csv_folder = "output_csvs"
    all_data = load_data(csv_folder)

    level = request.args.get("level", "gene_name")
    parent_name = request.args.get("parent", None)
    sample_name = request.args.get("sample", None)

    if not sample_name:
        return jsonify({"error": "Sample name is required"}), 400

    if sample_name not in all_data:
        return jsonify({"error": f"Sample '{sample_name}' not found"}), 404

    data = all_data[sample_name]

    plot_data, _ = generate_circle_packing(data, level, parent_name)
    return jsonify({"plot_data": plot_data, "level": level})

if __name__ == "__main__":
    app.run(debug=True)