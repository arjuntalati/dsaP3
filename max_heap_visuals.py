import os
import pandas as pd
import circlify
from flask import Flask, render_template_string, jsonify, request, redirect, url_for
import threading
import time
from collections import defaultdict

# data processing functions from max_heap.py
import max_heap

app = Flask(__name__)

# global variables to store processing status, times, and data
processing_status = {}
processing_times = {}
processed_data = {}
data_lock = threading.Lock()

def process_files_thread(csv_files):
    for filename in csv_files:
        sample_name = os.path.splitext(os.path.basename(filename))[0]
        with data_lock:
            processing_status[sample_name] = 'Processing'
        # process one file
        output_data, elapsed_time = max_heap.process_file(filename)
        with data_lock:
            processed_data[sample_name] = output_data
            processing_times[sample_name] = elapsed_time
            processing_status[sample_name] = 'Completed'

def load_data_from_memory():
    data = {}
    with data_lock:
        for sample_name, output_data in processed_data.items():
            df = pd.DataFrame(output_data)
            df['sample_name'] = sample_name
            data[sample_name] = df
    return data

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        csv_files = [
            "csvs/P42_Brain_Ribo_rep1.1_no_headers_gene_only_processed.csv",
            "csvs/P42_Brain_Ribo_rep2.1_no_headers_gene_only_processed.csv",
            "csvs/P42_Heart_Ribo_rep1.1_no_headers_gene_only_processed.csv",
            "csvs/P42_Heart_Ribo_rep2.1_no_headers_gene_only_processed.csv",
            "csvs/P42_Kidney_Ribo_rep1.1_no_headers_gene_only_processed.csv",
            "csvs/P42_Kidney_Ribo_rep2.1_no_headers_gene_only_processed.csv",
            "csvs/P42_Liver_Ribo_rep1.1_no_headers_gene_only_processed.csv",
            "csvs/P42_Liver_Ribo_rep2.1_no_headers_gene_only_processed.csv",
            "csvs/P42_Lung_Ribo_rep1.1_no_headers_gene_only_processed.csv",
            "csvs/P42_Lung_Ribo_rep2.1_no_headers_gene_only_processed.csv",
            "csvs/P42_Retina_Ribo_rep1.1_no_headers_gene_only_processed.csv",
            "csvs/P42_Retina_Ribo_rep2.1_no_headers_gene_only_processed.csv",
        ]
        threading.Thread(target=process_files_thread, args=(csv_files,), daemon=True).start()
        return redirect(url_for('processing_status_page'))

    return render_template_string("""
        <h1>Start Data Processing</h1>
        <form method="post">
            <input type="submit" value="Process Data">
        </form>
    """)

@app.route("/processing_status")
def processing_status_page():
    with data_lock:
        status = dict(processing_status)
        times = dict(processing_times)
    all_completed = all(s == 'Completed' for s in status.values()) and status != {}
    return render_template_string("""
        <h1>Data Processing Status</h1>
        <ul>
            {% for sample_name, sample_status in status.items() %}
                <li>{{ sample_name }}: {{ sample_status }}
                    {% if sample_status == 'Completed' %}
                        - Time taken: {{ times[sample_name]|round(2) }} seconds
                    {% endif %}
                </li>
            {% endfor %}
        </ul>
        {% if not all_completed %}
            <script>
                setTimeout(function(){
                    window.location.reload(1);
                }, 5000);
            </script>
        {% else %}
            <a href="{{ url_for('select_samples') }}">Proceed to Sample Selection</a>
        {% endif %}
    """, status=status, times=times, all_completed=all_completed)

@app.route("/select_samples", methods=["GET", "POST"])
def select_samples():
    with data_lock:
        samples = list(processed_data.keys())
    if request.method == "POST":
        selected_samples = request.form.getlist('samples')
        if len(selected_samples) != 2:
            return "Please select exactly two samples.", 400
        # go to comparison page with selected samples
        return redirect(url_for('compare', sample1=selected_samples[0], sample2=selected_samples[1]))
    return render_template_string("""
        <h1>Select Two Samples to Compare</h1>
        <form method="post">
            {% for sample in samples %}
                <input type="checkbox" name="samples" value="{{ sample }}"> {{ sample }}<br>
            {% endfor %}
            <br>
            <input type="submit" value="Compare">
            <button onclick="window.location.href='/'" type="button">Home</button>
        </form>
    """, samples=samples)

# circle packing
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

    # use counts or usage_rate based on level for circle sizes
    if level == "optimal_codon":
        grouped = filtered_data.groupby(group_col).agg({"usage_rate": "sum"}).reset_index()
        size_col = 'usage_rate'
    else:
        grouped = filtered_data.groupby(group_col).size().reset_index(name='count')
        size_col = 'count'

    # prepare for circlify
    circle_data = [
        {"id": row[group_col], "datum": row[size_col]}
        for _, row in grouped.iterrows()
    ]
    if not circle_data:
        return None, None

    circles = circlify.circlify(
        circle_data,
        show_enclosure=False,
        target_enclosure=circlify.Circle(x=0, y=0, r=1)
    )

    # prepare for Plotly
    plot_data = []
    for circle in circles:
        if circle.level != 1:
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

@app.route("/compare")
def compare():
    sample1 = request.args.get('sample1')
    sample2 = request.args.get('sample2')

    if not sample1 or not sample2:
        return "Two samples are required for comparison.", 400

    return render_template_string("""
        <!-- Changed the button to go back to sample selection -->
        <button onclick="window.location.href='{{ url_for('select_samples') }}'">Back to Sample Selection</button>
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
                const levelTitleMap = {
                    "gene_name": "Genes",
                    "amino_acid": "Amino Acids",
                    "optimal_codon": "Codons"
                };

                const titleText = `${levelTitleMap[level] || level}`;

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

                // Function to calculate font size based on circle radius and text length
                function calculateFontSize(r, textLength) {
                    const scalingFactor = 30; // Adjust this value as needed
                    const fontSize = (r * scalingFactor) / textLength;
                    const minFontSize = 6;
                    const maxFontSize = 18;
                    return Math.max(minFontSize, Math.min(fontSize, maxFontSize));
                }

                const annotations = plotData.map(p => {
                    let fontSize;
                    if (level === 'gene_name') {
                        fontSize = calculateFontSize(p.r, p.id.length);
                    } else {
                        fontSize = 12; // Fixed font size for other layers
                    }
                    return {
                        x: p.x,
                        y: p.y,
                        text: p.id,
                        showarrow: false,
                        font: {
                            size: fontSize,
                            color: 'black'
                        },
                    };
                });

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
                    title: titleText,
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
    level = request.args.get("level", "gene_name")
    parent_name = request.args.get("parent", None)
    sample_name = request.args.get("sample", None)

    if not sample_name:
        return jsonify({"error": "Sample name is required"}), 400

    with data_lock:
        if sample_name not in processed_data:
            return jsonify({"error": f"Sample '{sample_name}' not found"}), 404
        data_list = processed_data[sample_name]
        data = pd.DataFrame(data_list)

    plot_data, _ = generate_circle_packing(data, level, parent_name)
    return jsonify({"plot_data": plot_data, "level": level})

if __name__ == "__main__":
    app.run(debug=True)