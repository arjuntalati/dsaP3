import os
import pandas as pd
import circlify
from flask import Flask, render_template_string, jsonify, request, redirect, url_for
import threading
import json
import gzip
import base64
import time
from collections import defaultdict

# Import data processing functions from max_heap.py
import max_heap

app = Flask(__name__)

# Global variables to store processing status, times, and data
processing_status = {}
processing_times = {}
processed_data = {}
data_lock = threading.Lock()

# Function to process multiple files in a separate thread
def process_files_thread(csv_files):
    for filename in csv_files:
        sample_name = os.path.splitext(os.path.basename(filename))[0]
        with data_lock:
            processing_status[sample_name] = 'Processing'
        # Process one file
        output_data, elapsed_time = max_heap.process_file(filename)
        with data_lock:
            processed_data[sample_name] = output_data
            processing_times[sample_name] = elapsed_time
            processing_status[sample_name] = 'Completed'

# Load data from memory (if needed)
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
            # "csvs/P42_Brain_Ribo_rep1.csv",
            # "csvs/P42_Brain_Ribo_rep2.csv",
            # "csvs/P42_Heart_Ribo_rep1.csv",
            # "csvs/P42_Heart_Ribo_rep2.csv",
            # "csvs/P42_Kidney_Ribo_rep1.csv",
            # "csvs/P42_Kidney_Ribo_rep2.csv",
            # "csvs/P42_Liver_Ribo_rep1.csv",
            "csvs/P42_Lung_Ribo_rep1.csv",
            "csvs/P42_Lung_Ribo_rep2.csv",
            "csvs/P42_Retina_Ribo_rep2.csv",
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
        # Redirect to comparison page with selected samples
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

@app.route("/compare")
def compare():
    sample1 = request.args.get('sample1')
    sample2 = request.args.get('sample2')

    if not sample1 or not sample2:
        return "Two samples are required for comparison.", 400

    with data_lock:
        data_sample1 = processed_data.get(sample1)
        data_sample2 = processed_data.get(sample2)

    if not data_sample1 or not data_sample2:
        return "Sample data not found.", 404

    # Compress data and encode it to base64 for safe embedding
    data_bytes_sample1 = json.dumps(data_sample1).encode('utf-8')
    compressed_data_sample1 = base64.b64encode(gzip.compress(data_bytes_sample1)).decode('utf-8')

    data_bytes_sample2 = json.dumps(data_sample2).encode('utf-8')
    compressed_data_sample2 = base64.b64encode(gzip.compress(data_bytes_sample2)).decode('utf-8')

    return render_template_string("""
        <!-- Back to Sample Selection button -->
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

        <!-- Include necessary libraries -->
        <script src="https://cdnjs.cloudflare.com/ajax/libs/pako/2.0.4/pako.min.js"></script>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <script src="https://d3js.org/d3.v6.min.js"></script>

        <!-- Decompress and store data -->
        <script>
            const compressedDataSample1 = {{ compressed_data_sample1 | tojson }};
            const compressedDataSample2 = {{ compressed_data_sample2 | tojson }};

            function decompressData(compressedData) {
                const compressedBytes = Uint8Array.from(atob(compressedData), c => c.charCodeAt(0));
                const decompressedBytes = pako.inflate(compressedBytes);
                const decoder = new TextDecoder('utf-8');
                const jsonString = decoder.decode(decompressedBytes);
                return JSON.parse(jsonString);
            }

            const dataSample1 = decompressData(compressedDataSample1);
            const dataSample2 = decompressData(compressedDataSample2);
            const samplesData = {
                "{{ sample1 }}": dataSample1,
                "{{ sample2 }}": dataSample2
            };

            // Debugging: Log the data to ensure it's correctly decompressed
            console.log('Data for {{ sample1 }}:', dataSample1);
            console.log('Data for {{ sample2 }}:', dataSample2);

            const samples = [ "{{ sample1 }}", "{{ sample2 }}" ];
            const history = {};  // To keep track of navigation history for each sample

            samples.forEach(function(sampleName) {
                const containerId = 'circle-container-' + sampleName;
                history[sampleName] = [];  // Initialize history for each sample
                loadVisualization(sampleName, 'gene_name', null, containerId);
            });

            function loadVisualization(sampleName, level, parentName, containerId) {
                // Filter data on the client side
                const data = samplesData[sampleName];
                const filteredData = filterData(data, level, parentName);
                const plotData = generateCirclePacking(filteredData, level);

                if (!plotData) {
                    console.error("No data available for the selected level and parent.");
                    return;
                }

                plotVisualization(plotData, level, parentName, sampleName, containerId);
            }

            function filterData(data, level, parentName) {
                const nextLevelMap = {
                    "gene_name": "amino_acid",
                    "amino_acid": "optimal_codon"
                };
                const groupColMap = {
                    "gene_name": "gene_name",
                    "amino_acid": "amino_acid",
                    "optimal_codon": "optimal_codon"
                };
                const groupCol = groupColMap[level];
                const prevLevelMap = {};
                for (let key in nextLevelMap) {
                    prevLevelMap[nextLevelMap[key]] = key;
                }

                let filteredData = data;

                if (level !== "gene_name" && parentName != null) {
                    const previousLevel = prevLevelMap[level];
                    const prevGroupCol = groupColMap[previousLevel];
                    filteredData = data.filter(item => item[prevGroupCol] === parentName);
                }

                return filteredData;
            }

            function generateCirclePacking(data, level) {
                if (!data || data.length === 0) return null;

                const groupColMap = {
                    "gene_name": "gene_name",
                    "amino_acid": "amino_acid",
                    "optimal_codon": "optimal_codon"
                };
                const groupCol = groupColMap[level];

                // Group data
                const groupedData = d3.rollups(
                    data,
                    v => ({
                        count: v.length,
                        usage_rate: d3.sum(v, d => d.usage_rate || 0)
                    }),
                    d => d[groupCol]
                );

                // Convert to hierarchical data
                const rootData = {
                    name: "root",
                    children: groupedData.map(([key, value]) => ({
                        name: key,
                        value: level === 'optimal_codon' ? value.usage_rate : value.count
                    }))
                };

                // Create hierarchy
                const root = d3.hierarchy(rootData)
                    .sum(d => d.value)
                    .sort((a, b) => b.value - a.value);

                // Create circle packing layout
                const pack = d3.pack()
                    .size([600, 600]) // Adjust size as needed
                    .padding(3);

                const nodes = pack(root).leaves();

                // Prepare plot data
                const plotData = nodes.map(node => ({
                    x: node.x,
                    y: node.y,
                    r: node.r,
                    id: node.data.name,
                    datum: node.value,
                    level: level
                }));

                return plotData;
            }

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

                const annotations = plotData.map(p => ({
                    x: p.x,
                    y: p.y,
                    text: p.id,
                    showarrow: false,
                    font: {
                        size: Math.max(Math.min(p.r / 2, 18), 6),
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
                        const containerParent = document.getElementById(containerId).parentElement;
                        containerParent.insertBefore(backButton, containerParent.firstChild);
                    }
                } else {
                    // Remove Back button if at top level
                    const backButton = document.getElementById('back-button-' + sampleName);
                    if (backButton) {
                        backButton.remove();
                    }
                }
            }
        </script>
    """, sample1=sample1, sample2=sample2,
         compressed_data_sample1=compressed_data_sample1,
         compressed_data_sample2=compressed_data_sample2)

if __name__ == "__main__":
    app.run(port=5002, debug=True)