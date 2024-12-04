from flask import Flask, render_template_string

app = Flask(__name__)

# NOTE: to properly run this script, the hash_map_visuals.py script and the max_heap_visuals.py script must also be running on your machine
@app.route("/")
def index():
    return render_template_string("""
        <h1>Choose which data structure to look at:</h1>
        <div>
            <button onclick="window.location.href='/hash_map'">Hash Map Visualizations</button>
            <button onclick="window.location.href='/max_heap'">Max Heap Visualizations</button>
        </div>
    """)

@app.route("/hash_map")
def hash_map_redirect():
    return render_template_string("""
        <h1>Hash Map Visualizations</h1>
        <p>Go to Hash Map Visualizations</p>
        <script>
            window.location.href = 'http://127.0.0.1:5001';
        </script>
    """)

@app.route("/max_heap")
def max_heap_redirect():
    return render_template_string("""
        <h1>Max Heap</h1>
        <p>Go to Max Heap Visualizations</p>
        <script>
            window.location.href = 'http://127.0.0.1:5002';
        </script>
    """)

if __name__ == "__main__":
    app.run(port=5000, debug=True)