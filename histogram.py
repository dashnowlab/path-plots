"""
Annotated interactive histogram, built entirely in Python with Plotly.
Exports to a single self-contained HTML file -- no JS knowledge needed,
no build step, no server. Just open the .html file in a browser, or
commit it to a repo and open it from there.

Swap `data` below for your own (e.g. load from a CSV with pandas).
"""

import numpy as np
import plotly.graph_objects as go

# --- 1. Your data (replace this with pandas.read_csv(...) etc.) ---
np.random.seed(7)
data = np.concatenate([
    np.random.normal(loc=68, scale=9, size=1200),   # main cohort
    np.random.normal(loc=95, scale=5, size=150),    # a smaller secondary group
])

mean_val = float(np.mean(data))
p95_val = float(np.percentile(data, 95))

# --- 2. Build the histogram ---
fig = go.Figure()

fig.add_trace(go.Histogram(
    x=data,
    nbinsx=40,
    marker=dict(color="#4C78A8", line=dict(color="white", width=0.5)),
    hovertemplate="Range: %{x}<br>Count: %{y}<extra></extra>",
    name="Score distribution",
))

# --- 3. Reference lines ---
fig.add_vline(
    x=mean_val, line_width=2, line_dash="dash", line_color="#E45756",
)
fig.add_vline(
    x=p95_val, line_width=2, line_dash="dot", line_color="#54A24B",
)

# --- 4. Annotations (arrows + text callouts, done in plain Python) ---
fig.add_annotation(
    x=mean_val, y=1, yref="paper", yanchor="bottom",
    text=f"Mean = {mean_val:.1f}",
    showarrow=False,
    font=dict(color="#E45756", size=13),
)
fig.add_annotation(
    x=p95_val, y=0.92, yref="paper", yanchor="bottom",
    text=f"95th pct = {p95_val:.1f}",
    showarrow=False,
    font=dict(color="#54A24B", size=13),
)
fig.add_annotation(
    x=95, y=25,
    text="Secondary cluster<br>(worth investigating)",
    showarrow=True, arrowhead=2, ax=40, ay=-40,
    bgcolor="white", bordercolor="#666", borderwidth=1,
)

# --- 5. Layout / styling ---
fig.update_layout(
    title="Example annotated histogram (replace with your data)",
    xaxis_title="Value",
    yaxis_title="Count",
    bargap=0.02,
    template="plotly_white",
    width=900,
    height=520,
)

# --- 6. Export to a single, portable HTML file ---
fig.write_html(
    "annotated_histogram.html",
    include_plotlyjs="cdn",   # keeps file small; loads Plotly from CDN
    full_html=True,
)

print("Wrote annotated_histogram.html")