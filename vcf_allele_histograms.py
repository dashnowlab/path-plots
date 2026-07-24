"""
Per-locus histograms from a TRGT VCF, built with Plotly.
Reads the VCF directly (no intermediate file needed) with Python's stdlib,
pulls per-sample FORMAT fields for every locus, and writes a single
self-contained HTML file with a histogram switchable via two dropdowns:

  Locus   which tandem repeat (TRID) to show
  Metric  AL  Length of each allele
          MC  Motif count per allele, one histogram per motif
          MS  Motif span (bp) per allele, one histogram per motif
          AP  Allele purity per allele

MC and MS are broken out per motif (colored, with a legend) rather than
summed across a locus's motifs: for loci with more than one motif (e.g.
CANVAS_RFC1's pathogenic AAGGG vs benign ACAGG), summing would hide exactly
the distinction that matters.

All post-load rendering (which traces are visible, their data, the x-axis
range/title, the plot title) is driven by a single JS `redraw()` function
injected at the bottom, so there's one source of truth for both the initial
paint and every dropdown click -- Python only builds placeholder traces and
the dropdown/annotation skeleton.

Benign/intermediate/pathogenic ranges are always pulled from whichever
STRchive release is currently latest (never a pinned version) and shown as
shaded bands on the AL and MC metrics; the release used is printed to the
console and shown in the plot's footer. STRchive's ranges are repeat-unit
copy numbers, so the AL bands are an approximation (scaled by motif length).

Swap `vcf_path` below for your own TRGT VCF (plain-text or .vcf.gz).
"""

import gzip
import json
import re
import urllib.request

import plotly.graph_objects as go

vcf_path = "data/100HPRC.trgt-v0.8.0.STRchive.sorted.vcf"

METRICS = ["MC", "AL", "MS", "AP"]  # MC first: it's the default/initial view
INTEGER_METRICS = {"AL", "MC", "MS"}  # AP is the only continuous (0-1) metric
PER_MOTIF_METRICS = {"MC", "MS"}

STRCHIVE_LATEST_RELEASE_URL = "https://api.github.com/repos/dashnowlab/STRchive/releases/latest"
TIER_COLORS = {"benign": "#54A24B", "intermediate": "#EECA3B", "pathogenic": "#E45756"}

# Each motif is colored by its pathogenic/benign/unknown classification, but a
# locus can have more than one motif in the same class (e.g. CANVAS_RFC1's
# AAGGG and ACAGG are both pathogenic) -- cycling shades within a class keeps
# the classification legible while still telling individual motifs apart.
MOTIF_CLASS_SHADES = {
    "pathogenic": ["#E45756", "#8C2E2D", "#F2918F"],
    "benign": ["#54A24B", "#2E5C29", "#9FCB97"],
    "unknown": ["#999999", "#5C5C5C", "#C7C7C7"],
}

# Shown under the Metric/Locus selectors, updated whenever the metric changes.
METRIC_EXPLANATIONS = {
    "AL": "Total allele length, in base pairs (TRGT genotyper).",
    "MC": "Number of copies of each motif in the allele (TRGT genotyper).",
    "MS": "Number of base pairs in the allele spanned by each motif (TRGT genotyper).",
    "AP": "Purity: the fraction of the allele matching the "
          "expected motif sequence, ranging from 0 (no match) to 1 (perfect match) (TRGT genotyper).",
}


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_trid(info_field):
    match = re.search(r"(?:^|;)TRID=([^;]+)", info_field)
    return match.group(1) if match else None


def parse_motifs(info_field):
    """MOTIFS is an INFO field: one catalog of motifs per locus, shared by
    every sample/allele at that record (which motifs an allele actually
    contains is a separate, per-allele question -- see MC/MS)."""
    match = re.search(r"(?:^|;)MOTIFS=([^;]+)", info_field)
    return match.group(1).split(",") if match else []


def locus_title(locus, motifs, motif_classes):
    if not motifs:
        return locus
    annotated = [f"{m} [{c}]" for m, c in zip(motifs, motif_classes)]
    return f"{locus} ({', '.join(annotated)})"


def parse_format_labels(path):
    """Return {FORMAT id: description} from the ##FORMAT header lines."""
    labels = {}
    with open_vcf(path) as vcf:
        for line in vcf:
            if not line.startswith("##FORMAT"):
                if not line.startswith("#"):
                    break
                continue
            id_match = re.search(r"ID=([^,]+)", line)
            desc_match = re.search(r'Description="([^"]+)"', line)
            if id_match and desc_match:
                labels[id_match.group(1)] = desc_match.group(1)
    return labels


def parse_mc(token):
    """MC: integer count per motif, in MOTIFS order, e.g. "14_0" -> [14, 0]."""
    return [int(x) for x in token.split("_")]


def parse_ms(token, n_motifs):
    """MS: bp span per motif, keyed by the motif's own index, e.g.
    "0(0-27)_1(27-55)" -> motif 0 spans 27bp, motif 1 spans 28bp. A motif
    absent from this allele has no segment at all, so it defaults to 0."""
    spans = [0] * n_motifs
    for segment in token.split("_"):
        idx, start, end = re.search(r"(\d+)\((\d+)-(\d+)\)", segment).groups()
        spans[int(idx)] = int(end) - int(start)
    return spans


def parse_loci(path):
    """Return {locus_name: {"chrom": ..., "end": ..., "motifs": [...],
    "AL": [...], "AP": [...], "MC": [[per-motif values], ...],
    "MS": [[per-motif values], ...]}}."""
    loci = {}
    with open_vcf(path) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, info, fmt, samples = fields[0], fields[7], fields[8], fields[9:]
            locus = parse_trid(info)
            end = int(re.search(r"(?:^|;)END=(\d+)", info).group(1))
            motifs = parse_motifs(info)
            n_motifs = len(motifs)
            fmt_keys = fmt.split(":")
            al_idx, ap_idx = fmt_keys.index("AL"), fmt_keys.index("AP")
            mc_idx, ms_idx = fmt_keys.index("MC"), fmt_keys.index("MS")

            al, ap = [], []
            mc = [[] for _ in range(n_motifs)]
            ms = [[] for _ in range(n_motifs)]

            for sample in samples:
                sample_fields = sample.split(":")
                al_raw, ap_raw = sample_fields[al_idx], sample_fields[ap_idx]
                mc_raw, ms_raw = sample_fields[mc_idx], sample_fields[ms_idx]

                al.extend(int(x) for x in al_raw.split(",") if x != ".")
                ap.extend(float(x) for x in ap_raw.split(",") if x != ".")
                for token in mc_raw.split(","):
                    if token == ".":
                        continue
                    for m, count in enumerate(parse_mc(token)):
                        mc[m].append(count)
                for token in ms_raw.split(","):
                    if token == ".":
                        continue
                    for m, span in enumerate(parse_ms(token, n_motifs)):
                        ms[m].append(span)

            loci[locus] = {
                "chrom": chrom, "end": end, "motifs": motifs,
                "AL": al, "AP": ap, "MC": mc, "MS": ms,
            }
    return loci


def gene_name(locus):
    return locus.split("_", 1)[1] if "_" in locus else locus


def http_get_json(url):
    request = urllib.request.Request(url, headers={"User-Agent": "path-plots-script"})
    with urllib.request.urlopen(request) as response:
        return json.load(response)


def fetch_latest_strchive():
    """Always pull whatever STRchive release is current -- never pin a version."""
    release = http_get_json(STRCHIVE_LATEST_RELEASE_URL)
    version = release["tag_name"]
    asset = next(a for a in release["assets"] if re.fullmatch(r"STRchive-loci-.+\.json", a["name"]))
    records = http_get_json(asset["browser_download_url"])
    return version, records


def is_rotation(motif_a, motif_b):
    """STR motifs have no fixed reading frame: GAAAT and TGAAA are the same
    repeated unit read starting at a different position."""
    return len(motif_a) == len(motif_b) and motif_a in (motif_b + motif_b)


def classify_motif(motif, record):
    """Classify one of our VCF's tracked motifs as "pathogenic", "benign", or
    "unknown" per STRchive -- e.g. CANVAS_RFC1's AAGGG is pathogenic while
    ACAGG is also pathogenic-associated but AAAAG/AAAGGG (not tracked in our
    VCF) are benign. Matching is rotation-aware (see is_rotation)."""
    for tier, field in (
        ("pathogenic", "pathogenic_motif_reference_orientation"),
        ("benign", "benign_motif_reference_orientation"),
        ("unknown", "unknown_motif_reference_orientation"),
    ):
        if any(is_rotation(motif, m) for m in record.get(field) or []):
            return tier
    return "unknown"


def build_strchive_index(records):
    by_id = {r["id"]: r for r in records}
    by_coord = {(r["chrom"], r["stop_hg38"]): r for r in records if r.get("stop_hg38") is not None}
    by_gene = {}
    for r in records:
        by_gene.setdefault(r["gene"], []).append(r)
    return by_id, by_coord, by_gene


def match_strchive_record(trid, chrom, end, index):
    """Match a VCF locus to its STRchive record. Disease-code prefixes and
    exact coordinates can each drift a little between VCF annotation and the
    current STRchive release (e.g. our VCF's OPDM_ABCD3 is STRchive's
    OPDM5_ABCD3), so try, in order: exact TRID, exact coordinates, then the
    gene name -- picking the closest by coordinate if a gene has more than
    one STRchive locus (e.g. HOXA13's three separate repeats)."""
    by_id, by_coord, by_gene = index
    if trid in by_id:
        return by_id[trid]
    if (chrom, end) in by_coord:
        return by_coord[(chrom, end)]
    gene = gene_name(trid)
    base_gene = re.sub(r"-[IVX]+$", "", gene)
    candidates = [c for c in (by_gene.get(gene) or by_gene.get(base_gene) or []) if c["chrom"] == chrom]
    if not candidates:
        return None
    return min(candidates, key=lambda c: abs(c["stop_hg38"] - end))


OPEN_ENDED_MAX = 10 ** 7  # sentinel for a null upper bound; the plot's own x-range clips it


def strchive_tier(lo, hi):
    if lo is None and hi is None:
        return None
    return [lo if lo is not None else 0, hi if hi is not None else OPEN_ENDED_MAX]


def strchive_info_by_locus(loci, index):
    """{locus: {"motif_classes": [pathogenic/benign/unknown, ...] (always present,
    aligned with the locus's VCF motifs), "motif_len": N or None,
    "benign"/"intermediate"/"pathogenic": [lo, hi] or None}} for every VCF
    locus -- motif_classes defaults to all "unknown" when there's no
    matching STRchive record at all."""
    info = {}
    for locus, fields in loci.items():
        record = match_strchive_record(locus, fields["chrom"], fields["end"], index)
        if record is None:
            info[locus] = {
                "motif_classes": ["unknown"] * len(fields["motifs"]),
                "motif_len": None, "benign": None, "intermediate": None, "pathogenic": None,
            }
            continue
        info[locus] = {
            "motif_classes": [classify_motif(m, record) for m in fields["motifs"]],
            "motif_len": record["motif_len"],
            "benign": strchive_tier(record["benign_min"], record["benign_max"]),
            "intermediate": strchive_tier(record["intermediate_min"], record["intermediate_max"]),
            "pathogenic": strchive_tier(record["pathogenic_min"], record["pathogenic_max"]),
        }
    return info


def build_figure(loci, strchive_version):
    fig = go.Figure()
    locus_names = sorted(loci.keys(), key=gene_name)
    max_motifs = max(len(loci[locus]["motifs"]) for locus in locus_names)

    # One trace per (locus, motif slot); only the slot(s) relevant to the
    # currently selected locus+metric are ever made visible (by the JS
    # redraw() below). Content is a placeholder -- redraw() populates
    # everything, including the very first paint, so there's no risk of
    # Python's initial state and JS's dynamic state drifting apart.
    # Marker color is set on every redraw() (by pathogenic/benign/unknown
    # classification for MC/MS, neutral for AL/AP) -- the placeholder color
    # here is never actually seen.
    for _ in locus_names:
        for m in range(max_motifs):
            fig.add_trace(go.Histogram(
                x=[],
                marker=dict(
                    color="#cccccc",
                    line=dict(color="white", width=0.5),
                    opacity=0.75 if max_motifs > 1 else 1,
                ),
                visible=False,
            ))

    # Locus/Metric selection lives in a plain HTML <select> toolbar injected
    # by the JS below, not Plotly's native updatemenus dropdown: Plotly's
    # dropdown has no internal scroll region for 68 loci, so scrolling the
    # list scrolls the whole page along with it. A native <select>'s popup
    # scrolls independently, like any ordinary browser dropdown.
    fig.update_layout(
        annotations=[
            dict(
                text="", x=0.5, xanchor="center", y=1.18, yanchor="bottom", xref="paper", yref="paper",
                showarrow=False, font=dict(size=18),
            ),
            dict(
                text="", x=0.5, xanchor="center", y=-0.20, yanchor="top", xref="paper", yref="paper",
                showarrow=False, font=dict(size=12), align="left",
            ),
            dict(
                text=f"STRchive {strchive_version}",
                x=1.0, xanchor="right", y=-0.46, xref="paper", yref="paper",
                showarrow=False, font=dict(size=11, color="#888"),
            ),
        ],
        xaxis=dict(autorange=False),
        yaxis_title="Count",
        bargap=0.02,
        barmode="overlay",
        showlegend=False,
        template="plotly_white",
        width=900,
        height=520,
        margin=dict(t=150, b=190),
    )
    return fig, locus_names, max_motifs


def redraw_script(loci, locus_names, metric_labels, max_motifs, strchive_info):
    data = {
        locus: {
            "motifs": loci[locus]["motifs"],
            "AL": loci[locus]["AL"],
            "AP": loci[locus]["AP"],
            "MC": loci[locus]["MC"],
            "MS": loci[locus]["MS"],
        }
        for locus in locus_names
    }
    locus_titles = {
        locus: locus_title(locus, loci[locus]["motifs"], strchive_info[locus]["motif_classes"])
        for locus in locus_names
    }

    return """
    var gd = document.getElementById('{plot_id}');
    var DATA = %s;
    var LOCUS_NAMES = %s;
    var LOCUS_TITLES = %s;
    var METRICS = %s;
    var METRIC_LABELS = %s;
    var METRIC_EXPLANATIONS = %s;
    var INTEGER_METRICS = %s;
    var PER_MOTIF_METRICS = %s;
    var MAX_MOTIFS = %s;
    var TARGET_BINS = 30;
    var STRCHIVE = %s;
    var TIER_COLORS = %s;
    var MOTIF_CLASS_SHADES = %s;
    var OPEN_ENDED_MAX = %s;
    // STRchive's benign/intermediate/pathogenic ranges are repeat-unit
    // copy numbers. MC is already in those units; AL and MS are bp, so
    // their bands are approximated by scaling with the locus's motif length.
    var SHAPE_METRIC_SCALE = {AL: 'motif_len', MS: 'motif_len', MC: 1};

    var state = {locus: LOCUS_NAMES[0], metric: METRICS[0]};
    // Captured once, before any redraw touches the annotations array: index
    // 0 is the title, 1 the ranges-text block, 2 the static STRchive footer.
    var BASE_ANNOTATIONS = gd.layout.annotations.map(function(a) { return JSON.parse(JSON.stringify(a)); });

    function median(vals) {
        var sorted = vals.slice().sort(function(a, b) { return a - b; });
        var mid = Math.floor(sorted.length / 2);
        return (sorted.length & 1) ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
    }

    // Round a raw bin width to a "nice" step (1/2/5 x a power of ten), the
    // same convention d3/matplotlib use for readable histogram bins.
    function niceStep(rough) {
        if (rough <= 0) return null;
        var magnitude = Math.pow(10, Math.floor(Math.log10(rough)));
        var residual = rough / magnitude;
        var niceResidual = residual > 5 ? 10 : (residual > 2 ? 5 : (residual > 1 ? 2 : 1));
        return niceResidual * magnitude;
    }

    function computeBinSize(lo, hi, isInteger) {
        var step = niceStep((hi - lo) / TARGET_BINS);
        if (isInteger) {
            return Math.max(1, Math.round(step || 1));
        }
        return step || 0.01;
    }

    function computeShapes(locus, metric) {
        var ranges = STRCHIVE[locus];
        var scaleKey = SHAPE_METRIC_SCALE[metric];
        if (!ranges || ranges.motif_len == null || !scaleKey) return [];
        var scale = (scaleKey === 'motif_len') ? ranges.motif_len : scaleKey;
        var shapes = [];
        ['benign', 'intermediate', 'pathogenic'].forEach(function(tier) {
            var range = ranges[tier];
            if (!range) return;
            shapes.push({
                type: 'rect', xref: 'x', yref: 'paper',
                x0: range[0] * scale, x1: range[1] * scale, y0: 0, y1: 1,
                fillcolor: TIER_COLORS[tier], opacity: 0.15,
                line: {width: 0}, layer: 'below',
            });
        });
        return shapes;
    }

    function formatRange(range, scale) {
        if (!range) return 'not defined';
        var lo = range[0] * scale;
        if (range[1] >= OPEN_ENDED_MAX) return '≥' + lo;
        return lo + '–' + (range[1] * scale);
    }

    function rangesText(locus) {
        var s = STRCHIVE[locus];
        if (!s || s.motif_len == null) return 'No STRchive benign/intermediate/pathogenic ranges for this locus';
        var lines = ['benign', 'intermediate', 'pathogenic'].map(function(tier) {
            var r = s[tier];
            var label = tier.charAt(0).toUpperCase() + tier.slice(1);
            var swatch = '<span style="color:' + TIER_COLORS[tier] + '">■</span>';
            var text = r
                ? (formatRange(r, 1) + ' copies  /  ' + formatRange(r, s.motif_len) + ' bp')
                : 'not defined';
            return swatch + ' ' + label + ': ' + text;
        });
        return lines.join('<br>');
    }

    function redraw() {
        var locus = state.locus, metric = state.metric;
        var d = DATA[locus];
        var s = STRCHIVE[locus];
        var isPerMotif = PER_MOTIF_METRICS.indexOf(metric) !== -1;
        var n = isPerMotif ? d.motifs.length : 1;
        var base = LOCUS_NAMES.indexOf(locus) * MAX_MOTIFS;

        var visible = new Array(LOCUS_NAMES.length * MAX_MOTIFS).fill(false);
        var idx = [], xUpd = [], nameUpd = [], hoverUpd = [], colorUpd = [];
        var medianShapes = [], medianAnnotations = [];
        var lo = Infinity, hi = -Infinity;
        var classShadeCounts = {};

        for (var m = 0; m < n; m++) {
            var vals = isPerMotif ? d[metric][m] : d[metric];
            var color = '#4C78A8';
            if (isPerMotif) {
                var cls = s.motif_classes[m];
                var shadeIdx = classShadeCounts[cls] || 0;
                classShadeCounts[cls] = shadeIdx + 1;
                var shades = MOTIF_CLASS_SHADES[cls];
                color = shades[Math.min(shadeIdx, shades.length - 1)];
            }
            visible[base + m] = true;
            idx.push(base + m);
            xUpd.push(vals);
            var motifLabel = isPerMotif ? (d.motifs[m] + ' (' + s.motif_classes[m] + ')') : null;
            nameUpd.push(isPerMotif ? motifLabel : locus);
            colorUpd.push(color);
            var label = isPerMotif ? (motifLabel + ' ' + metric) : (METRIC_LABELS[metric] + ' (' + metric + ')');
            hoverUpd.push(label + ': %%{x}<br>Count: %%{y}<extra></extra>');
            for (var k = 0; k < vals.length; k++) {
                if (vals[k] < lo) lo = vals[k];
                if (vals[k] > hi) hi = vals[k];
            }

            var med = median(vals);
            medianShapes.push({
                type: 'line', xref: 'x', yref: 'paper',
                x0: med, x1: med, y0: 0, y1: 1,
                line: {color: color, width: 1.5, dash: 'dash'},
            });
            medianAnnotations.push({
                x: med, xref: 'x', y: 1.02, yref: 'paper', yanchor: 'bottom', xanchor: 'center',
                text: (isPerMotif ? d.motifs[m] + ' m' : 'M') + 'edian: ' + med,
                showarrow: false, font: {size: 11, color: color},
            });
        }

        var binSize = computeBinSize(lo, hi, INTEGER_METRICS.indexOf(metric) !== -1);
        var pad = Math.max((hi - lo) * 0.05, binSize);

        Plotly.restyle(gd, {visible: visible});
        Plotly.restyle(gd, {
            x: xUpd, name: nameUpd, hovertemplate: hoverUpd,
            'marker.color': colorUpd, 'xbins.size': binSize,
        }, idx);

        var titleAnnotation = JSON.parse(JSON.stringify(BASE_ANNOTATIONS[0]));
        titleAnnotation.text = LOCUS_TITLES[locus];
        var rangesAnnotation = JSON.parse(JSON.stringify(BASE_ANNOTATIONS[1]));
        rangesAnnotation.text = rangesText(locus);

        Plotly.relayout(gd, {
            'xaxis.range': [lo - pad, hi + pad],
            'xaxis.autorange': false,
            'xaxis.title.text': METRIC_LABELS[metric] + ' (' + metric + ')',
            'annotations': [titleAnnotation, rangesAnnotation, BASE_ANNOTATIONS[2]].concat(medianAnnotations),
            'showlegend': isPerMotif,
            'shapes': computeShapes(locus, metric).concat(medianShapes),
        });
        explanationEl.textContent = METRIC_EXPLANATIONS[metric] || '';
    }

    function addSelect(labelText, optionValues, optionLabels, initialValue, onChange) {
        var wrap = document.createElement('label');
        wrap.style.cssText = 'display:inline-flex; flex-direction:column; gap:4px; font:14px sans-serif; margin-right:24px;';
        var span = document.createElement('span');
        span.textContent = labelText;
        span.style.fontWeight = '600';
        var select = document.createElement('select');
        select.style.cssText = 'font-size:14px; padding:4px; max-width:340px;';
        optionValues.forEach(function(value, i) {
            var option = document.createElement('option');
            option.value = value;
            option.textContent = optionLabels[i];
            select.appendChild(option);
        });
        select.value = initialValue;
        select.addEventListener('change', function() { onChange(select.value); });
        wrap.appendChild(span);
        wrap.appendChild(select);
        return wrap;
    }

    var toolbar = document.createElement('div');
    toolbar.style.cssText = 'margin-bottom:8px;';
    toolbar.appendChild(addSelect(
        'Metric', METRICS, METRICS.map(function(m) { return m + ' (' + METRIC_LABELS[m] + ')'; }),
        state.metric, function(v) { state.metric = v; redraw(); },
    ));
    toolbar.appendChild(addSelect(
        'Locus', LOCUS_NAMES, LOCUS_NAMES,
        state.locus, function(v) { state.locus = v; redraw(); },
    ));
    gd.parentNode.insertBefore(toolbar, gd);

    var explanationEl = document.createElement('div');
    explanationEl.style.cssText = 'font:13px sans-serif; color:#444; max-width:900px; margin-bottom:16px; min-height:1.2em;';
    gd.parentNode.insertBefore(explanationEl, gd);

    redraw();
    """ % (
        json.dumps(data), json.dumps(locus_names), json.dumps(locus_titles),
        json.dumps(METRICS), json.dumps(metric_labels), json.dumps(METRIC_EXPLANATIONS), json.dumps(list(INTEGER_METRICS)),
        json.dumps(list(PER_MOTIF_METRICS)), json.dumps(max_motifs),
        json.dumps(strchive_info), json.dumps(TIER_COLORS), json.dumps(MOTIF_CLASS_SHADES), json.dumps(OPEN_ENDED_MAX),
    )


loci = parse_loci(vcf_path)
metric_labels = parse_format_labels(vcf_path)

strchive_version, strchive_records = fetch_latest_strchive()
print(f"Using STRchive {strchive_version} for benign/intermediate/pathogenic ranges")
strchive_index = build_strchive_index(strchive_records)
strchive_info = strchive_info_by_locus(loci, strchive_index)

fig, locus_names, max_motifs = build_figure(loci, strchive_version)

fig.write_html(
    "vcf_allele_histograms.html",
    include_plotlyjs="cdn",
    full_html=True,
    post_script=redraw_script(loci, locus_names, metric_labels, max_motifs, strchive_info),
)

print(f"Wrote vcf_allele_histograms.html with {len(loci)} loci x {len(METRICS)} metrics")
