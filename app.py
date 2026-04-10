"""
Streamlit app for primer design (cloning, deletion, protein tagging).

Wraps the existing CLI scripts via subprocess — no modifications to the
original programs.
"""

import json
import os
import re
import subprocess
import sys
import tempfile
import threading
import time
from datetime import date
from pathlib import Path

import pandas as pd
import streamlit as st
from Bio import SeqIO
from Bio.Restriction import CommOnly

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(page_title="Primer Designer", layout="wide")

# Tighten Streamlit's default top padding so the title sits near the top.
st.markdown(
    """
    <style>
    .block-container { padding-top: 2rem; }
    </style>
    """,
    unsafe_allow_html=True,
)

# NEB high-fidelity variants — map display name → BioPython base name.
# HF enzymes cut identically to the base enzyme, so we pass the base name
# to the script but let the user pick the HF variant in the UI.
HF_VARIANTS = {
    "AgeI-HF": "AgeI",
    "ApaLI-HF": "ApaLI",
    "AvrII-HF": "AvrII",
    "BamHI-HF": "BamHI",
    "BmtI-HF": "BmtI",
    "BsaI-HFv2": "BsaI",
    "BsiWI-HF": "BsiWI",
    "BsrGI-HF": "BsrGI",
    "BstEII-HF": "BstEII",
    "BstZ17I-HF": "BstZ17I",
    "DraIII-HF": "DraIII",
    "EagI-HF": "EagI",
    "EcoRI-HF": "EcoRI",
    "EcoRV-HF": "EcoRV",
    "HindIII-HF": "HindIII",
    "KpnI-HF": "KpnI",
    "MfeI-HF": "MfeI",
    "MluI-HF": "MluI",
    "NcoI-HF": "NcoI",
    "NheI-HF": "NheI",
    "NotI-HF": "NotI",
    "NruI-HF": "NruI",
    "NsiI-HF": "NsiI",
    "PstI-HF": "PstI",
    "PvuI-HF": "PvuI",
    "PvuII-HF": "PvuII",
    "SacI-HF": "SacI",
    "SalI-HF": "SalI",
    "SbfI-HF": "SbfI",
    "ScaI-HF": "ScaI",
    "SpeI-HF": "SpeI",
    "SphI-HF": "SphI",
    "SspI-HF": "SspI",
    "StyI-HF": "StyI",
}

# Searchable enzyme list: base enzymes from BioPython + HF variants.
_BASE_ENZYMES = {str(e) for e in CommOnly}
ENZYME_NAMES = sorted(
    _BASE_ENZYMES | {hf for hf, base in HF_VARIANTS.items() if base in _BASE_ENZYMES}
)


def _enzyme_index(name: str) -> int:
    try:
        return ENZYME_NAMES.index(name)
    except ValueError:
        return 0


def _resolve_enzyme(name: str) -> str:
    """Map an HF variant display name to the base enzyme name used by the script."""
    return HF_VARIANTS.get(name, name)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def save_uploaded(uploaded_file) -> str:
    """Write a Streamlit UploadedFile to a temp file and return its path."""
    suffix = Path(uploaded_file.name).suffix or ".fasta"
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, "wb") as f:
        f.write(uploaded_file.getvalue())
    return path


def run_script(script: str, args: list[str], cwd: Path) -> tuple[int, str, str]:
    """Run a primer design script via subprocess; return (returncode, csv_text, stderr)."""
    result = subprocess.run(
        [sys.executable, script] + args,
        capture_output=True, text=True, cwd=str(cwd),
    )
    return result.returncode, result.stdout, result.stderr


# Matches tqdm progress lines like "Designing primers:  50%|████| 5/10 [00:02<00:02]"
_TQDM_RE = re.compile(r"(\d+)%\|.*?\|\s*(\d+)/(\d+)")


def run_script_with_progress(
    script: str, args: list[str], cwd: Path, progress_bar, label: str = "Designing primers"
) -> tuple[int, str, str]:
    """Run a primer design script via subprocess, parsing tqdm stderr output to
    update a Streamlit progress bar. Returns (returncode, stdout, stderr_text)."""
    proc = subprocess.Popen(
        [sys.executable, script] + args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=str(cwd),
    )

    state = {"pct": 0, "n": 0, "total": 0}
    stderr_lines: list[str] = []

    def reader():
        buffer = b""
        while True:
            chunk = proc.stderr.read(64)
            if not chunk:
                break
            buffer += chunk
            # Split on both \n and \r so tqdm's in-place updates are parsed
            while True:
                i_n = buffer.find(b"\n")
                i_r = buffer.find(b"\r")
                if i_n == -1 and i_r == -1:
                    break
                if i_n == -1:
                    idx = i_r
                elif i_r == -1:
                    idx = i_n
                else:
                    idx = min(i_n, i_r)
                line_bytes = buffer[:idx]
                buffer = buffer[idx + 1 :]
                text = line_bytes.decode("utf-8", errors="replace")
                m = _TQDM_RE.search(text)
                if m:
                    state["pct"] = int(m.group(1))
                    state["n"] = int(m.group(2))
                    state["total"] = int(m.group(3))
                elif text.strip():
                    stderr_lines.append(text)
        # Flush remaining buffer
        if buffer.strip():
            stderr_lines.append(buffer.decode("utf-8", errors="replace"))

    thread = threading.Thread(target=reader, daemon=True)
    thread.start()

    while proc.poll() is None:
        total = state["total"]
        if total:
            text = f"{label}... {state['n']}/{total}"
        else:
            text = f"{label}..."
        progress_bar.progress(state["pct"] / 100.0, text=text)
        time.sleep(0.1)

    thread.join()
    # Final update
    progress_bar.progress(1.0, text=f"{label}... done")

    stdout = proc.stdout.read().decode("utf-8", errors="replace") if proc.stdout else ""
    return proc.returncode, stdout, "\n".join(stderr_lines)


def resolve_plasmid(selected_name: str | None, uploaded_file) -> str | None:
    """Return plasmid path: uploaded file if provided, else local selection."""
    if uploaded_file:
        return save_uploaded(uploaded_file)
    if selected_name:
        return str(SCRIPT_DIR / selected_name)
    return None


def resolve_genome(uploaded_files: list) -> list[str] | None:
    """Return genome paths: uploaded files if provided, else local defaults."""
    if uploaded_files:
        return [save_uploaded(g) for g in uploaded_files]
    if DEFAULT_GENOME:
        return [str(SCRIPT_DIR / f) for f in DEFAULT_GENOME]
    return None


def resolve_genes(uploaded_file) -> str | None:
    """Return genes path: uploaded file if provided, else local default."""
    if uploaded_file:
        return save_uploaded(uploaded_file)
    if DEFAULT_GENES:
        return str(SCRIPT_DIR / DEFAULT_GENES)
    return None


@st.cache_data(show_spinner=False)
def _parse_gene_ids_from_path(path: str) -> frozenset:
    return frozenset(rec.id for rec in SeqIO.parse(path, "fasta"))


@st.cache_data(show_spinner=False)
def _parse_gene_ids_from_bytes(data: bytes) -> frozenset:
    from io import StringIO
    return frozenset(rec.id for rec in SeqIO.parse(StringIO(data.decode("utf-8")), "fasta"))


def get_current_gene_ids(uploaded_file) -> frozenset:
    """Return the set of gene IDs in the currently-configured genes FASTA."""
    if uploaded_file is not None:
        try:
            return _parse_gene_ids_from_bytes(uploaded_file.getvalue())
        except Exception:
            return frozenset()
    if DEFAULT_GENES:
        try:
            return _parse_gene_ids_from_path(str(SCRIPT_DIR / DEFAULT_GENES))
        except Exception:
            return frozenset()
    return frozenset()


def _date_stamp() -> str:
    """Return today's date as mmddyy."""
    return date.today().strftime("%m%d%y")


def _safe_part(s: str) -> str:
    """Sanitize a string for use as a filename component."""
    return re.sub(r"[^\w.-]", "_", s)


def _plasmid_stem(plasmid_name, plasmid_upload) -> str:
    """Return the plasmid name without extension, preferring uploaded file."""
    if plasmid_upload is not None:
        return Path(plasmid_upload.name).stem
    if plasmid_name:
        return Path(plasmid_name).stem
    return "plasmid"


def show_stderr(stderr: str):
    """Parse stderr lines and display warnings / errors only."""
    for line in stderr.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.upper().startswith("WARNING"):
            st.warning(line)
        elif line.startswith("Error"):
            st.error(line)


def show_results(
    csv_path: str,
    stderr: str,
    download_filename: str | None = None,
    rename_columns: dict[str, str] | None = None,
    extra_downloads: list[dict] | None = None,
):
    """Read the output CSV, display it, and offer download buttons.

    ``extra_downloads`` is an optional list of dicts describing additional
    download buttons. Each dict must contain ``label``, ``data``,
    ``file_name``, ``mime`` and may optionally contain ``help``. When
    ``extra_downloads`` is non-empty, the CSV download button is suppressed
    entirely and only the extras are shown (centered in their own row);
    otherwise the CSV download button is shown on its own.
    """
    show_stderr(stderr)

    p = Path(csv_path)
    if not p.exists() or p.stat().st_size == 0:
        st.error("No output was produced. Check the warnings above.")
        return

    raw = p.read_text()

    if rename_columns:
        # Rewrite just the header row with renamed columns; leave the rest
        # of the CSV byte-for-byte intact so multi-row layouts (e.g. the
        # protein-tag Avg row) survive.
        import csv as _csv
        from io import StringIO as _StringIO

        first_newline = raw.find("\n")
        if first_newline != -1:
            header_line = raw[:first_newline]
            rest = raw[first_newline + 1:]
            header = next(_csv.reader([header_line]))
            new_header = [rename_columns.get(c, c) for c in header]
            buf = _StringIO()
            _csv.writer(buf).writerow(new_header)
            raw = buf.getvalue().rstrip("\r\n") + "\n" + rest

    _extras = extra_downloads or []
    if _extras:
        # Single-gene assembly mode: the .dna (and optionally .gbk) file
        # carries everything the user needs, so we suppress the CSV download
        # and only show the extras, centered in their own row. Extras get
        # wider columns so their long labels ("Download SnapGene file
        # (.dna)", "Download annotated GenBank (.gbk)") don't wrap.
        _cols = st.columns(
            [1] + [3] * len(_extras) + [1], gap="small"
        )
        for _col, _btn in zip(_cols[1:-1], _extras):
            with _col:
                st.download_button(
                    label=_btn["label"],
                    data=_btn["data"],
                    file_name=_btn["file_name"],
                    mime=_btn["mime"],
                    help=_btn.get("help"),
                    use_container_width=True,
                )
    else:
        st.download_button(
            label="Download CSV",
            data=raw,
            file_name=download_filename or p.name,
            mime="text/csv",
        )

    try:
        df = pd.read_csv(csv_path)
        df = df.dropna(how="all")
        # Replace all NaN/NA with empty strings up front so nothing can render
        # as "nan" or "<NA>". Numeric columns become object dtype after this,
        # so the formatters below use try/except to handle both numbers and
        # blanks / non-numeric labels (like the protein-tag "Avg" row).
        df = df.fillna("")
        if rename_columns:
            df = df.rename(columns=rename_columns)

        def _fmt_int(v):
            try:
                return f"{int(float(v))}"
            except (ValueError, TypeError):
                return str(v)

        def _fmt_float(digits):
            def _format(v):
                try:
                    return f"{float(v):.{digits}f}"
                except (ValueError, TypeError):
                    return str(v)
            return _format

        int_cols = [
            c for c in df.columns
            if any(k in c.lower() for k in ("length", "_bp", "start", "end"))
        ]
        tm_cols = [c for c in df.columns if "tm" in c.lower()]
        penalty_cols = [c for c in df.columns if "penalty" in c.lower()]

        fmt = {c: _fmt_int for c in int_cols}
        fmt.update({c: _fmt_float(1) for c in tm_cols})
        fmt.update({c: _fmt_float(3) for c in penalty_cols})

        styled = (
            df.style
            .format(fmt)
            .set_table_styles([
                {"selector": "th", "props": [("white-space", "nowrap")]},
                {"selector": "td", "props": [("white-space", "nowrap")]},
            ])
        )
        st.dataframe(styled, hide_index=True, use_container_width=True)
    except Exception:
        st.code(raw)


# ---------------------------------------------------------------------------
# Scripts
# ---------------------------------------------------------------------------

SCRIPTS = {
    "Cloning primers": "design_cloning_primers_2.0.py",
    "Deletion primers": "design_deletion_primers.py",
    "Protein tag primers": "design_protein_tag_primers.py",
}

# ---------------------------------------------------------------------------
# Script directory (where the .py files live)
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).parent

# ---------------------------------------------------------------------------
# Auto-detect local genome and genes files
# ---------------------------------------------------------------------------

DEFAULT_GENOME = sorted(
    p.name for p in SCRIPT_DIR.glob("chr*.fasta")
)
DEFAULT_GENES = "genes.fasta" if (SCRIPT_DIR / "genes.fasta").exists() else None
LOCAL_PLASMIDS = sorted(
    p.name for p in SCRIPT_DIR.glob("p*.fasta")
)

# ---------------------------------------------------------------------------
# Main area
# ---------------------------------------------------------------------------

_title_col, _issue_col = st.columns([5, 1])
with _title_col:
    st.title("Primer Designer")
with _issue_col:
    st.markdown(
        '<div style="text-align: right; padding-top: 1.8rem;">'
        '<a href="mailto:julek@batory.pl?subject=Primers%20app%20feedback" '
        'style="font-size: 0.9rem;">Report an issue</a>'
        '</div>',
        unsafe_allow_html=True,
    )

tab_cloning, tab_deletion, tab_tag = st.tabs(["Cloning", "Deletion", "Protein Tag"])

# ========================== CLONING ========================================
with tab_cloning:
    st.caption(
        "Designs a forward + reverse primer pair per gene for amplifying genes "
        "and inserting them into a vector via HiFi assembly or restriction cloning."
    )

    st.subheader("Input files")
    with st.expander("Genome and genes (default: Vibrio cholerae O1 El Tor A1552)"):
        genome_files = st.file_uploader("Genome FASTA(s) — upload to override", type=["fasta", "fa", "fna"], accept_multiple_files=True, key="cloning_genome")
        genes_file = st.file_uploader("Genes FASTA — upload to override", type=["fasta", "fa", "fna"], key="cloning_genes")
    col1, col2, _spacer_plasmid = st.columns([1, 3, 4])
    with col1:
        _pmm_idx = LOCAL_PLASMIDS.index("pMMB67EH.fasta") if "pMMB67EH.fasta" in LOCAL_PLASMIDS else 0
        plasmid_name = st.selectbox("Plasmid", LOCAL_PLASMIDS, index=_pmm_idx, format_func=lambda p: Path(p).stem, key="cloning_plasmid_sel") if LOCAL_PLASMIDS else None
    with col2:
        plasmid_upload = st.file_uploader("Or upload a plasmid FASTA", type=["fasta", "fa", "fna"], key="cloning_plasmid_up")

    st.subheader("Parameters")
    col1, col2, _spacer = st.columns([1, 1, 2])
    with col1:
        five_prime_enzyme = st.selectbox(
            "5' End",
            ENZYME_NAMES,
            index=_enzyme_index("BamHI-HF"),
            key="cloning_five_prime_enzyme",
            help="Enzyme at the vector backbone's 5' end (where the insert's 3' end attaches). Matches NEBuilder's 5' enzyme.",
        )
    with col2:
        three_prime_enzyme = st.selectbox(
            "3' End",
            ENZYME_NAMES,
            index=_enzyme_index("BamHI-HF"),
            key="cloning_three_prime_enzyme",
            help="Enzyme at the vector backbone's 3' end (where the insert's 5' end attaches). Matches NEBuilder's 3' enzyme.",
        )

    col1, col2, _spacer = st.columns([1, 1, 2])
    with col1:
        tail_mode = st.selectbox(
            "Tail mode",
            ["plasmid_overlaps", "restriction_sites"],
            key="cloning_tail_mode",
            help=(
                "**plasmid_overlaps**: primer tails are homologous to the plasmid "
                "sequences flanking the cut sites (for HiFi/Gibson assembly).\n\n"
                "**restriction_sites**: tails contain the restriction enzyme "
                "recognition sites for classic restriction cloning."
            ),
        )
    with col2:
        replace_arc = st.selectbox(
            "Replace arc",
            ["shorter_arc", "longer_arc", "three_to_five", "five_to_three"],
            index=0,
            disabled=(tail_mode != "plasmid_overlaps"),
            key="cloning_replace_arc",
            help=(
                "For circular plasmids with two distinct cut sites, chooses which "
                "arc between the cuts is replaced by the insert.\n\n"
                "**shorter_arc** / **longer_arc**: pick by length (usually what you want — "
                "shorter_arc replaces the small MCS fragment).\n\n"
                "**three_to_five** / **five_to_three**: pick by direction regardless of "
                "length (escape hatch for unusual cases).\n\n"
                "Ignored for single-cut or linear plasmids."
            ),
        )

    col1, col2, _spacer = st.columns([1, 3, 4])
    with col1:
        overlap_length = st.number_input(
            "Overlap length (bp)",
            value=20,
            min_value=0,
            key="cloning_overlap",
            help="Length of the 5' plasmid overlap tail added to each primer (for HiFi assembly into the vector).",
        )
    with col2:
        gene_ids = st.text_input("Gene IDs (space-separated, leave blank for all)", value="", key="cloning_gene_ids")

    # Restriction-site mode extras
    if tail_mode == "restriction_sites":
        col1, col2, col3, col4, _spacer = st.columns([1, 1, 1, 1, 4])
        with col1:
            five_prime_clamp = st.text_input("5' clamp", value="GCGC", key="cloning_five_prime_clamp")
        with col2:
            three_prime_clamp = st.text_input("3' clamp", value="GCGC", key="cloning_three_prime_clamp")
        with col3:
            five_prime_extra = st.text_input("5' extra", value="", key="cloning_five_prime_extra")
        with col4:
            three_prime_extra = st.text_input("3' extra", value="", key="cloning_three_prime_extra")
    else:
        five_prime_clamp = "GCGC"
        three_prime_clamp = "GCGC"
        five_prime_extra = ""
        three_prime_extra = ""

    with st.expander("Advanced parameters"):
        col1, col2, col3 = st.columns(3)
        with col1:
            opt_primer_size = st.number_input("Optimal primer size", value=22, min_value=10, key="cloning_opt_size")
        with col2:
            min_primer_size = st.number_input("Min primer size", value=18, min_value=10, key="cloning_min_size")
        with col3:
            max_primer_size = st.number_input("Max primer size", value=28, min_value=10, key="cloning_max_size")

        col1, col2, col3 = st.columns(3)
        with col1:
            opt_tm = st.number_input("Optimal Tm", value=60.0, key="cloning_opt_tm")
        with col2:
            min_tm = st.number_input("Min Tm", value=57.0, key="cloning_min_tm")
        with col3:
            max_tm = st.number_input("Max Tm", value=63.0, key="cloning_max_tm")

        col1, col2, col3 = st.columns(3)
        with col1:
            mv_conc = st.number_input("Monovalent cation conc (mM)", value=500.0, key="cloning_mv_conc")
        with col2:
            gc_clamp = st.number_input("GC clamp", value=1, min_value=0, key="cloning_gc_clamp")
        with col3:
            min_gc = st.number_input("Min GC%", value=35.0, key="cloning_min_gc")

        max_gc = st.number_input("Max GC%", value=65.0, key="cloning_max_gc")
        linear_plasmid = st.checkbox("Linear plasmid", value=False, key="cloning_linear")
        allow_unmatched = st.checkbox("Allow unmatched genes", value=False, key="cloning_allow_unmatched")

    if st.button("Design primers", type="primary", use_container_width=True, key="cloning_run"):
        plasmid_path = resolve_plasmid(plasmid_name, plasmid_upload)
        genome_paths = resolve_genome(genome_files)
        genes_path = resolve_genes(genes_file)
        if not plasmid_path:
            st.error("Select or upload a plasmid FASTA.")
        elif not genome_paths:
            st.error("No genome files found. Upload genome FASTA(s).")
        elif not genes_path:
            st.error("No genes file found. Upload a genes FASTA.")
        elif not five_prime_enzyme or not three_prime_enzyme:
            st.error("Both enzymes are required.")
        else:
            fd, out_path = tempfile.mkstemp(suffix=".csv")
            os.close(fd)

            # Only request the assembled plasmid outputs when exactly one gene
            # is being designed (they are single-assembly views).
            _gene_id_list = gene_ids.split()
            gbk_path: str | None = None
            dna_path: str | None = None
            if len(_gene_id_list) == 1 and tail_mode == "plasmid_overlaps" and not linear_plasmid:
                gbk_fd, gbk_path = tempfile.mkstemp(suffix=".gbk")
                os.close(gbk_fd)
                dna_fd, dna_path = tempfile.mkstemp(suffix=".dna")
                os.close(dna_fd)

            args = [
                "--plasmid", plasmid_path,
                "--genome", *genome_paths,
                "--genes", genes_path,
                "--output", out_path,
                "--three-prime-enzyme", _resolve_enzyme(three_prime_enzyme),
                "--five-prime-enzyme", _resolve_enzyme(five_prime_enzyme),
                "--tail-mode", tail_mode,
                "--replace-arc", replace_arc,
                "--overlap-length", str(overlap_length),
                "--opt-primer-size", str(opt_primer_size),
                "--min-primer-size", str(min_primer_size),
                "--max-primer-size", str(max_primer_size),
                "--opt-tm", str(opt_tm),
                "--min-tm", str(min_tm),
                "--max-tm", str(max_tm),
                "--min-gc", str(min_gc),
                "--max-gc", str(max_gc),
                "--mv-conc", str(mv_conc),
                "--gc-clamp", str(gc_clamp),
                "--three-prime-clamp", three_prime_clamp,
                "--five-prime-clamp", five_prime_clamp,
                "--three-prime-extra", three_prime_extra,
                "--five-prime-extra", five_prime_extra,
            ]
            if gbk_path:
                args.extend(["--gbk-output", gbk_path])
            if dna_path:
                args.extend(["--dna-output", dna_path])
            if linear_plasmid:
                args.append("--linear-plasmid")
            if allow_unmatched:
                args.append("--allow-unmatched-genes")
            if gene_ids.strip():
                args.extend(["--gene-ids"] + gene_ids.strip().split())

            progress_bar = st.progress(0.0, text="Designing primers...")
            rc, _, stderr = run_script_with_progress(
                "design_cloning_primers_2.0.py", args, SCRIPT_DIR, progress_bar
            )
            progress_bar.empty()
            if rc != 0 and not Path(out_path).exists():
                st.error(f"Script exited with code {rc}")
                show_stderr(stderr)
            else:
                download_name = (
                    f"cloning_{_safe_part(_plasmid_stem(plasmid_name, plasmid_upload))}"
                    f"_{_safe_part(five_prime_enzyme)}_{_safe_part(three_prime_enzyme)}"
                    f"_{_date_stamp()}.csv"
                )
                rename_map = None
                if len(_gene_id_list) == 1:
                    _gid = _gene_id_list[0]
                    rename_map = {
                        "forward_primer_full_5to3": f"{_gid}_fwd",
                        "reverse_primer_full_5to3": f"{_gid}_rev",
                    }

                # Collect the assembled-plasmid download buttons (single gene
                # only) to display next to the CSV download button.
                _extra_downloads: list[dict] = []
                if dna_path and Path(dna_path).exists() and Path(dna_path).stat().st_size > 0:
                    dna_download_name = (
                        f"{_safe_part(_plasmid_stem(plasmid_name, plasmid_upload))}"
                        f"_{_safe_part(_gene_id_list[0])}"
                        f"_{_date_stamp()}.dna"
                    )
                    with open(dna_path, "rb") as _df:
                        _dna_bytes = _df.read()
                    _extra_downloads.append({
                        "label": "Download SnapGene file (.dna)",
                        "data": _dna_bytes,
                        "file_name": dna_download_name,
                        "mime": "application/octet-stream",
                        "help": "Native SnapGene file — primers land in the Primers panel as proper arrows (not feature rectangles).",
                    })
                if gbk_path and Path(gbk_path).exists() and Path(gbk_path).stat().st_size > 0:
                    gbk_download_name = (
                        f"{_safe_part(_plasmid_stem(plasmid_name, plasmid_upload))}"
                        f"_{_safe_part(_gene_id_list[0])}"
                        f"_{_date_stamp()}.gbk"
                    )
                    with open(gbk_path, "rb") as _gf:
                        _gbk_bytes = _gf.read()
                    _extra_downloads.append({
                        "label": "Download annotated GenBank (.gbk)",
                        "data": _gbk_bytes,
                        "file_name": gbk_download_name,
                        "mime": "chemical/seq-na-genbank",
                        "help": "Open in Benchling or any GenBank viewer to see the assembled plasmid.",
                    })

                show_results(
                    out_path, stderr,
                    download_filename=download_name,
                    rename_columns=rename_map,
                    extra_downloads=_extra_downloads,
                )

# ========================== DELETION =======================================
with tab_deletion:
    st.caption(
        "Designs four primers per gene (AB_fwd, AB_rev, CD_fwd, CD_rev) across two amplicons "
        "for in-frame chromosomal deletion via HiFi assembly."
    )

    st.subheader("Input files")
    with st.expander("Genome and genes (default: Vibrio cholerae O1 El Tor A1552)"):
        genome_files = st.file_uploader("Genome FASTA(s) — upload to override", type=["fasta", "fa", "fna"], accept_multiple_files=True, key="deletion_genome")
        genes_file = st.file_uploader("Genes FASTA — upload to override", type=["fasta", "fa", "fna"], key="deletion_genes")
    col1, col2, _spacer_plasmid = st.columns([1, 3, 4])
    with col1:
        plasmid_name = st.selectbox("Plasmid", LOCAL_PLASMIDS, format_func=lambda p: Path(p).stem, key="deletion_plasmid_sel") if LOCAL_PLASMIDS else None
    with col2:
        plasmid_upload = st.file_uploader("Or upload a plasmid FASTA", type=["fasta", "fa", "fna"], key="deletion_plasmid_up")

    st.subheader("Parameters")
    col1, col2, _spacer = st.columns([1, 1, 2])
    with col1:
        five_prime_enzyme = st.selectbox(
            "5' End",
            ENZYME_NAMES,
            index=_enzyme_index("SacI-HF"),
            key="deletion_five_prime_enzyme",
            help="Enzyme at the vector backbone's 5' end (where the insert's 3' end attaches). Matches NEBuilder's 5' enzyme.",
        )
    with col2:
        three_prime_enzyme = st.selectbox(
            "3' End",
            ENZYME_NAMES,
            index=_enzyme_index("NcoI-HF"),
            key="deletion_three_prime_enzyme",
            help="Enzyme at the vector backbone's 3' end (where the insert's 5' end attaches). Matches NEBuilder's 3' enzyme.",
        )

    col1, col2, col3, _spacer = st.columns([1, 1, 1.3, 3.7])
    with col1:
        overlap_length = st.number_input(
            "Vector overlap (bp)",
            value=20,
            min_value=0,
            key="deletion_overlap",
            help="Length of the vector overlap tail prepended to Primers A and D (for HiFi assembly into the digested vector).",
        )
    with col2:
        junction_overlap = st.number_input(
            "Junction overlap (bp)",
            value=10,
            min_value=0,
            key="deletion_junction",
            help="Length of the AB↔CD junction overlap tail prepended to Primers B and C so the two amplicons stitch together during HiFi assembly.",
        )
    with col3:
        flank_length = st.text_input(
            "Flank length (blank = auto)",
            value="",
            help="Auto-scaled by gene length: <1500 bp -> 509 | 1500-3000 -> 709 | >3000 -> 909",
            key="deletion_flank",
        )

    col1, _spacer = st.columns([2, 2])
    with col1:
        gene_ids = st.text_input("Gene IDs (space-separated, leave blank for all)", value="", key="deletion_gene_ids")

    with st.expander("Advanced parameters"):
        col1, col2, col3 = st.columns(3)
        with col1:
            opt_primer_size = st.number_input("Optimal primer size", value=20, min_value=10, key="deletion_opt_size")
        with col2:
            min_primer_size = st.number_input("Min primer size", value=18, min_value=10, key="deletion_min_size")
        with col3:
            max_primer_size = st.number_input("Max primer size", value=28, min_value=10, key="deletion_max_size")

        col1, col2, col3 = st.columns(3)
        with col1:
            opt_tm = st.number_input("Optimal Tm", value=60.0, key="deletion_opt_tm")
        with col2:
            mv_conc = st.number_input("Monovalent cation conc (mM)", value=500.0, key="deletion_mv_conc")
        with col3:
            gc_clamp = st.number_input("GC clamp", value=1, min_value=0, key="deletion_gc_clamp")

        linear_plasmid = st.checkbox("Linear plasmid", value=False, key="deletion_linear")
        allow_unmatched = st.checkbox("Allow unmatched genes", value=False, key="deletion_allow_unmatched")

    if st.button("Design primers", type="primary", use_container_width=True, key="deletion_run"):
        plasmid_path = resolve_plasmid(plasmid_name, plasmid_upload)
        genome_paths = resolve_genome(genome_files)
        genes_path = resolve_genes(genes_file)
        if not plasmid_path:
            st.error("Select or upload a plasmid FASTA.")
        elif not genome_paths:
            st.error("No genome files found. Upload genome FASTA(s).")
        elif not genes_path:
            st.error("No genes file found. Upload a genes FASTA.")
        elif not five_prime_enzyme or not three_prime_enzyme:
            st.error("Both enzymes are required.")
        else:
            fd, out_path = tempfile.mkstemp(suffix=".csv")
            os.close(fd)

            # Only request the assembled deletion plasmid outputs when exactly
            # one gene is being designed against a circular plasmid (they are
            # single-assembly views).
            _gene_id_list = gene_ids.split()
            dna_path: str | None = None
            gbk_path: str | None = None
            if len(_gene_id_list) == 1 and not linear_plasmid:
                dna_fd, dna_path = tempfile.mkstemp(suffix=".dna")
                os.close(dna_fd)
                gbk_fd, gbk_path = tempfile.mkstemp(suffix=".gbk")
                os.close(gbk_fd)

            args = [
                "--plasmid", plasmid_path,
                "--genome", *genome_paths,
                "--genes", genes_path,
                "--output", out_path,
                "--three-prime-enzyme", _resolve_enzyme(three_prime_enzyme),
                "--five-prime-enzyme", _resolve_enzyme(five_prime_enzyme),
                "--overlap-length", str(overlap_length),
                "--junction-overlap", str(junction_overlap),
                "--opt-primer-size", str(opt_primer_size),
                "--min-primer-size", str(min_primer_size),
                "--max-primer-size", str(max_primer_size),
                "--opt-tm", str(opt_tm),
                "--mv-conc", str(mv_conc),
                "--gc-clamp", str(gc_clamp),
            ]
            if dna_path:
                args.extend(["--dna-output", dna_path])
            if gbk_path:
                args.extend(["--gbk-output", gbk_path])
            if flank_length.strip():
                args.extend(["--flank-length", flank_length.strip()])
            if linear_plasmid:
                args.append("--linear-plasmid")
            if allow_unmatched:
                args.append("--allow-unmatched-genes")
            if gene_ids.strip():
                args.extend(["--gene-ids"] + gene_ids.strip().split())

            progress_bar = st.progress(0.0, text="Designing primers...")
            rc, _, stderr = run_script_with_progress(
                "design_deletion_primers.py", args, SCRIPT_DIR, progress_bar
            )
            progress_bar.empty()
            if rc != 0 and not Path(out_path).exists():
                st.error(f"Script exited with code {rc}")
                show_stderr(stderr)
            else:
                download_name = (
                    f"deletion_{_safe_part(_plasmid_stem(plasmid_name, plasmid_upload))}"
                    f"_{_safe_part(five_prime_enzyme)}_{_safe_part(three_prime_enzyme)}"
                    f"_{_date_stamp()}.csv"
                )
                rename_map = None
                if len(_gene_id_list) == 1:
                    _gid = _gene_id_list[0]
                    rename_map = {
                        "AB_fwd": f"{_gid}_AB_fwd",
                        "AB_rev": f"{_gid}_AB_rev",
                        "CD_fwd": f"{_gid}_CD_fwd",
                        "CD_rev": f"{_gid}_CD_rev",
                    }

                # Collect the assembled deletion plasmid download buttons
                # (single gene only) to display next to the CSV download.
                _extra_downloads: list[dict] = []
                if dna_path and Path(dna_path).exists() and Path(dna_path).stat().st_size > 0:
                    dna_download_name = (
                        f"{_safe_part(_plasmid_stem(plasmid_name, plasmid_upload))}"
                        f"_delta_{_safe_part(_gene_id_list[0])}"
                        f"_{_date_stamp()}.dna"
                    )
                    with open(dna_path, "rb") as _df:
                        _dna_bytes = _df.read()
                    _extra_downloads.append({
                        "label": "Download SnapGene file (.dna)",
                        "data": _dna_bytes,
                        "file_name": dna_download_name,
                        "mime": "application/octet-stream",
                        "help": "Native SnapGene file — deletion primers land in the Primers panel as proper arrows (not feature rectangles).",
                    })
                if gbk_path and Path(gbk_path).exists() and Path(gbk_path).stat().st_size > 0:
                    gbk_download_name = (
                        f"{_safe_part(_plasmid_stem(plasmid_name, plasmid_upload))}"
                        f"_delta_{_safe_part(_gene_id_list[0])}"
                        f"_{_date_stamp()}.gbk"
                    )
                    with open(gbk_path, "rb") as _gf:
                        _gbk_bytes = _gf.read()
                    _extra_downloads.append({
                        "label": "Download annotated GenBank (.gbk)",
                        "data": _gbk_bytes,
                        "file_name": gbk_download_name,
                        "mime": "chemical/seq-na-genbank",
                        "help": "Open in Benchling or any GenBank viewer to see the assembled deletion plasmid.",
                    })

                show_results(
                    out_path, stderr,
                    download_filename=download_name,
                    rename_columns=rename_map,
                    extra_downloads=_extra_downloads,
                )

# ========================== PROTEIN TAG ====================================
with tab_tag:
    st.caption(
        "Designs six primers across three amplicons to fuse a protein tag to the "
        "N- or C-terminus of a single gene via HiFi/Gibson assembly."
    )

    st.subheader("Input files")
    with st.expander("Genome and genes (default: Vibrio cholerae O1 El Tor A1552)"):
        genome_files = st.file_uploader("Genome FASTA(s) — upload to override", type=["fasta", "fa", "fna"], accept_multiple_files=True, key="tag_genome")
        genes_file = st.file_uploader("Genes FASTA — upload to override", type=["fasta", "fa", "fna"], key="tag_genes")
    col1, col2, _spacer_plasmid = st.columns([1, 3, 4])
    with col1:
        plasmid_name = st.selectbox("Plasmid", LOCAL_PLASMIDS, format_func=lambda p: Path(p).stem, key="tag_plasmid_sel") if LOCAL_PLASMIDS else None
    with col2:
        plasmid_upload = st.file_uploader("Or upload a plasmid FASTA", type=["fasta", "fa", "fna"], key="tag_plasmid_up")

    st.subheader("Parameters")
    col1, col2, col3, _spacer = st.columns([1, 1, 1, 4])
    with col1:
        terminus = st.radio(
            "Terminus",
            ["C", "N"],
            index=0,
            horizontal=True,
            key="tag_terminus",
        )
    with col2:
        five_prime_enzyme = st.selectbox(
            "5' End",
            ENZYME_NAMES,
            index=_enzyme_index("SacI-HF"),
            key="tag_five_prime_enzyme",
            help="Enzyme at the vector backbone's 5' end (where the insert's 3' end attaches). Matches NEBuilder's 5' enzyme.",
        )
    with col3:
        three_prime_enzyme = st.selectbox(
            "3' End",
            ENZYME_NAMES,
            index=_enzyme_index("NcoI-HF"),
            key="tag_three_prime_enzyme",
            help="Enzyme at the vector backbone's 3' end (where the insert's 5' end attaches). Matches NEBuilder's 3' enzyme.",
        )

    # Tag selection
    if terminus == "C":
        tag_options = ["GGGGG_GFP", "GGSS_Halo", "Custom"]
    else:
        tag_options = ["Custom"]

    tag_label = "Linker + Tag" if terminus == "C" else "Tag + Linker"
    col_tag, col_gene, _spacer_tag = st.columns([1, 2, 4])
    with col_tag:
        tag_choice = st.selectbox(tag_label, tag_options, key="tag_choice")
    with col_gene:
        gene_id = st.text_input("Gene ID (single, required)", value="", key="tag_gene_id")

    _current_gene_ids = get_current_gene_ids(genes_file)
    _valid_ids_json = json.dumps(list(_current_gene_ids))
    st.components.v1.html(
        f"""
        <script>
        (function() {{
            const VALID = new Set({_valid_ids_json});
            const STYLE_ID = 'tag-gene-id-validator-style';
            const doc = window.parent.document;

            if (!doc.getElementById(STYLE_ID)) {{
                const style = doc.createElement('style');
                style.id = STYLE_ID;
                style.textContent = `
                    .st-key-tag_gene_id.invalid-gene [data-baseweb="input"],
                    .st-key-tag_gene_id.invalid-gene [data-baseweb="base-input"] {{
                        border-color: rgb(255, 75, 75) !important;
                    }}
                    .st-key-tag_gene_id.invalid-gene label,
                    .st-key-tag_gene_id.invalid-gene label p {{
                        color: rgb(255, 75, 75) !important;
                    }}
                    .st-key-tag_gene_id.valid-gene [data-baseweb="input"],
                    .st-key-tag_gene_id.valid-gene [data-baseweb="base-input"] {{
                        border-color: rgba(49, 51, 63, 0.2) !important;
                        box-shadow: none !important;
                    }}
                `;
                doc.head.appendChild(style);
            }}

            function check(input) {{
                const wrapper = input.closest('.st-key-tag_gene_id');
                if (!wrapper) return;
                if (VALID.has(input.value.trim())) {{
                    wrapper.classList.remove('invalid-gene');
                    wrapper.classList.add('valid-gene');
                }} else {{
                    wrapper.classList.remove('valid-gene');
                    wrapper.classList.add('invalid-gene');
                }}
            }}

            function attach() {{
                const input = doc.querySelector('.st-key-tag_gene_id input');
                if (!input) {{ setTimeout(attach, 100); return; }}
                check(input);
                if (input.dataset.geneValidatorAttached) return;
                input.dataset.geneValidatorAttached = '1';
                input.addEventListener('input', () => check(input));
            }}
            attach();
        }})();
        </script>
        """,
        height=0,
    )

    tag_fasta_file = None
    if tag_choice == "Custom":
        if terminus == "N":
            tag_fasta_file = st.file_uploader(
                "Tag + Linker FASTA (required — must be [ATG + protein + linker] with no stop codon)",
                type=["fasta", "fa", "fna"], key="tag_fasta",
            )
        else:
            tag_fasta_file = st.file_uploader(
                "Linker + Tag FASTA (required — must be [linker + protein + stop])",
                type=["fasta", "fa", "fna"], key="tag_fasta",
            )

    col1, col2, col3, _spacer_bottom = st.columns([1, 1, 1, 4])
    with col1:
        overlap_length = st.number_input(
            "Vector overlap (bp)",
            value=20,
            min_value=0,
            key="tag_overlap",
            help="Length of the vector overlap tail prepended to AB_fwd and CD_rev (for HiFi assembly into the digested vector).",
        )
    with col2:
        junction_overlap = st.number_input(
            "Junction overlap (bp)",
            value=10,
            min_value=0,
            key="tag_junction",
            help="Length contributed by each side to the AB↔LT and LT↔CD junction overlaps (total junction width is 2× this value).",
        )
    with col3:
        flank_length = st.number_input(
            "Flank length (bp)",
            value=650,
            min_value=1,
            step=50,
            key="tag_flank",
            help="Length of the upstream and downstream genomic flanks amplified alongside the gene for HiFi assembly.",
        )

    with st.expander("Advanced parameters"):
        col1, col2, col3 = st.columns(3)
        with col1:
            opt_primer_size = st.number_input("Optimal primer size", value=20, min_value=10, key="tag_opt_size")
        with col2:
            min_primer_size = st.number_input("Min primer size", value=18, min_value=10, key="tag_min_size")
        with col3:
            max_primer_size = st.number_input("Max primer size", value=28, min_value=10, key="tag_max_size")

        col1, col2, col3 = st.columns(3)
        with col1:
            opt_tm = st.number_input("Optimal Tm", value=60.0, key="tag_opt_tm")
        with col2:
            mv_conc = st.number_input("Monovalent cation conc (mM)", value=500.0, key="tag_mv_conc")
        with col3:
            gc_clamp = st.number_input("GC clamp", value=1, min_value=0, key="tag_gc_clamp")

        linear_plasmid = st.checkbox("Linear plasmid", value=False, key="tag_linear")
        allow_unmatched = st.checkbox("Allow unmatched genes", value=False, key="tag_allow_unmatched")

    if st.button("Design primers", type="primary", use_container_width=True, key="tag_run"):
        plasmid_path = resolve_plasmid(plasmid_name, plasmid_upload)
        genome_paths = resolve_genome(genome_files)
        genes_path = resolve_genes(genes_file)
        if not plasmid_path:
            st.error("Select or upload a plasmid FASTA.")
        elif not genome_paths:
            st.error("No genome files found. Upload genome FASTA(s).")
        elif not genes_path:
            st.error("No genes file found. Upload a genes FASTA.")
        elif not five_prime_enzyme or not three_prime_enzyme:
            st.error("Both enzymes are required.")
        elif not gene_id.strip():
            st.error("Gene ID is required for protein tagging.")
        elif tag_choice == "Custom" and not tag_fasta_file:
            st.error("Upload a tag FASTA file above or select a hardcoded tag.")
        else:
            with st.spinner("Designing primers..."):
                fd, out_path = tempfile.mkstemp(suffix=".csv")
                os.close(fd)

                args = [
                    "--plasmid", plasmid_path,
                    "--genome", *genome_paths,
                    "--genes", genes_path,
                    "--output", out_path,
                    "--three-prime-enzyme", _resolve_enzyme(three_prime_enzyme),
                    "--five-prime-enzyme", _resolve_enzyme(five_prime_enzyme),
                    "--terminus", terminus,
                    "--overlap-length", str(overlap_length),
                    "--junction-overlap", str(junction_overlap),
                    "--flank-length", str(flank_length),
                    "--opt-primer-size", str(opt_primer_size),
                    "--min-primer-size", str(min_primer_size),
                    "--max-primer-size", str(max_primer_size),
                    "--opt-tm", str(opt_tm),
                    "--mv-conc", str(mv_conc),
                    "--gc-clamp", str(gc_clamp),
                    "--gene-ids", gene_id.strip(),
                ]
                if tag_choice == "Custom":
                    tag_path = save_uploaded(tag_fasta_file)
                    args.extend(["--tag", tag_path])
                elif tag_choice != "Custom":
                    args.extend(["--tag", tag_choice])

                if linear_plasmid:
                    args.append("--linear-plasmid")
                if allow_unmatched:
                    args.append("--allow-unmatched-genes")

                rc, _, stderr = run_script("design_protein_tag_primers.py", args, SCRIPT_DIR)
                if rc != 0 and not Path(out_path).exists():
                    st.error(f"Script exited with code {rc}")
                    show_stderr(stderr)
                else:
                    if tag_choice == "Custom" and tag_fasta_file is not None:
                        tag_label = Path(tag_fasta_file.name).stem
                    else:
                        tag_label = tag_choice
                    download_name = (
                        f"tag_{_safe_part(gene_id.strip())}_{terminus}term"
                        f"_{_safe_part(tag_label)}_{_date_stamp()}.csv"
                    )
                    show_results(out_path, stderr, download_filename=download_name)
