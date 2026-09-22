"""Generate v6.1 manuscript and supplementary DOCX from Markdown."""
import subprocess, os

PROJECT = r"/path/to/xelox_project"
REPORTS = os.path.join(PROJECT, "reports")

# Use pandoc to convert MD to DOCX
files = [
    ("manuscript_draft_v6.1_submission.md", "manuscript_draft_v6.1_submission.docx"),
    ("supplementary_materials_v6.1.md", "supplementary_materials_v6.1.docx"),
]

for md_file, docx_file in files:
    md_path = os.path.join(REPORTS, md_file)
    docx_path = os.path.join(REPORTS, docx_file)
    
    if not os.path.exists(md_path):
        print(f"SKIP: {md_file} not found")
        continue
    
    cmd = [
        "pandoc", md_path,
        "-o", docx_path,
        "--reference-doc=reference.docx" if os.path.exists(os.path.join(REPORTS, "reference.docx")) else "",
        "--toc", "--toc-depth=3",
        "-f", "markdown",
        "-t", "docx",
    ]
    cmd = [c for c in cmd if c]  # remove empty strings
    
    print(f"Converting {md_file} -> {docx_file}...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=REPORTS)
        if result.returncode == 0:
            size = os.path.getsize(docx_path) / 1024
            print(f"  OK: {docx_file} ({size:.0f} KB)")
        else:
            print(f"  ERROR: {result.stderr[:200]}")
    except FileNotFoundError:
        print("  pandoc not found, trying python-docx fallback...")
        # Fallback: just copy the md as a note
        print(f"  MD file ready at: {md_path}")
        print(f"  Convert manually with: pandoc {md_file} -o {docx_file}")

print("\nDone.")
