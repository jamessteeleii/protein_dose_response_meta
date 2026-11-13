"""Render, validate and package the preprint without running the analysis."""

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import html
import json
import os
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
from urllib.parse import unquote, urlsplit
import zipfile

from pypdf import PdfReader, PdfWriter
from pypdf.generic import NameObject, TextStringObject


PAGES_URL = "https://jamessteeleii.github.io/protein_dose_response_meta/"
REPOSITORY_URL = "https://github.com/jamessteeleii/protein_dose_response_meta"
REVIEW_PORT = 8766
REVIEW_URL = f"http://127.0.0.1:{REVIEW_PORT}/"
MANUSCRIPT = Path("manuscript/manuscript.pdf")
WORKBOOKS = (
    ("data/protein_meta_data_all.xlsx", "Download canonical extraction workbook (.xlsx)"),
    ("data/update_2026/protein_meta_data_update_2026.xlsx", "Download 2026 search-update workbook (.xlsx)"),
)
SUPPLEMENTS = (
    ("prisma-reporting", "prisma-reporting", "PRISMA reporting and study characteristics"),
    ("systematic-review-data-audit", "systematic-review-data-audit", "Systematic review data audit"),
    ("primary-dose-response", "primary-dose-response", "Primary dose–response figures"),
    ("relative-intake-sensitivity-analyses", "sensitivity-analyses", "Relative-intake sensitivities"),
    ("absolute-intake-sensitivity-analyses", "absolute-intake-sensitivity-analyses", "Absolute total intake"),
    ("absolute-supplementation-sensitivity-analyses", "absolute-supplementation-sensitivity-analyses", "Prescribed supplementation"),
    ("muscle-measurement-sensitivity-analyses", "muscle-measurement-sensitivity-analyses", "Muscle measurement"),
    ("age-duration-interactions", "age-duration-interactions", "Age and intervention duration"),
    ("training-status", "training-status", "Training status"),
    ("risk-of-bias", "risk-of-bias", "Risk of bias"),
    ("grade-appraisal", "grade-appraisal", "GRADE certainty"),
    ("publication-bias", "publication-bias", "Publication bias"),
    ("preregistration-deviations", "preregistration-deviations", "Preregistration deviations"),
    ("protein_calculations", "protein_calculations", "Protein calculations"),
    ("software-citations", "grateful-report", "Analysis software and citations"),
)


class DocumentLinks(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.anchors = set()
        self.links = set()
        self.feed(source)

    def handle_starttag(self, tag, attributes):
        attrs = dict(attributes)
        if attrs.get("id"):
            self.anchors.add(attrs["id"])
        if tag == "a" and attrs.get("name"):
            self.anchors.add(attrs["name"])
        key = {"a": "href", "link": "href", "img": "src", "script": "src"}.get(tag)
        if key and attrs.get(key) and not attrs[key].startswith("data:"):
            self.links.add(attrs[key])


def target_of(root, origin, link, index):
    if link.startswith((PAGES_URL, REVIEW_URL)):
        prefix = PAGES_URL if link.startswith(PAGES_URL) else REVIEW_URL
        origin = root / "index.html"
        link = link.removeprefix(prefix)
    url = urlsplit(link)
    if url.scheme == "file" or url.hostname in {"localhost", "127.0.0.1"}:
        raise ValueError(f"Unexpected machine-local link in {origin.name}: {link}")
    if url.scheme or url.netloc:
        return None
    path = unquote(url.path)
    target = (origin.parent / PurePosixPath(path)).resolve() if path else origin
    if not target.is_relative_to(root) or (target != index and not target.is_file()):
        raise ValueError(f"Missing or invalid link in {origin.name}: {link}")
    return target, unquote(url.fragment)


def pdf_links(document):
    for page in document.pages:
        for reference in page.get("/Annots", []):
            annotation = reference.get_object()
            action = annotation.get("/A")
            if action and action.get("/S") == "/URI":
                yield action


def collect_files(root, index_html):
    manuscript = root / MANUSCRIPT
    index = root / "index.html"
    supplements = {root / "supplement" / f"{output}.html" for _, output, _ in SUPPLEMENTS}
    documents = {path: DocumentLinks(path.read_text(encoding="utf-8")) for path in supplements}
    documents[index] = DocumentLinks(index_html)
    files = {manuscript, *supplements}
    workbooks = {root / path for path, _ in WORKBOOKS}

    def check(origin, link):
        match = target_of(root, origin, link, index)
        if match is None:
            return
        target, anchor = match
        if target.suffix == ".html":
            if target not in documents:
                raise ValueError(f"Unlisted HTML file linked by {origin.name}: {link}")
            if anchor and anchor not in documents[target].anchors:
                raise ValueError(f"Missing anchor in {target.name}: {anchor}")
        elif not target.is_relative_to(root / "manuscript" / "generated") and not (
            target.is_relative_to(root / "supplement" / "data")
            or target.is_relative_to(root / "supplement" / "figures")
            or target == manuscript
            or target in workbooks
        ):
            raise ValueError(f"Link outside publication inputs in {origin.name}: {link}")
        if target != index:
            files.add(target)

    for origin, document in documents.items():
        for link in sorted(document.links):
            check(origin, link)
    pdf = PdfReader(manuscript)
    links = set(str(action["/URI"]) for action in pdf_links(pdf))
    if PAGES_URL not in links:
        raise ValueError("The manuscript PDF is missing the supplementary landing page URL")
    supplement_links = {link for link in links if link.startswith(PAGES_URL)}
    for link in sorted(links):
        if link.startswith((PAGES_URL, "file:", REVIEW_URL)):
            check(manuscript, link)
    print(f"Validated {len(pdf.pages)} PDF pages, {len(supplement_links)} distinct "
          f"manuscript supplement URLs, and {len(files) - len(supplements) - 1} linked assets.")
    return files


def check_public_assets(files):
    private_fields = re.compile(
        r"\bcorresponding_author_email\b|"
        r"\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b|"
        r"\b[A-Z]:[/\\](?:Users|projects)[/\\]|file:///",
        re.IGNORECASE,
    )
    for path in sorted(files):
        if path.name.endswith(".csv.gz"):
            contents = gzip.decompress(path.read_bytes()).decode("utf-8-sig")
        elif path.suffix in {".csv", ".txt"}:
            contents = path.read_text(encoding="utf-8-sig")
        else:
            continue
        if private_fields.search(contents):
            raise ValueError(f"Review correspondence fields or local paths in {path.name}")


def render(root):
    for source in ["manuscript/manuscript.qmd", *(
        f"supplement/{source}.qmd" for source, _, _ in SUPPLEMENTS
    )]:
        command = ["quarto", "render", source]
        if source.startswith("supplement/"):
            command += ["--to", "html"]
        print("Rendering", source, flush=True)
        subprocess.run(command, cwd=root, check=True)


def landing_page(root):
    citation_file = root / "supplement" / "includes" / "preprint-citation.md"
    citation = re.search(
        r"^> (.+) <(https://doi\.org/[^<>\s]+)>$",
        citation_file.read_text(encoding="utf-8"),
        re.MULTILINE,
    )
    if citation is None:
        raise ValueError(f"Missing preprint citation and DOI in {citation_file}")
    reference, doi = map(html.escape, citation.groups())

    downloads = "\n".join(
        f'<li><a href="{path}" download>{title}</a></li>'
        for path, title in WORKBOOKS
    )
    rows = "\n".join(
        f'<li><a href="supplement/{output}.html">{html.escape(title)}</a></li>'
        for _, output, title in SUPPLEMENTS
    )
    return (
        '<!doctype html>\n<html lang="en"><meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<title>Protein dose–response meta-analysis</title>\n'
        '<style>body{max-width:44rem;margin:3rem auto;padding:0 1rem;'
        'font:1.05rem/1.6 system-ui,sans-serif;color:#18232b}'
        'a{color:#075b85}</style>\n'
        '<h1>Protein dose–response meta-analysis</h1>\n'
        '<p><a href="manuscript/manuscript.pdf">Manuscript PDF</a></p>\n'
        '<p><strong>Please cite as:</strong></p>\n'
        f'<blockquote>{reference} <a href="{doi}">{doi}</a></blockquote>\n'
        f'<p>The <a href="{REPOSITORY_URL}">GitHub repository</a> contains the analysis code, '
        'data, manuscript sources, and frozen reporting files. See the '
        f'<a href="{REPOSITORY_URL}/blob/master/REPRODUCIBILITY.md">reproducibility guide</a> '
        'for instructions to rebuild the documents or rerun the analyses.</p>\n'
        '<h2>Extraction data</h2>\n'
        '<p>Open either workbook in a spreadsheet app. Use both for the updated review.</p>\n'
        f'<ul>\n{downloads}\n</ul>\n'
        f'<h2>Supplementary materials</h2>\n<ol>\n{rows}\n</ol>\n</html>\n'
    )


def rewrite_review_pdf(path):
    document = PdfWriter(clone_from=path)
    count = 0
    for action in pdf_links(document):
        uri = str(action["/URI"])
        if uri.startswith(PAGES_URL):
            action[NameObject("/URI")] = TextStringObject(uri.replace(PAGES_URL, REVIEW_URL, 1))
            count += 1
    temporary = path.with_suffix(".review.pdf")
    document.write(temporary)
    temporary.replace(path)
    links = [str(action["/URI"]) for action in pdf_links(PdfReader(path))]
    if sum(link.startswith(REVIEW_URL) for link in links) != count:
        raise ValueError("Could not verify review PDF links")
    return count


def write_manifest(destination, commit, tag, source_root):
    entries = []
    for path in sorted(destination.rglob("*")):
        if path.is_file():
            entries.append({
                "file": path.relative_to(destination).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            })
    tools = subprocess.check_output(["quarto", "--version"], text=True).strip()
    checksums = {}
    for relative in ("renv.lock", "config/data-manifest.yml"):
        path = source_root / relative
        if path.is_file():
            checksums[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    r_version = subprocess.check_output(
        [str(Path(os.environ["QUARTO_R"]) / "Rscript") if "QUARTO_R" in os.environ else "Rscript",
         "--vanilla", "-e", "cat(R.version.string)"], text=True
    ).strip()
    manifest = {"git_commit": commit, "tag": tag,
                "build_date_utc": datetime.now(timezone.utc).isoformat(),
                "quarto": tools, "r": r_version,
                "source_checksums": checksums, "files": entries}
    (destination / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def package(root, destination, review, commit, tag):
    index_html = landing_page(root)
    files = collect_files(root, index_html)
    check_public_assets(files)
    destination.mkdir(parents=True, exist_ok=False)
    for path in sorted(files):
        target = destination / path.relative_to(root)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    (destination / "index.html").write_text(index_html, encoding="utf-8")
    if review:
        count = rewrite_review_pdf(destination / MANUSCRIPT)
        for path in destination.glob("supplement/*.html"):
            original = path.read_text(encoding="utf-8")
            if PAGES_URL in original:
                path.write_text(original.replace(PAGES_URL, REVIEW_URL), encoding="utf-8")
        (destination / "open-review.py").write_text(
            'from functools import partial\n'
            'from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer\n'
            'from pathlib import Path\n'
            'import webbrowser\n\n'
            'class ReviewServer(ThreadingHTTPServer):\n'
            '    allow_reuse_address = False\n\n'
            'handler = partial(SimpleHTTPRequestHandler, directory=str(Path(__file__).resolve().parent))\n'
            f'with ReviewServer(("127.0.0.1", {REVIEW_PORT}), handler) as server:\n'
            f'    print("Serving the review copy on {REVIEW_URL}. Close this window when finished.")\n'
            f'    webbrowser.open("{REVIEW_URL}")\n'
            '    server.serve_forever()\n',
            encoding="ascii",
        )
        (destination / "Open review.cmd").write_text(
            '@echo off\ncd /d "%~dp0"\n'
            'python open-review.py\n'
            'pause\n',
            encoding="ascii",
        )
        print(f"Verified {count} PDF link annotations for local review.")
    write_manifest(destination, commit, tag, root)
    archive = destination.with_suffix(".zip")
    with zipfile.ZipFile(archive, mode="x", compression=zipfile.ZIP_DEFLATED) as bundle:
        for path in sorted(destination.rglob("*")):
            if path.is_file():
                bundle.write(path, path.relative_to(destination).as_posix())
    print(f"Packaged {destination} and {archive}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--render", action="store_true")
    parser.add_argument("--review", action="store_true")
    parser.add_argument("--commit", default=os.environ.get("GITHUB_SHA"))
    parser.add_argument("--tag", help="Stage tag for the release manifest")
    args = parser.parse_args()
    root = args.root.resolve()
    if args.render:
        render(root)
    package(root, args.destination.resolve(), args.review, args.commit, args.tag)


if __name__ == "__main__":
    main()
