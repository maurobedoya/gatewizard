#!/usr/bin/env python3
"""Write the GateWizard manual as one PDF.

The page order is the ``nav`` list in ``mkdocs.yml``. Snippet includes
(``--8<--``) are expanded the same way the HTML site expands them. Microsoft
Edge or Chrome prints the HTML. The cover uses the GUI's black emblem and
wordmark. Versions come from this checkout and the GUI checkout. Contributors
come from ``contributors.yml`` in each repository listed in
``docs/manual/repositories.yml``.

Install the docs tools, then run this from the repository root::

    pip install -r requirements-docs.txt
    python scripts/build_docs_pdf.py

The PDF is ``docs/manual/GateWizard-manual.pdf``. The README links to that
file. Pass ``-o`` for another path. Set ``DOCS_PDF_BROWSER`` to the browser
executable if Edge and Chrome are not in the usual locations. Set
``GATEWIZARD_GUI_ROOT`` if the GUI checkout is not next to this repository.
"""

from __future__ import annotations

import argparse
import base64
import html
import json
import os
import re
import shutil
import struct
import subprocess
import tempfile
import urllib.request
import zlib
from datetime import date
from pathlib import Path

import markdown
import yaml
from markdown.extensions.toc import TocExtension
from pypdf import PdfReader, PdfWriter
from pypdf.generic import (
    ArrayObject,
    DecodedStreamObject,
    DictionaryObject,
    NameObject,
    TextStringObject,
)
from pygments.formatters import HtmlFormatter

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
DEFAULT_PDF = DOCS / "manual" / "GateWizard-manual.pdf"
REPOSITORIES = DOCS / "manual" / "repositories.yml"

ATTR_RE = re.compile(
    r"""(?P<attr>\b(?:src|href))=(?P<q>["'])(?P<url>.*?)(?P=q)""",
    re.IGNORECASE,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_PDF,
        help=f"PDF path (default: {DEFAULT_PDF})",
    )
    parser.add_argument(
        "--html-only",
        action="store_true",
        help="Write the print HTML and skip the browser",
    )
    args = parser.parse_args()
    pages = flatten_nav(load_nav())
    if not pages:
        raise SystemExit("mkdocs.yml has no documentation pages in nav")
    api_version = read_version()
    gui_version = read_gui_version()
    body, toc = render_pages(pages)
    contributors = render_contributors()
    sections = [("doc-contributors", "Contributors", 1), *section_index(pages)]
    html_path = ROOT / "site" / "documentation.print.html"
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_text = render_document(body, toc, api_version, gui_version, contributors)
    html_path.write_text(html_text, encoding="utf-8")
    print(f"Wrote {html_path}")
    if args.html_only:
        return
    browser = find_browser()
    pdf_path = args.output.resolve()
    draft = pdf_path.with_suffix(".draft.pdf")
    print_pdf(browser, html_path, draft)
    page_numbers = read_section_pages(draft, sections)
    html_path.write_text(fill_toc_pages(html_text, page_numbers), encoding="utf-8")
    print_pdf(browser, html_path, pdf_path)
    draft.unlink(missing_ok=True)
    print(f"Wrote {pdf_path} ({pdf_path.stat().st_size} bytes)")


def load_nav() -> list:
    config = yaml.safe_load((ROOT / "mkdocs.yml").read_text(encoding="utf-8"))
    nav = config.get("nav")
    if not isinstance(nav, list):
        raise SystemExit("mkdocs.yml has no nav list")
    return nav


def flatten_nav(nav: list, depth: int = 1) -> list[dict]:
    """Return nav entries in site order, including section groups."""
    entries: list[dict] = []
    for item in nav:
        if not isinstance(item, dict):
            continue
        for title, target in item.items():
            name = str(title)
            if isinstance(target, str):
                entries.append(
                    {"kind": "page", "depth": depth, "title": name, "path": DOCS / target}
                )
            elif isinstance(target, list):
                entries.append({"kind": "group", "depth": depth, "title": name, "path": None})
                entries.extend(flatten_nav(target, depth + 1))
    return entries


def read_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'(?m)^version\s*=\s*"([^"]+)"', text)
    return match.group(1) if match else ""


def load_repositories() -> list[tuple[str, Path]]:
    """API and GUI checkouts, plus any later repository listed for the manual."""
    if REPOSITORIES.is_file():
        payload = yaml.safe_load(REPOSITORIES.read_text(encoding="utf-8")) or {}
        found: list[tuple[str, Path]] = []
        for item in payload.get("repositories") or []:
            if not isinstance(item, dict):
                continue
            name = str(item.get("name") or "").strip()
            raw = str(item.get("path") or "").strip()
            if not name or not raw:
                continue
            found.append((name, (ROOT / raw).resolve()))
        if found:
            return found
    return [("API", ROOT), ("GUI", (ROOT.parent / "gatewizard-gui").resolve())]


def read_gui_version() -> str:
    """Version from the GUI checkout's package.json."""
    override = os.environ.get("GATEWIZARD_GUI_ROOT", "").strip()
    candidates = []
    if override:
        candidates.append(Path(override) / "package.json")
    for name, path in load_repositories():
        if name.lower() == "gui":
            candidates.append(path / "package.json")
    candidates.append(ROOT.parent / "gatewizard-gui" / "package.json")
    for path in candidates:
        if not path.is_file():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        version = str(payload.get("version") or "").strip()
        if version:
            return version
    print("GUI package.json was not found. The cover will omit the GUI version.")
    return ""


def read_contributors(repo: Path) -> list[str]:
    path = repo / "contributors.yml"
    if not path.is_file():
        print(f"No contributors.yml in {repo}")
        return []
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    names = payload.get("contributors") or []
    return [str(item).strip() for item in names if str(item).strip()]


def render_contributors() -> str:
    blocks = [
        '<div class="doc-page doc-break">',
        '<h1 id="doc-contributors">Contributors</h1>',
    ]
    for name, path in load_repositories():
        people = read_contributors(path)
        accounts = github_accounts(path)
        blocks.append(f"<h2>{html.escape(name)}</h2>")
        if people:
            items = "\n".join(f"<li>{html.escape(person)}</li>" for person in people)
            blocks.append(f"<ul>\n{items}\n</ul>")
        if accounts:
            blocks.append("<h3>GitHub</h3>")
            links = "\n".join(
                f'<li><a href="{html.escape(url)}">{html.escape(login)}</a></li>'
                for login, url in accounts
            )
            blocks.append(f"<ul>\n{links}\n</ul>")
        if not people and not accounts:
            blocks.append("<p>No contributors are listed for this repository yet.</p>")
    blocks.append("</div>")
    return "\n".join(blocks)


def github_accounts(repo: Path) -> list[tuple[str, str]]:
    """Usernames and profile links from the public GitHub contributor list."""
    slug = github_slug(repo)
    if not slug:
        print(f"No GitHub remote for {repo}. GitHub usernames were skipped.")
        return []
    url = f"https://api.github.com/repos/{slug}/contributors?per_page=100&anon=0"
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "gatewizard-manual",
            "Accept": "application/vnd.github+json",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            payload = json.load(response)
    except Exception as exc:
        print(f"GitHub contributors for {slug} were not loaded: {exc}")
        return []
    if not isinstance(payload, list):
        print(f"GitHub contributors for {slug} were not loaded.")
        return []
    accounts: list[tuple[str, str]] = []
    for item in payload:
        if not isinstance(item, dict):
            continue
        if str(item.get("type") or "") == "Bot":
            continue
        login = str(item.get("login") or "").strip()
        profile = str(item.get("html_url") or "").strip()
        if not login or login.endswith("[bot]"):
            continue
        if not profile:
            profile = f"https://github.com/{login}"
        accounts.append((login, profile))
    return accounts


def github_slug(repo: Path) -> str:
    if not (repo / ".git").exists():
        return ""
    completed = subprocess.run(
        ["git", "-C", str(repo), "remote", "get-url", "origin"],
        capture_output=True,
        text=True,
        check=False,
    )
    remote = (completed.stdout or "").strip()
    match = re.search(r"github\.com[:/]([^/]+)/([^/]+?)(?:\.git)?$", remote)
    if match is None:
        return ""
    return f"{match.group(1)}/{match.group(2)}"


def brand_logo_uris() -> tuple[str, str]:
    """Black emblem and wordmark on transparent backgrounds, from the GUI checkout."""
    gui = next((path for name, path in load_repositories() if name.lower() == "gui"), None)
    if gui is None:
        gui = (ROOT.parent / "gatewizard-gui").resolve()
    logos = gui / "resources" / "brand" / "logos"
    emblem = logos / "app-window-light.png"
    wordmark = logos / "wordmark-light.png"
    missing = [path for path in (emblem, wordmark) if not path.is_file()]
    if missing:
        raise SystemExit(
            "Black GUI logos were not found:\n"
            + "\n".join(str(path) for path in missing)
        )
    return data_uri(emblem), data_uri(wordmark)


def data_uri(path: Path) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def make_markdown() -> markdown.Markdown:
    return markdown.Markdown(
        extensions=[
            "admonition",
            "attr_list",
            "footnotes",
            "md_in_html",
            "tables",
            "pymdownx.highlight",
            "pymdownx.inlinehilite",
            "pymdownx.snippets",
            "pymdownx.superfences",
            "pymdownx.tabbed",
            "pymdownx.details",
            "pymdownx.arithmatex",
            "pymdownx.critic",
            "pymdownx.caret",
            "pymdownx.keys",
            "pymdownx.mark",
            "pymdownx.tilde",
            TocExtension(permalink=False),
        ],
        extension_configs={
            "pymdownx.snippets": {"base_path": [str(ROOT)], "check_paths": True},
            "pymdownx.tabbed": {"alternate_style": True},
            "pymdownx.arithmatex": {"generic": True},
            "pymdownx.highlight": {"anchor_linenums": False},
        },
    )


def render_pages(pages: list[dict]) -> tuple[str, str]:
    converter = make_markdown()
    sections: list[str] = []
    toc: list[str] = ['<nav class="toc">', '<p class="toc-title">Contents</p>', "<ol>"]
    open_depth = 1
    for entry in pages:
        depth = int(entry["depth"])
        title = str(entry["title"])
        while open_depth < depth:
            toc.append("<ol>")
            open_depth += 1
        while open_depth > depth:
            toc.append("</ol></li>")
            open_depth -= 1
        if entry["kind"] == "group":
            slug = group_slug(title)
            fragment = ""
            css = "doc-group"
        else:
            path = entry["path"]
            if not path.is_file():
                raise SystemExit(f"Documentation page is missing: {path.relative_to(ROOT)}")
            slug = page_slug(path)
            fragment = converter.convert(source_without_white_splash(path.read_text(encoding="utf-8")))
            converter.reset()
            fragment = strip_matching_heading(fragment, title)
            fragment = demote_headings(fragment, depth)
            fragment = scope_ids(fragment, slug)
            fragment = rewrite_links(fragment, path, slug)
            css = "doc-page"
        level = min(depth, 6)
        break_class = " doc-break" if depth == 1 else ""
        sections.append(
            f'<div class="{css}{break_class}">\n'
            f'<h{level} id="{slug}">{html.escape(title)}</h{level}>\n'
            f"{fragment}\n</div>"
        )
        toc.append(
            f'<li><a class="toc-entry" href="#{slug}">'
            f'<span class="toc-label">{html.escape(title)}</span>'
            f'<span class="toc-leader"></span>'
            f'<span class="toc-page" data-toc="{slug}">—</span>'
            f"</a>"
        )
        print(f"  {title}")
    while open_depth > 1:
        toc.append("</ol></li>")
        open_depth -= 1
    toc.append("</ol></nav>")
    return "\n".join(sections), "\n".join(toc)


def source_without_white_splash(source: str) -> str:
    """Use the black splash on the home page. The white one disappears on paper."""
    splash = data_uri(light_splash_path())
    hero = (
        '<div class="doc-hero">'
        f'<img class="doc-hero-splash" src="{splash}" alt="GateWizard" />'
        "</div>"
    )
    source = re.sub(r'<div class="gw-hero".*?</div>', hero, source, flags=re.DOTALL)
    return re.sub(
        r"!\[[^\]]*]\(([^)]*splash\.png)\)",
        hero,
        source,
        flags=re.IGNORECASE,
    )


def light_splash_path() -> Path:
    """Black splash for a white page, kept next to the GUI logos."""
    gui = next((path for name, path in load_repositories() if name.lower() == "gui"), None)
    if gui is None:
        gui = (ROOT.parent / "gatewizard-gui").resolve()
    logos = gui / "resources" / "brand" / "logos"
    destination = logos / "splash-light.png"
    source = logos / "splash.png"
    if not source.is_file():
        source = DOCS / "images" / "splash.png"
    if not source.is_file():
        raise SystemExit("splash.png was not found in the GUI logos or docs/images")
    if destination.is_file() and destination.stat().st_mtime >= source.stat().st_mtime:
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    _write_black_splash(source, destination)
    return destination


def _write_black_splash(source: Path, destination: Path) -> None:
    """Turn the white splash artwork black and keep the transparent background."""
    width, height, rows = _read_rgba(source)
    black_rows = []
    for row in rows:
        out = bytearray(len(row))
        for index in range(0, len(row), 4):
            out[index + 3] = row[index + 3]
        black_rows.append(out)
    _write_rgba(destination, width, height, black_rows)


def _read_rgba(path: Path) -> tuple[int, int, list[bytearray]]:
    data = path.read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise SystemExit(f"Not a PNG: {path}")
    pos = 8
    width = height = color = 0
    idat = b""
    while pos < len(data):
        length = struct.unpack(">I", data[pos : pos + 4])[0]
        kind = data[pos + 4 : pos + 8]
        chunk = data[pos + 8 : pos + 8 + length]
        if kind == b"IHDR":
            width, height, _bit, color = struct.unpack(">IIBB", chunk[:10])
        elif kind == b"IDAT":
            idat += chunk
        elif kind == b"IEND":
            break
        pos += 12 + length
    if color != 6:
        raise SystemExit(f"Splash must be an RGBA PNG: {path}")
    raw = zlib.decompress(idat)
    channels = 4
    stride = width * channels
    rows: list[bytearray] = []
    cursor = 0
    previous = bytearray(stride)
    for _y in range(height):
        filter_type = raw[cursor]
        cursor += 1
        row = bytearray(raw[cursor : cursor + stride])
        cursor += stride
        if filter_type == 1:
            for index in range(stride):
                left = row[index - channels] if index >= channels else 0
                row[index] = (row[index] + left) & 255
        elif filter_type == 2:
            for index in range(stride):
                row[index] = (row[index] + previous[index]) & 255
        elif filter_type == 3:
            for index in range(stride):
                left = row[index - channels] if index >= channels else 0
                row[index] = (row[index] + ((left + previous[index]) // 2)) & 255
        elif filter_type == 4:
            for index in range(stride):
                left = row[index - channels] if index >= channels else 0
                up = previous[index]
                upper_left = previous[index - channels] if index >= channels else 0
                estimate = left + up - upper_left
                nearest = (
                    left
                    if abs(estimate - left) <= abs(estimate - up)
                    and abs(estimate - left) <= abs(estimate - upper_left)
                    else up
                    if abs(estimate - up) <= abs(estimate - upper_left)
                    else upper_left
                )
                row[index] = (row[index] + nearest) & 255
        elif filter_type != 0:
            raise SystemExit(f"Unsupported PNG filter {filter_type} in {path}")
        rows.append(row)
        previous = row
    return width, height, rows


def _write_rgba(path: Path, width: int, height: int, rows: list[bytearray]) -> None:
    def chunk(tag: bytes, payload: bytes) -> bytes:
        check = zlib.crc32(tag + payload) & 0xFFFFFFFF
        return struct.pack(">I", len(payload)) + tag + payload + struct.pack(">I", check)

    raw = b"".join(b"\x00" + bytes(row) for row in rows)
    header = struct.pack(">IIBBBBB", width, height, 8, 6, 0, 0, 0)
    path.write_bytes(
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", header)
        + chunk(b"IDAT", zlib.compress(raw, 9))
        + chunk(b"IEND", b"")
    )


def section_index(pages: list[dict]) -> list[tuple[str, str, int]]:
    rows: list[tuple[str, str, int]] = []
    for entry in pages:
        title = str(entry["title"])
        depth = int(entry["depth"])
        if entry["kind"] == "group":
            rows.append((group_slug(title), title, depth))
        else:
            rows.append((page_slug(entry["path"]), title, depth))
    return rows


def fill_toc_pages(html_text: str, pages: dict[str, int]) -> str:
    def replace(match: re.Match[str]) -> str:
        slug = match.group(1)
        number = pages.get(slug)
        label = str(number) if number else "—"
        return f'<span class="toc-page" data-toc="{slug}">{label}</span>'

    return re.sub(
        r'<span class="toc-page" data-toc="([^"]*)">.*?</span>',
        replace,
        html_text,
    )


def read_section_pages(pdf_path: Path, sections: list[tuple[str, str, int]]) -> dict[str, int]:
    """Page number of each contents entry, from the PDF bookmarks.

    A page can contain a heading with the same words as a later chapter
    ("Installation" inside Home). Only a bookmark at that chapter's depth counts.
    """
    reader = PdfReader(str(pdf_path))
    found: dict[str, int] = {}
    index = 0

    def walk(items, depth: int = 0) -> None:
        nonlocal index
        for item in items or []:
            if isinstance(item, list):
                walk(item, depth + 1)
                continue
            if index >= len(sections):
                return
            title = undouble_title(str(item.title))
            slug, expected, section_depth = sections[index]
            if title != expected or depth != section_depth - 1:
                continue
            found[slug] = reader.get_destination_page_number(item) + 1
            index += 1

    walk(reader.outline)
    missing = [title for _slug, title, _depth in sections if _slug not in found]
    if missing:
        print("Contents entries without a page number: " + ", ".join(missing))
    return found


def undouble_title(title: str) -> str:
    size = len(title)
    if size >= 2 and size % 2 == 0 and title[: size // 2] == title[size // 2 :]:
        return title[: size // 2]
    return title


def group_slug(title: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
    return f"doc-group-{slug}"


def strip_matching_heading(fragment: str, title: str) -> str:
    """Drop the page's first heading when it repeats the navigation title."""
    match = re.search(r"<h([1-6])([^>]*)>(.*?)</h\1>", fragment, flags=re.IGNORECASE | re.DOTALL)
    if match is None:
        return fragment
    text = re.sub(r"<[^>]+>", "", match.group(3))
    text = html.unescape(text).strip()
    if text.lower() != title.strip().lower():
        return fragment
    return fragment[: match.start()] + fragment[match.end() :]


def demote_headings(fragment: str, depth: int) -> str:
    """Place a page's headings under the navigation heading of that page."""
    ranks = [int(rank) for rank in re.findall(r"<h([1-6])\b", fragment, flags=re.IGNORECASE)]
    if not ranks:
        return fragment
    shift = (depth + 1) - min(ranks)
    if shift <= 0:
        return fragment

    def replace(match: re.Match[str]) -> str:
        rank = min(6, int(match.group(2)) + shift)
        return f"{match.group(1)}{rank}"

    return re.sub(r"(</?h)([1-6])\b", replace, fragment, flags=re.IGNORECASE)


def page_slug(path: Path) -> str:
    relative = path.resolve().relative_to(DOCS.resolve()).with_suffix("")
    return "doc-" + relative.as_posix().replace("/", "-")


def scope_ids(fragment: str, slug: str) -> str:
    """Keep heading anchors unique across pages."""
    fragment = re.sub(
        r'\bid="([^"]*)"',
        lambda match: f'id="{slug}-{match.group(1)}"',
        fragment,
    )
    return re.sub(
        r'href="#([^"]*)"',
        lambda match: f'href="#{slug}-{match.group(1)}"',
        fragment,
    )


def rewrite_links(fragment: str, source: Path, slug: str) -> str:
    def replace(match: re.Match[str]) -> str:
        attr = match.group("attr")
        quote = match.group("q")
        url = match.group("url")
        if url.startswith(("#", "http://", "https://", "mailto:", "data:")):
            return match.group(0)
        path_part, _sep, fragment_id = url.partition("#")
        if not path_part:
            return match.group(0)
        resolved = (source.parent / path_part).resolve()
        try:
            relative = resolved.relative_to(DOCS.resolve())
        except ValueError:
            relative = None
        if relative is not None and relative.suffix == ".md":
            target = page_slug(DOCS / relative)
            anchor = f"#{target}-{fragment_id}" if fragment_id else f"#{target}"
            return f"{attr}={quote}{anchor}{quote}"
        if attr.lower() == "src":
            return f"{attr}={quote}{embed_local_file(resolved)}{quote}"
        return f"{attr}={quote}{file_uri(resolved)}{quote}"

    return ATTR_RE.sub(replace, fragment)


def render_document(
    body: str,
    toc: str,
    api_version: str,
    gui_version: str,
    contributors: str,
) -> str:
    pygments_css = HtmlFormatter(style="default").get_style_defs(".highlight")
    emblem, wordmark = brand_logo_uris()
    api_line = (
        f"<p>API {html.escape(api_version)}</p>" if api_version else ""
    )
    gui_line = (
        f"<p>GUI {html.escape(gui_version)}</p>" if gui_version else ""
    )
    today = date.today().isoformat()
    contents = toc.replace(
        "<ol>",
        "<ol>\n"
        '<li><a class="toc-entry" href="#doc-contributors">'
        '<span class="toc-label">Contributors</span>'
        '<span class="toc-leader"></span>'
        '<span class="toc-page" data-toc="doc-contributors">—</span>'
        "</a></li>",
        1,
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>GateWizard Manual</title>
<style>
{PRINT_CSS}
{pygments_css}
</style>
</head>
<body>
<section class="cover">
  <img class="cover-emblem" src="{emblem}" alt="" />
  <img class="cover-wordmark" src="{wordmark}" alt="GateWizard" />
  <p class="kicker">Manual</p>
  <div class="versions">
    {api_line}
    {gui_line}
  </div>
  <ul class="cover-links">
    <li><a href="https://github.com/maurobedoya/gatewizard">API repository</a></li>
    <li><a href="https://github.com/franciscoadasme/gatewizard-gui">GUI repository</a></li>
    <li><a href="https://maurobedoya.github.io/gatewizard/">Documentation</a></li>
  </ul>
  <p class="generated">Generated {today}</p>
</section>
{contents}
{contributors}
{body}
</body>
</html>
"""


PRINT_CSS = """
@page {
  size: A4;
  margin: 16mm 15mm 18mm 15mm;
}
* { box-sizing: border-box; }
html { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
body {
  margin: 0;
  color: #1c1c1c;
  font-family: "Segoe UI", "Liberation Sans", sans-serif;
  font-size: 10.5pt;
  line-height: 1.45;
}
.cover { break-after: page; padding-top: 28mm; text-align: center; }
.cover-emblem { width: 38mm; height: auto; display: block; margin: 0 auto 7mm; }
.cover-wordmark { width: 92mm; height: auto; display: block; margin: 0 auto 12mm; }
.kicker {
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: #0d7377;
  font-size: 12pt;
  margin: 0 0 8mm;
}
.versions p { margin: 0.15em 0; font-size: 12pt; }
.cover-links {
  list-style: none;
  margin: 8mm 0 0;
  padding: 0;
  font-size: 11pt;
}
.cover-links li { margin: 0.2em 0; }
.cover-links a { color: #0b5f73; }
.generated { margin-top: 10mm; color: #555; font-size: 10pt; }
.toc { break-after: page; }
.toc-title { font-size: 18pt; font-weight: 650; margin: 0 0 0.6em; }
.doc-break { break-before: page; }
.toc ol { list-style: none; margin: 0; padding: 0; }
.toc ol ol { margin-left: 1.15em; }
.toc li { margin: 0.22em 0; }
.toc a, .toc-entry {
  display: flex;
  align-items: baseline;
  color: inherit;
  text-decoration: none;
}
.toc-label { flex: 0 1 auto; }
.toc-leader {
  flex: 1 1 auto;
  border-bottom: 1px dotted #888;
  margin: 0 0.45em 0.2em;
  min-width: 0.8em;
}
.toc-page { flex: 0 0 auto; font-variant-numeric: tabular-nums; }
.doc-hero { margin: 0 0 1em; text-align: center; }
.doc-hero-splash { width: 78mm; height: auto; display: block; margin-left: auto; margin-right: auto; }
h1, h2, h3, h4 { break-after: avoid; line-height: 1.25; }
h1 { font-size: 18pt; }
h2 { font-size: 14pt; margin-top: 1.4em; }
h3 { font-size: 12pt; }
p, li { orphans: 3; widows: 3; }
.doc-page p,
.doc-page li,
.doc-group p,
.doc-group li {
  text-align: justify;
  hyphens: auto;
}
.doc-page h1, .doc-page h2, .doc-page h3, .doc-page h4,
.doc-page pre, .doc-page code, .doc-page th, .doc-page td,
.doc-group h1, .doc-group h2, .doc-group h3, .doc-group h4 {
  text-align: left;
  hyphens: manual;
}
a { color: #0b5f73; }
img { max-width: 100%; height: auto; }
table { border-collapse: collapse; width: 100%; margin: 0.8em 0; font-size: 9pt; }
th, td { border: 1px solid #ccc; padding: 3px 6px; vertical-align: top; text-align: left; }
th { background: #f3f4f4; }
pre, code { font-family: "Cascadia Mono", "Consolas", monospace; font-size: 8.5pt; }
code { background: #f4f5f5; padding: 0 0.15em; }
pre {
  white-space: pre-wrap;
  word-break: break-word;
  background: #f6f7f7;
  border: 1px solid #e1e3e3;
  padding: 8px 10px;
  break-inside: avoid;
}
pre code { background: transparent; padding: 0; }
.admonition, details {
  border-left: 4px solid #0d7377;
  background: #f4fbfa;
  padding: 8px 12px;
  margin: 0.8em 0;
}
.admonition-title, summary { font-weight: 700; }
.tabbed-set > input { display: none; }
.tabbed-set > label {
  display: block;
  font-weight: 650;
  margin-top: 0.8em;
}
.tabbed-block { display: block !important; }
.gw-hero { margin: 0 0 1em; }
"""


IMAGE_TYPES = {
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".gif": "image/gif",
    ".webp": "image/webp",
    ".svg": "image/svg+xml",
}


def embed_local_file(path: Path) -> str:
    """Inline an image so the browser can print it without reading other files."""
    mime = IMAGE_TYPES.get(path.suffix.lower())
    if mime and path.is_file() and path.stat().st_size <= 8_000_000:
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        return f"data:{mime};base64,{encoded}"
    return file_uri(path)


def file_uri(path: Path) -> str:
    windows = to_windows(path)
    if re.match(r"^[A-Za-z]:[\\/]", windows):
        return "file:///" + windows.replace("\\", "/")
    return Path(windows).as_uri()


def to_windows(path: Path) -> str:
    """Map a WSL ``/mnt/c/...`` path to ``C:/...`` for Edge and Chrome."""
    text = path.resolve().as_posix()
    match = re.fullmatch(r"/mnt/([A-Za-z])/(.*)", text)
    if match:
        return f"{match.group(1).upper()}:/{match.group(2)}"
    return str(path.resolve())


def find_browser() -> Path:
    override = os.environ.get("DOCS_PDF_BROWSER", "").strip()
    if override:
        browser = Path(override)
        if browser.is_file():
            return browser
        raise SystemExit(f"DOCS_PDF_BROWSER is not a file: {override}")
    names = (
        "Microsoft/Edge/Application/msedge.exe",
        "Google/Chrome/Application/chrome.exe",
    )
    roots = (
        Path("/mnt/c/Program Files (x86)"),
        Path("/mnt/c/Program Files"),
        Path(r"C:\Program Files (x86)"),
        Path(r"C:\Program Files"),
    )
    for root in roots:
        for name in names:
            candidate = root / name
            if candidate.is_file():
                return candidate
    for command in ("msedge", "google-chrome", "chromium", "chromium-browser"):
        found = shutil.which(command)
        if found:
            return Path(found)
    raise SystemExit(
        "Edge or Chrome was not found. Set DOCS_PDF_BROWSER to the executable."
    )


def print_pdf(browser: Path, html_path: Path, pdf_path: Path) -> None:
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    if pdf_path.is_file():
        pdf_path.unlink()
    profile = Path(tempfile.mkdtemp(prefix="gw-docs-pdf-"))
    pdf_arg = to_windows(pdf_path)
    command = [
        str(browser),
        "--headless=new",
        "--disable-gpu",
        "--no-first-run",
        "--disable-extensions",
        f"--user-data-dir={to_windows(profile)}",
        "--no-pdf-header-footer",
        "--generate-pdf-document-outline",
        "--virtual-time-budget=20000",
        f"--print-to-pdf={pdf_arg}",
        file_uri(html_path),
    ]
    try:
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
    finally:
        shutil.rmtree(profile, ignore_errors=True)
    if not pdf_path.is_file() or pdf_path.stat().st_size < 1000:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise SystemExit(
            "The browser did not write a PDF"
            + (f"\n{detail}" if detail else "")
        )
    undouble_outline_titles(pdf_path)
    stamp_page_numbers(pdf_path)
    prefer_outline_bookmarks(pdf_path)


def undouble_outline_titles(pdf_path: Path) -> None:
    """Edge repeats some heading titles in the bookmark tree. Keep one copy."""
    reader = PdfReader(str(pdf_path))
    if not reader.outline:
        return
    writer = PdfWriter()
    writer.append(reader)

    def retitle(items) -> None:
        for item in items or []:
            if isinstance(item, list):
                retitle(item)
                continue
            title = undouble_title(str(item.title))
            if title != item.title:
                item[NameObject("/Title")] = TextStringObject(title)

    retitle(writer.outline)
    temporary = pdf_path.with_suffix(".tmp.pdf")
    with temporary.open("wb") as handle:
        writer.write(handle)
    temporary.replace(pdf_path)


def stamp_page_numbers(pdf_path: Path) -> None:
    """Write a page number in the bottom margin of every page.

    The number is a separate page merged on top. Drawing into the printed
    page keeps the browser's transformation, which places the number off the sheet.
    """
    reader = PdfReader(str(pdf_path))
    writer = PdfWriter()
    writer.append(reader)
    for index, page in enumerate(writer.pages, start=1):
        page.merge_page(_number_stamp(page, index))
    temporary = pdf_path.with_suffix(".tmp.pdf")
    with temporary.open("wb") as handle:
        writer.write(handle)
    temporary.replace(pdf_path)


def _number_stamp(page, number: int):
    box = page.cropbox or page.mediabox
    width = float(box.width)
    height = float(box.height)
    stamp_writer = PdfWriter()
    stamp = stamp_writer.add_blank_page(width=width, height=height)
    _draw_page_number(stamp_writer, stamp, number, width)
    return stamp


def prefer_outline_bookmarks(pdf_path: Path) -> None:
    """Keep the bookmark tree and drop the print tag tree.

    Chrome's viewer follows the bookmarks. Edge follows the tag tree from
    the HTML headings, which lists every heading inside a page and looks
    unordered. Without that tag tree, Edge uses the same bookmarks.
    """
    reader = PdfReader(str(pdf_path))
    writer = PdfWriter()
    writer.append(reader)
    root = writer.root_object
    for key in ("/StructTreeRoot", "/MarkInfo"):
        if key in root:
            del root[key]
    root[NameObject("/PageMode")] = NameObject("/UseOutlines")
    for page in writer.pages:
        if "/StructParents" in page:
            del page["/StructParents"]
    temporary = pdf_path.with_suffix(".tmp.pdf")
    with temporary.open("wb") as handle:
        writer.write(handle)
    temporary.replace(pdf_path)


def _draw_page_number(writer: PdfWriter, page, number: int, width: float) -> None:
    resources = page.get("/Resources")
    if resources is None:
        resources = DictionaryObject()
        page[NameObject("/Resources")] = resources
    else:
        resources = resources.get_object()
    fonts = resources.get("/Font")
    if fonts is None:
        fonts = DictionaryObject()
        resources[NameObject("/Font")] = fonts
    else:
        fonts = fonts.get_object()
    fonts[NameObject("/FPage")] = DictionaryObject(
        {
            NameObject("/Type"): NameObject("/Font"),
            NameObject("/Subtype"): NameObject("/Type1"),
            NameObject("/BaseFont"): NameObject("/Helvetica"),
        }
    )
    label = str(number)
    size = 10
    x = width / 2 - (size * 0.5 * len(label)) / 2
    stream = DecodedStreamObject()
    stream.set_data(
        f"q BT /FPage {size} Tf 0 0 0 rg {x:.2f} 28 Td ({label}) Tj ET Q".encode()
    )
    new_ref = writer._add_object(stream)
    existing = page.get("/Contents")
    if existing is None:
        page[NameObject("/Contents")] = new_ref
    elif isinstance(existing, ArrayObject):
        existing.append(new_ref)
    else:
        page[NameObject("/Contents")] = ArrayObject([existing, new_ref])


if __name__ == "__main__":
    main()
