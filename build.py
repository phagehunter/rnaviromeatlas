#!/usr/bin/env python3
"""Render index.html from src/template.html and data/data.json.

    python3 build.py            # render the page, sitemap and llms.txt
    python3 build.py --check    # report the figures, write nothing

data/data.json is produced by aggregate.py from the RVMT v1.4 release. Every
number printed on this page is read from it here, so the prose and the charts
can never disagree, and a re-run of the aggregator updates both at once.

The one figure not taken from the data is the count of metatranscriptomes. The
released tables do not carry it, so it is cited to the paper in the page copy
rather than presented as computed.
"""
import argparse
import datetime
import html
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent
DATA = ROOT / "data" / "data.json"
TEMPLATE = ROOT / "src" / "template.html"

SITE_URL = "https://curtishoffmann.com/rnaviromeatlas/"
AUTHOR = "Curtis Hoffmann"
AUTHOR_URL = "https://curtishoffmann.com"

PAPER = ("Neri U, Wolf YI, Roux S, et al. Expansion of the global RNA virome reveals "
         "diverse clades of bacteriophages. Cell. 2022;185(21):4023-4037. "
         "doi:10.1016/j.cell.2022.08.023")
PAPER_DOI = "https://doi.org/10.1016/j.cell.2022.08.023"

# The paper's figure, not ours: the released tables do not record it.
METATRANSCRIPTOMES = "5,150"

PHY_COLOR = {
    "Lenarviricota": "#6cf2ff",
    "Pisuviricota": "#ff7aa2",
    "Kitrinoviricota": "#ffd58a",
    "Duplornaviricota": "#7dffb4",
    "Negarnaviricota": "#ff9f6e",
    "Unclassified": "#8d9bb0",
}
OTHER_COLOR = "#c599ff"


def e(s):
    return html.escape("" if s is None else str(s), quote=True)


def n(v):
    return f"{v:,}" if isinstance(v, int) else str(v)


def fill(template, mapping):
    def sub(m):
        key = m.group(1)
        if key not in mapping:
            raise KeyError(f"no value for placeholder {{{{{key}}}}}")
        return str(mapping[key])
    return re.sub(r"\{\{(\w+)\}\}", sub, template)


# ------------------------------------------------------------- text fallbacks

def table(summary, caption, headers, rows, colors=None):
    """One chart's figures as a table, present in the HTML without JavaScript."""
    cls = ' class="n"'
    head = "".join("<th" + (cls if i else "") + ">" + e(h) + "</th>"
                   for i, h in enumerate(headers))
    body = []
    for i, row in enumerate(rows):
        cells = []
        for j, cell in enumerate(row):
            if j == 0:
                swatch = ""
                if colors:
                    c = colors(row[0]) if callable(colors) else colors.get(row[0])
                    if c:
                        swatch = f'<span class="swatch" style="background:{e(c)}"></span>'
                cells.append(f"<td>{swatch}{e(cell)}</td>")
            else:
                cells.append(f'<td class="n">{e(n(cell)) if isinstance(cell, int) else e(cell)}</td>')
        body.append("<tr>" + "".join(cells) + "</tr>")
    return (f'<details class="data-table">\n'
            f"          <summary>{e(summary)}</summary>\n"
            f"          <table><caption>{e(caption)}</caption>\n"
            f"            <thead><tr>{head}</tr></thead>\n"
            f"            <tbody>{''.join(body)}</tbody>\n"
            f"          </table>\n        </details>")


def pct(part, whole):
    return f"{100 * part / whole:.1f}%" if whole else "—"


def build_tables(d):
    total = d["totals"]["viruses"]
    phy = {k: v for k, v in sorted(d["phylum_count"].items(), key=lambda kv: -kv[1])}
    color = lambda k: PHY_COLOR.get(k, OTHER_COLOR)

    t = {}

    t["table_phylum"] = table(
        "The same figures as a table",
        "Sequences per phylum in RVMT v1.4, and each phylum's share of the catalogue.",
        ["Phylum", "Sequences", "Share"],
        [[k, v, pct(v, total)] for k, v in phy.items()],
        colors=color)

    hosts = sorted(d["host_kingdom"].items(), key=lambda kv: -kv[1])
    t["table_hosts"] = table(
        "The same figures as a table",
        "Host of record for each sequence. Most sequences have no host assigned.",
        ["Host kingdom", "Sequences", "Share"],
        [[k, v, pct(v, total)] for k, v in hosts])

    euk = sorted(d["euk_clade"].items(), key=lambda kv: -kv[1])
    euk_total = sum(v for _, v in euk)
    t["table_euk"] = table(
        "The same figures as a table",
        "Eukaryotic host clade, for the sequences whose host of record is a eukaryote.",
        ["Clade", "Sequences", "Share of eukaryotic hosts"],
        [[k, v, pct(v, euk_total)] for k, v in euk])

    rows = []
    for k in phy:
        c = d["novelty"].get(k, {})
        nov, ref = c.get("TRUE", 0), c.get("FALSE", 0)
        if nov or ref:
            rows.append([k, nov, ref, pct(nov, nov + ref)])
    t["table_novelty"] = table(
        "The same figures as a table",
        "Sequences flagged novel against those matching a reference, per phylum.",
        ["Phylum", "Novel", "Reference", "Novel share"],
        rows, colors=color)

    rows = []
    for k in phy:
        c = d["segmentation"].get(k, {})
        seg = c.get("Segmented", 0) or c.get("TRUE", 0)
        mono = sum(v for kk, v in c.items() if kk not in ("Segmented", "TRUE"))
        if seg or mono:
            rows.append([k, seg, mono])
    t["table_segmentation"] = table(
        "The same figures as a table",
        "Segmented against monopartite genome architecture, per phylum.",
        ["Phylum", "Segmented", "Other"],
        rows, colors=color)

    rows = []
    for p, fams in sorted(d["top_families_per_phylum"].items(),
                          key=lambda kv: -phy.get(kv[0], 0)):
        if p.startswith("p.") or p == "Unclassified":
            continue
        for f in fams[:8]:
            rows.append([p, f["family"], f["count"]])
    t["table_families"] = table(
        "The same figures as a table",
        "The eight most populous families in each established phylum.",
        ["Phylum", "Family", "Sequences"],
        rows, colors=color)

    src = sorted(d["source"].items(), key=lambda kv: -kv[1])
    t["table_sources"] = table(
        "The same figures as a table",
        "Study of record for each sequence in the catalogue.",
        ["Source study", "Sequences", "Share"],
        [[k, v, pct(v, total)] for k, v in src])

    gc = sorted(d["genetic_code"].items(), key=lambda kv: -kv[1])
    t["table_gcode"] = table(
        "The same figures as a table",
        "Translation table inferred for each sequence.",
        ["Genetic code", "Sequences", "Share"],
        [["Unspecified" if k == "NA" else f"Code {k}", v, pct(v, total)] for k, v in gc])

    afl = sorted(d["afflvl"].items(), key=lambda kv: -kv[1])
    t["table_afl"] = table(
        "The same figures as a table",
        "Affiliation level, the confidence with which each sequence is placed in the taxonomy.",
        ["Affiliation level", "Sequences", "Share"],
        [[k, v, pct(v, total)] for k, v in afl])

    t["table_clans"] = table(
        "The same figures as a table",
        "The thirty most frequently matched protein clans across the catalogue.",
        ["Protein clan", "Profile hits"],
        [[c["name"], c["count"]] for c in d["top_clans"][:30]])

    return t


# ------------------------------------------------------------------- JSON-LD

def jsonld(d, description):
    tot = d["totals"]
    prov = d.get("provenance", {})
    return json.dumps({
        "@context": "https://schema.org",
        "@type": "Dataset",
        "name": "RNA Virome Atlas — aggregated counts from RVMT v1.4",
        "description": description,
        "url": SITE_URL,
        "creator": {"@type": "Person", "name": AUTHOR, "url": AUTHOR_URL},
        "isBasedOn": PAPER_DOI,
        "citation": PAPER,
        "dateModified": prov.get("aggregated", datetime.date.today().isoformat()),
        "isAccessibleForFree": True,
        "distribution": [{
            "@type": "DataDownload",
            "encodingFormat": "application/json",
            "contentUrl": SITE_URL + "data/data.json",
        }],
        "variableMeasured": [
            {"@type": "PropertyValue", "name": "Sequences", "value": tot["viruses"]},
            {"@type": "PropertyValue", "name": "Phyla", "value": tot["phyla"]},
            {"@type": "PropertyValue", "name": "Families", "value": tot["families"]},
            {"@type": "PropertyValue", "name": "Genera", "value": tot["genera"]},
            {"@type": "PropertyValue", "name": "Genera with no previously known member",
             "value": tot["novel_genera"]},
            {"@type": "PropertyValue", "name": "Protein-profile hits", "value": tot["hmm_hits"]},
        ],
    }, ensure_ascii=False, indent=2)


# --------------------------------------------------------------------- render

def render(d):
    tot = d["totals"]
    phy = d["phylum_count"]
    description = (
        f'Interactive atlas of the {n(tot["viruses"])} RNA virus sequences in the RVMT v1.4 '
        f'catalogue, assembled from {METATRANSCRIPTOMES} metatranscriptomes: taxonomy, hosts, '
        f'novelty, genome architecture and proteins.')
    mapping = {
        "title": f'RNA Virome Atlas — {n(tot["viruses"])} RNA viruses from RVMT v1.4',
        "og_title": f'RNA Virome Atlas — {n(tot["viruses"])} RNA viruses, visualized',
        "description": e(description),
        "og_image_alt": e(f'RNA Virome Atlas: {n(tot["viruses"])} RNA virus sequences from '
                          f'{METATRANSCRIPTOMES} metatranscriptomes'),
        "canonical": SITE_URL,
        "jsonld": jsonld(d, description),
        "v": (d.get("provenance", {}).get("aggregated") or
              datetime.date.today().isoformat()).replace("-", ""),
        "viruses": n(tot["viruses"]), "viruses_raw": str(tot["viruses"]),
        "genera": n(tot["genera"]), "genera_raw": str(tot["genera"]),
        "hmm": n(tot["hmm_hits"]), "hmm_raw": str(tot["hmm_hits"]),
        "novel_pct": str(tot["novel_pct"]), "novel_pct_raw": str(int(tot["novel_pct"])),
        "phy_lenar": n(phy["Lenarviricota"]),
        "phy_pisu": n(phy["Pisuviricota"]),
        "DATA": ('<script id="data" type="application/json">'
                 + json.dumps(d, separators=(",", ":")) + "</script>"),
        "SCRIPT": ('<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.min.js" '
                   'integrity="sha512-vc58qvvBdrDR4etbxMdlTt4GBQk1qjvyORR2nrsPsFPyrs+/u5c3+1Ct6upOgdZoIl7eq6k3a1UPDSNAQi/32A==" '
                   'crossorigin="anonymous" referrerpolicy="no-referrer"></script>\n'
                   '<script src="https://cdnjs.cloudflare.com/ajax/libs/d3-sankey/0.12.3/d3-sankey.min.js" '
                   'integrity="sha512-KK15oKpabNDaLpWinMtNfTqy/V7pzlc2FRG174PfASes7RRx6TAsua8HJdRTKo8+BLvPBKNIkL7kXWcz5HoqqA==" '
                   'crossorigin="anonymous" referrerpolicy="no-referrer"></script>\n'
                   f'<script src="assets/app.js?v={{v}}" defer></script>'),
    }
    mapping["SCRIPT"] = mapping["SCRIPT"].replace("{v}", mapping["v"])
    mapping.update(build_tables(d))
    return fill(TEMPLATE.read_text(encoding="utf-8"), mapping), description


def write_sitemap(d):
    day = d.get("provenance", {}).get("aggregated") or datetime.date.today().isoformat()
    (ROOT / "sitemap.xml").write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
        f"  <url><loc>{SITE_URL}</loc><lastmod>{day}</lastmod>"
        "<changefreq>yearly</changefreq><priority>1.0</priority></url>\n"
        "</urlset>\n", encoding="utf-8")


def write_llms(d):
    tot, prov = d["totals"], d.get("provenance", {})
    phy = "\n".join(f"- {k}: {n(v)} sequences"
                    for k, v in sorted(d["phylum_count"].items(), key=lambda kv: -kv[1]))
    src = "\n".join(f"- {k}: {n(v)} sequences"
                    for k, v in sorted(d["source"].items(), key=lambda kv: -kv[1]))
    (ROOT / "llms.txt").write_text(f"""# RNA Virome Atlas

> A visualization of the RVMT v1.4 catalogue: {n(tot["viruses"])} RNA virus sequences assembled by
> Neri and colleagues from {METATRANSCRIPTOMES} globally distributed metatranscriptomes, about a
> five-fold expansion of the known RNA virosphere. Taxonomy, hosts, novelty, genome architecture
> and protein content, with every chart backed by a table of the same figures.

Site: {SITE_URL}
Maintainer: {AUTHOR} ({AUTHOR_URL})
Source repository: https://github.com/phagehunter/rnaviromeatlas
Aggregated: {prov.get("aggregated", "unknown")}

## The data

- Catalogue: RVMT, release RiboV1.4, February 2022.
- Paper: {PAPER}
- The aggregated counts this page draws are committed at data/data.json and produced by
  aggregate.py from three RVMT tables that are too large to commit. See README.md.

## Shape of the catalogue

- {n(tot["viruses"])} sequences, {tot["phyla"]} established phyla, {n(tot["families"])} families,
  {n(tot["genera"])} genera
- {n(tot["novel_sequences"])} sequences ({tot["novel_pct"]}%) carry no match to a reference sequence
- {n(tot["novel_genera"])} genera have no previously known member
- {n(tot["hmm_hits"])} protein-profile hits across the catalogue

### Sequences per phylum

{phy}

### Study of record

{src}

## Two figures readers get wrong

- The {n(tot["genera"])} genera are not all novel. Novelty is recorded per sequence, not per genus.
  {n(tot["novel_genera"])} genera have no previously known member; that is the "novel genera" number.
- The count of {METATRANSCRIPTOMES} metatranscriptomes comes from the paper. The released tables do
  not carry it, so this site cites it rather than computing it.
""", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="report and write nothing")
    args = ap.parse_args()

    if not DATA.exists():
        raise SystemExit(f"missing {DATA}; run aggregate.py first (see README.md)")
    d = json.loads(DATA.read_text(encoding="utf-8"))
    page, description = render(d)

    tot = d["totals"]
    print(f'RVMT v1.4 — aggregated {d.get("provenance", {}).get("aggregated", "date unknown")}')
    for k, v in tot.items():
        print(f"  {k:18} {n(v)}")
    print(f'  {"metatranscriptomes":18} {METATRANSCRIPTOMES} (cited to the paper, not computed)')

    if args.check:
        print(f"\nwould write index.html ({len(page)/1024:.1f} KB), sitemap.xml, llms.txt")
        return

    (ROOT / "index.html").write_text(page, encoding="utf-8")
    write_sitemap(d)
    write_llms(d)
    print(f"\n  wrote index.html   {len(page)/1024:7.1f} KB")
    for f in ("sitemap.xml", "llms.txt"):
        print(f"  wrote {f:12} {(ROOT / f).stat().st_size/1024:7.1f} KB")


if __name__ == "__main__":
    main()
