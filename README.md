# RNA Virome Atlas

A visualization of the RVMT v1.4 catalogue: 378,253 RNA virus sequences
assembled by Neri and colleagues from 5,150 globally distributed
metatranscriptomes, roughly a five-fold expansion of the known RNA virosphere.

**Live site:** <https://curtishoffmann.com/rnaviromeatlas/>

Taxonomy, hosts, novelty, genome architecture, protein content and provenance,
with every chart backed by a table of the same figures.

## How it works

```
aggregate.py          Reads three RVMT tables and writes data/data.json
build.py              Renders index.html, sitemap.xml and llms.txt from it
src/template.html     The page, with a placeholder for every computed number
data/data.json        Counts and summary statistics; the only data committed
assets/styles.css     Stylesheet
assets/app.js         The charts, drawn with d3
index.html            Build output; do not edit by hand
```

Standard library only on the Python side. The page loads one stylesheet, d3
and d3-sankey pinned by version with subresource integrity, and one script.

### Build

```bash
python3 build.py            # render the page from data/data.json
python3 build.py --check    # print the figures, write nothing
```

Every number on the page is read from `data/data.json` in `build.py`: the hero
figures, the per-phylum counts in the prose, the aggregate tables under each
chart, the structured data and `llms.txt`. Nothing about the catalogue is typed
into the template, so re-running the aggregator updates the prose and the
charts together.

### Re-aggregating from the RVMT release

The three source tables total about 620 MB, so they are not committed. Download
the RiboV1.4 release, then:

```bash
export RVMT_DIR=/path/to/RiboV1.4
python3 aggregate.py
python3 build.py
```

`aggregate.py` needs `RiboV1.4_Info.tsv`, `RiboV1.4_HMMatches.tsv` and
`RiboV1.4_IMG_Scaffold_cart.tsv`, and stops with the expected path if one is
missing. It records the size and MD5 of each input, plus the date, into the
`provenance` block of `data/data.json`, so the committed figures can always be
traced to the files that produced them.

## Two corrections worth knowing about

**The source column.** An earlier version of `aggregate.py` read column 23
("Source"), which holds a per-sample accession and has 12,813 distinct values.
The provenance chart on this page is about the study of record, which is column
24 ("Set") and has six values. The committed `data.json` carried the output of
the wrong column while the published page carried the right one, so the two
disagreed. `aggregate.py` now reads column 24 and the comment in the file says
why.

**"Novel genera."** The page used to label 4,569 as the count of novel genera.
That figure is the number of genera in the catalogue, all of them named; the
release carries no "Unclassified" genus rows at all. Novelty is recorded per
sequence, not per genus: 348,762 sequences of 378,253 carry no match to a
reference, which is the 92.2% the page reports. The number of genera in which
no member has a reference match is 267, and that is the figure this site now
reports as novel genera.

## The one figure that is cited rather than computed

The count of 5,150 metatranscriptomes comes from the paper. The released tables
do not record it, and the nearest thing in them is 12,813 distinct sample
accessions, which is a different quantity. The page says so where it uses the
number.

## Preview locally

```bash
python3 -m http.server 8000
# open http://localhost:8000
```

## Licensing

Code is MIT; see [LICENSE](LICENSE). The catalogue is the work of Neri and
colleagues and its terms are set by the RVMT release. No RVMT table is
redistributed here. Anyone reusing the figures should cite the paper:

> Neri U, Wolf YI, Roux S, et al. Expansion of the global RNA virome reveals
> diverse clades of bacteriophages. *Cell*. 2022;185(21):4023-4037.
> [doi:10.1016/j.cell.2022.08.023](https://doi.org/10.1016/j.cell.2022.08.023)
