#!/usr/bin/env python3
"""Aggregate RVMT TSVs into a compact JSON for the index.html visualization."""
import csv, json, math, sys
from collections import defaultdict, Counter

csv.field_size_limit(sys.maxsize)

INFO  = "RiboV1.4_Info.tsv"
HMM   = "RiboV1.4_HMMatches.tsv"
IMG   = "RiboV1.4_IMG_Scaffold_cart.tsv"

def safe_int(x):
    try: return int(x)
    except: return 0

# ---- Info.tsv pass ----
phylum_ct  = Counter()
tree       = {}                      # phylum -> class -> order -> family -> count
host_kingdom_ct = Counter()
host_phylum_pairs = Counter()        # (kingdom, phylum)
euk_clade_ct = Counter()
novel_by_phy = defaultdict(Counter)
seg_by_phy   = defaultdict(Counter)
source_ct    = Counter()
length_by_phy = defaultdict(list)
top_families = Counter()
genus_by_fam = defaultdict(Counter)
genetic_code_ct = Counter()
afflvl_ct = Counter()

EUK_TIERS = {"Metazoa","Viridiplantae","Fungi","Stramenopiles","Alveolata",
             "Rhizaria","Amoebozoa","Apicomplexa","Discoba","Haptista","Cryptophyta"}

print("Reading Info.tsv...", file=sys.stderr)
with open(INFO, newline="") as f:
    rdr = csv.reader(f, delimiter="\t")
    header = next(rdr)
    for row in rdr:
        if len(row) < 24: continue
        seg     = row[3]
        host    = row[13] or "NA"
        length  = safe_int(row[14])
        phylum  = row[16] or "Unclassified"
        cls     = row[17] or "Unclassified"
        order   = row[18] or "Unclassified"
        family  = row[19] or "Unclassified"
        genus   = row[20] or "Unclassified"
        novel   = row[21]
        source  = row[22] or "Unknown"
        gcode   = row[6] or "NA"
        afl     = (row[11] or "NA").split(" - ")[0]

        phylum_ct[phylum] += 1
        p = tree.setdefault(phylum, {"_count":0,"children":{}})
        p["_count"] += 1
        c = p["children"].setdefault(cls, {"_count":0,"children":{}})
        c["_count"] += 1
        o = c["children"].setdefault(order, {"_count":0,"children":{}})
        o["_count"] += 1
        fa = o["children"].setdefault(family, {"_count":0})
        fa["_count"] += 1

        if host == "NA" or not host:
            kingdom = "Unknown"
        else:
            parts = [s.strip() for s in host.split(";")]
            kingdom = parts[0]
            if kingdom == "Eukaryota":
                tier = next((t for t in parts if t in EUK_TIERS), "Other Eukaryota")
                euk_clade_ct[tier] += 1
        host_kingdom_ct[kingdom] += 1
        host_phylum_pairs[(kingdom, phylum)] += 1

        novel_by_phy[phylum][novel] += 1
        seg_by_phy[phylum][seg] += 1
        source_ct[source] += 1
        if length > 0: length_by_phy[phylum].append(length)
        top_families[(phylum, family)] += 1
        genus_by_fam[family][genus] += 1
        genetic_code_ct[gcode] += 1
        afflvl_ct[afl] += 1

# ---- length histograms (log bins) ----
def loghist(vals, bins=24):
    vs = [math.log10(v) for v in vals if v > 0]
    if not vs: return {"edges":[],"counts":[]}
    lo, hi = min(vs), max(vs)
    if hi == lo: hi = lo + 1
    w = (hi-lo)/bins
    cnt = [0]*bins
    for v in vs:
        i = min(int((v-lo)/w), bins-1)
        cnt[i] += 1
    edges = [10**(lo + i*w) for i in range(bins+1)]
    return {"edges":edges, "counts":cnt}

length_hist = {p: loghist(v) for p,v in length_by_phy.items()}

# ---- HMM pass ----
print("Reading HMMatches.tsv (large)...", file=sys.stderr)
profile_ct = Counter()
profile_phy = defaultdict(Counter)
clan_ct = Counter()
analysis_ct = Counter()
hmm_total = 0
classified_ct = Counter()

with open(HMM, newline="") as f:
    rdr = csv.reader(f, delimiter="\t")
    header = next(rdr)
    for i, row in enumerate(rdr):
        if len(row) < 22: continue
        new_name = row[9] or "Unknown"
        analysis = row[10] or "NA"
        clan     = row[12] or "Unknown"
        classed  = row[14] or "NA"
        phylum   = row[17] or "Unclassified"
        profile_ct[new_name] += 1
        profile_phy[new_name][phylum] += 1
        clan_ct[clan] += 1
        analysis_ct[analysis] += 1
        classified_ct[classed] += 1
        hmm_total += 1
        if i % 500000 == 0 and i:
            print(f"  {i:,} rows", file=sys.stderr)

# IMG scaffold count
print("Reading IMG cart...", file=sys.stderr)
img_total = 0
with open(IMG) as f:
    next(f)
    for _ in f: img_total += 1

# ---- shape outputs ----
def tree_to_d3(name, node):
    d = {"name": name, "value": node["_count"]}
    if "children" in node and node["children"]:
        d["children"] = [tree_to_d3(k, v) for k,v in node["children"].items()]
    return d

sunburst = {"name":"RVMT", "children":[tree_to_d3(k,v) for k,v in tree.items()]}

# Sankey: host kingdom -> phylum
sankey_nodes = []
sankey_idx = {}
def add_node(label, group):
    if label not in sankey_idx:
        sankey_idx[label] = len(sankey_nodes)
        sankey_nodes.append({"name":label, "group":group})
    return sankey_idx[label]

sankey_links = []
for (kingdom, phylum), v in host_phylum_pairs.items():
    a = add_node(kingdom, "host")
    b = add_node(phylum, "phylum")
    sankey_links.append({"source":a, "target":b, "value":v})

# Top families per phylum
fam_per_phy = defaultdict(list)
for (phy, fam), v in top_families.items():
    fam_per_phy[phy].append({"family":fam, "count":v})
for phy in fam_per_phy:
    fam_per_phy[phy].sort(key=lambda x: -x["count"])
    fam_per_phy[phy] = fam_per_phy[phy][:14]

# Top profiles overall
top_profiles = []
for prof, ct in profile_ct.most_common(40):
    top_profiles.append({
        "name": prof, "count": ct,
        "phylum_dist": dict(profile_phy[prof].most_common(6))
    })

# Top clans
top_clans = [{"name":n,"count":c} for n,c in clan_ct.most_common(30)]

novelty = {p: dict(c) for p,c in novel_by_phy.items()}
segmentation = {p: dict(c) for p,c in seg_by_phy.items()}

out = {
    "totals": {
        "viruses": sum(phylum_ct.values()),
        "phyla": len([p for p in phylum_ct if not p.startswith("p.")]),
        "families": sum(1 for k in {(p,f) for (p,f) in top_families}),
        "genera": len({g for f in genus_by_fam for g in genus_by_fam[f]}),
        "hmm_hits": hmm_total,
        "img_scaffolds": img_total,
        "novel_pct": round(100 * sum(novel_by_phy[p]["TRUE"] for p in novel_by_phy) / max(1, sum(phylum_ct.values())), 1),
    },
    "phylum_count": dict(phylum_ct),
    "sunburst": sunburst,
    "host_kingdom": dict(host_kingdom_ct),
    "euk_clade": dict(euk_clade_ct),
    "sankey": {"nodes": sankey_nodes, "links": sankey_links},
    "novelty": novelty,
    "segmentation": segmentation,
    "source": dict(source_ct),
    "length_hist": length_hist,
    "top_families_per_phylum": dict(fam_per_phy),
    "top_profiles": top_profiles,
    "top_clans": top_clans,
    "analysis_type": dict(analysis_ct),
    "classified": dict(classified_ct),
    "genetic_code": dict(genetic_code_ct),
    "afflvl": dict(afflvl_ct),
}

with open("data.json", "w") as f:
    json.dump(out, f, separators=(",",":"))
print(f"Wrote data.json ({sum(1 for _ in open('data.json'))} lines)", file=sys.stderr)
