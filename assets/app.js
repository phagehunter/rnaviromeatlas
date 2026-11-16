/* RNA Virome Atlas — charts, drawn with d3 from the aggregated catalogue.
 *
 * The data block is written into the page by build.py from data/data.json, so
 * the figures here and the figures in the prose come from the same run of
 * aggregate.py. Taxon, family and profile names come from the RVMT tables, so
 * every one of them is escaped before it reaches the DOM.
 *
 * Each chart runs inside guard(), which replaces a failed chart with a short
 * message rather than leaving an empty box, and every chart has a table of the
 * same figures underneath it for keyboard and screen-reader users.
 */

const D = JSON.parse(document.getElementById('data').textContent);
const tip = document.getElementById('tip');

/* Catalogue text is third-party data, so it is escaped on the way into markup. */
function esc(s){
  return String(s === null || s === undefined ? '' : s).replace(/[&<>"']/g, c => (
    {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
}

/* A legend row built from nodes, so a name can never be read as markup. */
function legendItem(color, label, trailing){
  const el = document.createElement('div');
  el.className = 'legend-item';
  const dot = document.createElement('div');
  dot.className = 'legend-dot';
  dot.style.background = color;
  el.appendChild(dot);
  el.appendChild(document.createTextNode(label));
  if (trailing){
    const t = document.createElement('span');
    t.style.opacity = '.6';
    t.style.marginLeft = '4px';
    t.textContent = trailing;
    el.appendChild(t);
  }
  return el;
}

/* Run one chart. A failure leaves a message in its container, not a blank. */
function guard(name, containerId, fn){
  try { fn(); }
  catch (err){
    const host = document.getElementById(containerId);
    if (host){
      host.textContent = '';
      const box = document.createElement('div');
      box.className = 'chart-error';
      const b = document.createElement('b');
      b.textContent = 'This chart could not be drawn.';
      box.appendChild(b);
      box.appendChild(document.createTextNode(
        ' The same figures are in the table below it, and in data/data.json in the repository.'));
      host.appendChild(box);
    }
    if (window.console) console.error('chart failed: ' + name, err);
  }
}

/* The tooltip is the one place markup is assembled; callers escape their text. */
function showTip(html, ev){tip.innerHTML=html;tip.style.opacity=1;moveTip(ev)}
function moveTip(ev){const x=ev.clientX+14,y=ev.clientY+14;tip.style.left=Math.min(x,innerWidth-260)+'px';tip.style.top=Math.min(y,innerHeight-140)+'px'}
function hideTip(){tip.style.opacity=0}
document.addEventListener('mousemove',ev=>{if(tip.style.opacity==='1')moveTip(ev)});

const PHY_COLOR = {
  'Lenarviricota':'#4cc9f0',
  'Pisuviricota':'#ff6b9d',
  'Kitrinoviricota':'#ffd166',
  'Duplornaviricota':'#80ffae',
  'Negarnaviricota':'#ff7a45',
  'Unclassified':'#6b7280',
  'p.0001.base-Kitrino':'#ada0ff',
  'p.0002.base-Kitrino':'#c8b8ff',
  'NA':'#3a3f55'
};
const phyColor = p => PHY_COLOR[p] || '#7a7f99';
const fmt = d3.format(',');
const fmtPct = d3.format('.1%');


/* A horizontal bar scale whose right edge is always past its left edge, so a
   narrow or not-yet-laid-out container cannot produce negative bar widths. */
function barScale(max, left, width, rightPad){
  const right = Math.max(left + 12, width - rightPad);
  return d3.scaleLinear().domain([0, max || 1]).range([left, right]);
}

/* Does this reader want motion at all? Used by the counters and the reveal. */
const wantsMotion = !window.matchMedia || !window.matchMedia('(prefers-reduced-motion: reduce)').matches;

/* Counters. build.py already wrote the final figure into the markup, so the
   animation only replaces text it is certain to finish on, and it is skipped
   entirely for readers who ask for reduced motion or whose browser lacks the
   observer. A half-finished count must never be left on screen. */
document.querySelectorAll('.stat .num').forEach(el => {
  const target = +el.dataset.count;
  const suffix = el.dataset.suffix || '';
  const settled = (target >= 1000 ? fmt(target) : String(target)) + suffix;
  if (!isFinite(target)){ return; }
  if (!wantsMotion || !('IntersectionObserver' in window)){
    el.textContent = settled;
    return;
  }
  let start = null;
  const dur = 2200;
  function tick(t){
    if (!start) start = t;
    const p = Math.min(1, (t - start) / dur);
    if (p >= 1){ el.textContent = settled; return; }
    const eased = 1 - Math.pow(1 - p, 3);
    const v = Math.floor(eased * target);
    el.textContent = (target >= 1000 ? fmt(v) : v) + suffix;
    requestAnimationFrame(tick);
  }
  const obs = new IntersectionObserver(es => es.forEach(e => {
    if (e.isIntersecting){ requestAnimationFrame(tick); obs.disconnect(); }
  }));
  obs.observe(el);
});

/* Reveal. The .js class on <html> is what switches the hidden starting state on,
   so a reader without JavaScript sees every section rather than a blank page.
   Readers who ask for reduced motion get the content without the movement. */
if (wantsMotion && 'IntersectionObserver' in window){
  document.documentElement.classList.add('js');
  const ro = new IntersectionObserver(
    es => es.forEach(e => { if (e.isIntersecting) e.target.classList.add('in'); }),
    {threshold: 0.08});
  document.querySelectorAll('.reveal').forEach(el => ro.observe(el));
  // Anything already on screen when the script runs should not wait for a scroll.
  requestAnimationFrame(() => document.querySelectorAll('.reveal').forEach(el => {
    const r = el.getBoundingClientRect();
    if (r.top < innerHeight && r.bottom > 0) el.classList.add('in');
  }));
}

/* ---- Sunburst ---- */
guard('sunburst', 'sunburst', function(){
  const W = 640, R = W/2;
  const root = d3.hierarchy(D.sunburst)
    .sum(d => d.children ? 0 : (d.value||0))
    .sort((a,b)=> b.value - a.value);
  d3.partition().size([2*Math.PI, root.height + 1])(root);
  root.each(d => d.current = d);

  const totalAll = root.value;

  const svg = d3.select('#sunburst').append('svg')
    .attr('viewBox', `${-R} ${-R} ${W} ${W}`)
    .attr('width', W).attr('height', W);

  const arc = d3.arc()
    .startAngle(d=>d.x0).endAngle(d=>d.x1)
    .padAngle(d=>Math.min((d.x1-d.x0)/2,0.0035))
    .padRadius(R*0.5)
    .innerRadius(d=>d.y0*R/(root.height+1))
    .outerRadius(d=>Math.max(d.y0*R/(root.height+1), d.y1*R/(root.height+1)-1));

  function colorFor(d){
    let n=d;while(n.depth>1)n=n.parent;
    const base = phyColor(n.data.name);
    const c = d3.color(base);
    if(!c) return base;
    const t = (d.depth-1)/3;
    return c.brighter(0).formatHex && t>0 ? d3.interpolate(base,'#0b0f2a')(t*0.45) : base;
  }

  const path = svg.append('g')
    .selectAll('path')
    .data(root.descendants().slice(1))
    .join('path')
      .attr('fill', d=>colorFor(d))
      .attr('fill-opacity', d=>arcVisible(d.current)? (d.children?0.85:0.65):0)
      .attr('pointer-events', d=>arcVisible(d.current)?'auto':'none')
      .attr('d', d=>arc(d.current));

  path.filter(d=>d.children).style('cursor','pointer').on('click', clicked);

  path.on('mouseover', (ev,d)=>{
    showTip(`<strong>${esc(d.data.name)}</strong><div class="small">${levelOf(d)} · ${fmt(d.value)} viruses (${fmtPct(d.value/totalAll)})</div>`,ev);
    updateInfo(d);
    path.attr('fill-opacity',o=>arcVisible(o.current)? (o===d?1:(o.children?0.55:0.35)):0);
  }).on('mouseout',()=>{
    hideTip();updateInfo(root);
    path.attr('fill-opacity',o=>arcVisible(o.current)? (o.children?0.85:0.65):0);
  });

  const center = svg.append('circle')
    .attr('r', R/(root.height+1))
    .attr('fill','rgba(11,15,42,0.95)')
    .attr('stroke','rgba(255,255,255,0.1)')
    .attr('pointer-events','all').style('cursor','pointer')
    .on('click', (ev)=>clicked(ev,root));

  const centerLabel = svg.append('text')
    .attr('text-anchor','middle').attr('dy','-0.1em')
    .style('font-family',"'Space Grotesk',sans-serif").style('fill','#e8eaf6')
    .style('font-size','13px').style('font-weight',500).style('pointer-events','none')
    .text('RVMT');
  const centerSub = svg.append('text')
    .attr('text-anchor','middle').attr('dy','1.2em')
    .style('font-family','var(--mono)').style('fill','#8a90a8')
    .style('font-size','10px').style('letter-spacing','0.06em').style('pointer-events','none')
    .text('378,253');

  function arcVisible(d){return d.y1<=4 && d.y0>=1 && d.x1>d.x0}
  function levelOf(d){return ['','phylum','class','order','family','genus'][d.depth]||'taxon'}

  function clicked(event,p){
    centerLabel.text(p.data.name==='RVMT'?'RVMT':p.data.name);
    centerSub.text(p.data.name==='RVMT'?'378,253':fmt(p.value));
    root.each(d=>d.target={
      x0:Math.max(0,Math.min(1,(d.x0-p.x0)/(p.x1-p.x0)))*2*Math.PI,
      x1:Math.max(0,Math.min(1,(d.x1-p.x0)/(p.x1-p.x0)))*2*Math.PI,
      y0:Math.max(0,d.y0-p.depth),
      y1:Math.max(0,d.y1-p.depth)
    });
    const t = svg.transition().duration(750);
    path.transition(t)
      .tween('data', d=>{const i=d3.interpolate(d.current,d.target);return t=>{d.current=i(t)}})
      .attr('fill-opacity', d=>arcVisible(d.target)? (d.children?0.85:0.65):0)
      .attr('pointer-events', d=>arcVisible(d.target)?'auto':'none')
      .attrTween('d', d=>()=>arc(d.current));
  }

  /* Side panel */
  const sbName=document.getElementById('sb-name'),sbLevel=document.getElementById('sb-level'),
        sbCount=document.getElementById('sb-count'),sbShare=document.getElementById('sb-share'),
        sbChild=document.getElementById('sb-children'),sbCrumb=document.getElementById('sb-crumb');
  function updateInfo(d){
    sbName.textContent = d.data.name === 'RVMT' ? 'All RVMT' : d.data.name;
    sbLevel.textContent = d.data.name === 'RVMT' ? 'root' : levelOf(d);
    sbCount.textContent = fmt(d.value);
    sbShare.textContent = fmtPct(d.value/totalAll);
    sbChild.textContent = d.children ? d.children.length : '0';
    const path = d.ancestors().reverse().map(a=>a.data.name).join(' › ');
    sbCrumb.textContent = '';
    path.split(' › ').forEach((part, i) => {
      if (i) sbCrumb.appendChild(document.createTextNode(' › '));
      if (part === 'RVMT'){
        const span = document.createElement('span');
        span.textContent = part;
        sbCrumb.appendChild(span);
      } else {
        sbCrumb.appendChild(document.createTextNode(part));
      }
    });
  }
  updateInfo(root);

  /* Legend */
  const leg = document.getElementById('phylum-legend');
  Object.entries(D.phylum_count)
    .filter(([k])=>!k.startsWith('p.'))
    .sort((a,b)=>b[1]-a[1])
    .forEach(([k,v])=>{
      leg.appendChild(legendItem(phyColor(k), k, fmt(v)));
    });
});

/* ---- Sankey ---- */
guard('host sankey', 'sankey', function(){
  const cont = document.getElementById('sankey');
  const W = cont.clientWidth || 600, H = 480;
  const svg = d3.select(cont).append('svg').attr('viewBox',`0 0 ${W} ${H}`).attr('width','100%').attr('height',H);

  const filtered = D.sankey.links.filter(l=>{
    const t = D.sankey.nodes[l.target];
    return !t.name.startsWith('p.');
  });
  const nodes = D.sankey.nodes.map(n=>({...n}));
  const links = filtered.map(l=>({...l}));

  const sankey = d3.sankey()
    .nodeAlign(d3.sankeyJustify)
    .nodeWidth(10).nodePadding(14)
    .extent([[8,18],[W-8,H-18]]);

  const g = sankey({nodes:nodes.map(d=>({...d})), links:links.map(d=>({...d}))});

  svg.append('g').selectAll('rect')
    .data(g.nodes).join('rect')
      .attr('x',d=>d.x0).attr('y',d=>d.y0)
      .attr('height',d=>d.y1-d.y0).attr('width',d=>d.x1-d.x0)
      .attr('fill',d=>d.group==='host'? '#7afcff' : phyColor(d.name))
      .attr('opacity',0.9);

  svg.append('g').attr('fill','none').selectAll('path')
    .data(g.links).join('path')
      .attr('class','sankey-link')
      .attr('d',d3.sankeyLinkHorizontal())
      .attr('stroke',d=>phyColor(d.target.name))
      .attr('stroke-opacity',0.28)
      .attr('stroke-width',d=>Math.max(1,d.width))
      .on('mouseover',(ev,d)=>{showTip(`<strong>${esc(d.source.name)} → ${esc(d.target.name)}</strong><div class="small">${fmt(d.value)} viruses</div>`,ev);d3.select(ev.currentTarget).attr('stroke-opacity',0.7)})
      .on('mouseout',(ev)=>{hideTip();d3.select(ev.currentTarget).attr('stroke-opacity',0.28)});

  svg.append('g').selectAll('text')
    .data(g.nodes).join('text')
      .attr('class','sankey-label')
      .attr('x',d=>d.x0<W/2 ? d.x1+8 : d.x0-8)
      .attr('y',d=>(d.y0+d.y1)/2)
      .attr('dy','0.35em')
      .attr('text-anchor',d=>d.x0<W/2 ? 'start' : 'end')
      .text(d=>d.name);

  svg.append('g').selectAll('text')
    .data(g.nodes).join('text')
      .attr('class','sankey-count')
      .attr('x',d=>d.x0<W/2 ? d.x1+8 : d.x0-8)
      .attr('y',d=>(d.y0+d.y1)/2 + 14)
      .attr('text-anchor',d=>d.x0<W/2 ? 'start' : 'end')
      .text(d=>fmt(d.value));
});

/* ---- Eukaryote treemap ---- */
guard('eukaryote treemap', 'euk-treemap', function(){
  const cont = document.getElementById('euk-treemap');
  const W = cont.clientWidth || 500, H = 480;
  const data = {name:'Eukaryota',children:Object.entries(D.euk_clade).map(([k,v])=>({name:k,value:v}))};
  const root = d3.hierarchy(data).sum(d=>d.value).sort((a,b)=>b.value-a.value);
  d3.treemap().size([W,H]).paddingInner(3).round(true)(root);
  const palette = ['#80ffae','#4cc9f0','#ffd166','#ff6b9d','#ff7a45','#ada0ff','#7afcff','#9d8cff','#ff5e9c','#67e8f9','#fbbf24'];
  const svg = d3.select(cont).append('svg').attr('viewBox',`0 0 ${W} ${H}`).attr('width','100%').attr('height',H);

  const cell = svg.selectAll('g').data(root.leaves()).join('g')
    .attr('transform',d=>`translate(${d.x0},${d.y0})`).attr('class','tm-cell');

  cell.append('rect')
    .attr('width',d=>d.x1-d.x0).attr('height',d=>d.y1-d.y0)
    .attr('fill',(d,i)=>palette[i%palette.length]).attr('rx',4)
    .on('mouseover',(ev,d)=>showTip(`<strong>${esc(d.data.name)}</strong><div class="small">${fmt(d.value)} viruses</div>`,ev))
    .on('mouseout',hideTip);

  cell.append('text')
    .attr('x',8).attr('y',18).text(d=>d.data.name)
    .each(function(d){if((d.x1-d.x0)<60||(d.y1-d.y0)<32)d3.select(this).remove()});
  cell.append('text').attr('class','v')
    .attr('x',8).attr('y',34).text(d=>fmt(d.value))
    .each(function(d){if((d.x1-d.x0)<60||(d.y1-d.y0)<46)d3.select(this).remove()});
});

/* ---- Novelty stacked bars ---- */
guard('novelty', 'novelty-chart', function(){
  const cont = document.getElementById('novelty-chart');
  const phyla = Object.keys(D.novelty).filter(p=>!p.startsWith('p.')).sort((a,b)=>{
    const sa=Object.values(D.novelty[a]).reduce((x,y)=>x+y,0);
    const sb=Object.values(D.novelty[b]).reduce((x,y)=>x+y,0);
    return sb-sa;
  });
  const W=cont.clientWidth||520, H=Math.max(300,phyla.length*48+40);
  const svg=d3.select(cont).append('svg').attr('viewBox',`0 0 ${W} ${H}`).attr('width','100%').attr('height',H);
  const margin={top:10,right:80,bottom:10,left:140};
  const max=d3.max(phyla,p=>Object.values(D.novelty[p]).reduce((x,y)=>x+y,0));
  const x=d3.scaleLinear().domain([0,max]).range([margin.left,W-margin.right]);
  const y=d3.scaleBand().domain(phyla).range([margin.top,H-margin.bottom]).padding(0.35);
  const g=svg.append('g');

  phyla.forEach(p=>{
    const counts = D.novelty[p];
    const novel = +counts['TRUE']||0;
    const ref = +counts['FALSE']||0;
    const total = novel+ref;
    const yy = y(p);
    g.append('text').attr('x',margin.left-12).attr('y',yy+y.bandwidth()/2).attr('dy','0.35em')
      .attr('text-anchor','end').style('fill','var(--text)').style('font-size','13px').text(p);
    g.append('rect').attr('x',x(0)).attr('y',yy).attr('height',y.bandwidth()).attr('width',x(novel)-x(0))
      .attr('fill',phyColor(p)).attr('opacity',0.9)
      .on('mouseover',ev=>showTip(`<strong>${esc(p)}</strong><div class="small">${fmt(novel)} novel · ${(novel/total*100).toFixed(1)}%</div>`,ev))
      .on('mouseout',hideTip);
    g.append('rect').attr('x',x(novel)).attr('y',yy).attr('height',y.bandwidth()).attr('width',x(novel+ref)-x(novel))
      .attr('fill',phyColor(p)).attr('opacity',0.28)
      .on('mouseover',ev=>showTip(`<strong>${esc(p)}</strong><div class="small">${fmt(ref)} reference</div>`,ev))
      .on('mouseout',hideTip);
    g.append('text').attr('x',W-margin.right+8).attr('y',yy+y.bandwidth()/2).attr('dy','0.35em')
      .style('fill','var(--muted)').style('font-size','11px').style('font-family','var(--mono)')
      .text(`${(novel/total*100).toFixed(0)}%`);
  });
});

/* ---- Segmentation ---- */
guard('segmentation', 'seg-chart', function(){
  const cont = document.getElementById('seg-chart');
  const phyla = Object.keys(D.segmentation).filter(p=>!p.startsWith('p.')).sort((a,b)=>{
    const sa=Object.values(D.segmentation[a]).reduce((x,y)=>x+y,0);
    const sb=Object.values(D.segmentation[b]).reduce((x,y)=>x+y,0);
    return sb-sa;
  });
  const W=cont.clientWidth||520, H=Math.max(300,phyla.length*48+40);
  const svg=d3.select(cont).append('svg').attr('viewBox',`0 0 ${W} ${H}`).attr('width','100%').attr('height',H);
  const margin={top:10,right:80,bottom:10,left:140};
  const max=d3.max(phyla,p=>Object.values(D.segmentation[p]).reduce((x,y)=>x+y,0));
  const x=d3.scaleLinear().domain([0,max]).range([margin.left,W-margin.right]);
  const y=d3.scaleBand().domain(phyla).range([margin.top,H-margin.bottom]).padding(0.35);
  const g=svg.append('g');
  phyla.forEach(p=>{
    const c=D.segmentation[p];
    const f=+c['FALSE']||0,t=+c['TRUE']||0,total=f+t;
    const yy=y(p);
    g.append('text').attr('x',margin.left-12).attr('y',yy+y.bandwidth()/2).attr('dy','0.35em')
      .attr('text-anchor','end').style('fill','var(--text)').style('font-size','13px').text(p);
    g.append('rect').attr('x',x(0)).attr('y',yy).attr('height',y.bandwidth()).attr('width',x(f)-x(0))
      .attr('fill',phyColor(p)).attr('opacity',0.85)
      .on('mouseover',ev=>showTip(`<strong>${esc(p)}</strong><div class="small">${fmt(f)} monopartite</div>`,ev)).on('mouseout',hideTip);
    g.append('rect').attr('x',x(f)).attr('y',yy).attr('height',y.bandwidth()).attr('width',x(f+t)-x(f))
      .attr('fill','#9d8cff').attr('opacity',0.85)
      .on('mouseover',ev=>showTip(`<strong>${esc(p)}</strong><div class="small">${fmt(t)} segmented</div>`,ev)).on('mouseout',hideTip);
    g.append('text').attr('x',W-margin.right+8).attr('y',yy+y.bandwidth()/2).attr('dy','0.35em')
      .style('fill','var(--muted)').style('font-size','11px').style('font-family','var(--mono)')
      .text(total?`${(t/total*100).toFixed(0)}% seg`:'—');
  });
  /* mini legend */
  const legend=svg.append('g').attr('transform',`translate(${margin.left},${H-2})`);
});

/* ---- Length small multiples ---- */
guard('genome length', 'length-charts', function(){
  const cont = document.getElementById('length-multiples');
  const phyla = Object.keys(D.length_hist).filter(p=>!p.startsWith('p.') && D.length_hist[p].edges.length);
  phyla.sort((a,b)=> (D.phylum_count[b]||0)-(D.phylum_count[a]||0));
  const cols=3, rows=Math.ceil(phyla.length/cols);
  const cellW=280, cellH=160, gap=24;
  const W = cellW*cols + gap*(cols-1) + 40;
  const H = cellH*rows + gap*(rows-1) + 40;
  const svg=d3.select(cont).append('svg').attr('viewBox',`0 0 ${W} ${H}`).attr('width','100%').attr('height',H);
  phyla.forEach((p,i)=>{
    const r=Math.floor(i/cols),c=i%cols;
    const tx = 20 + c*(cellW+gap), ty = 20 + r*(cellH+gap);
    const g=svg.append('g').attr('transform',`translate(${tx},${ty})`);
    const h=D.length_hist[p];
    const x = d3.scaleLog().domain([h.edges[0],h.edges[h.edges.length-1]]).range([0,cellW]);
    const y = d3.scaleLinear().domain([0,d3.max(h.counts)]).range([cellH-30,10]);
    g.append('text').attr('x',0).attr('y',-2).style('fill','var(--text)').style('font-size','13px')
      .style('font-family',"'Space Grotesk',sans-serif").text(p);
    g.append('text').attr('x',cellW).attr('y',-2).attr('text-anchor','end')
      .style('fill','var(--muted)').style('font-size','10px').style('font-family','var(--mono)')
      .text(`n=${fmt(D.phylum_count[p])}`);
    /* gridlines */
    [100,1000,10000].forEach(t=>{
      if(t>=h.edges[0]&&t<=h.edges[h.edges.length-1]){
        g.append('line').attr('class','gridline').attr('x1',x(t)).attr('x2',x(t)).attr('y1',10).attr('y2',cellH-30);
        g.append('text').attr('x',x(t)).attr('y',cellH-16).attr('text-anchor','middle')
          .style('fill','var(--dim)').style('font-size','9px').style('font-family','var(--mono)')
          .text(t>=1000?(t/1000)+'k':t+'');
      }
    });
    /* bars */
    g.selectAll('rect.b').data(h.counts).join('rect')
      .attr('x',(d,j)=>x(h.edges[j])+1)
      .attr('y',d=>y(d))
      .attr('width',(d,j)=>Math.max(1,x(h.edges[j+1])-x(h.edges[j])-1))
      .attr('height',d=>cellH-30-y(d))
      .attr('fill',phyColor(p)).attr('opacity',0.85)
      .on('mouseover',(ev,d,j)=>{const idx=h.counts.indexOf(d);showTip(`<strong>${esc(p)}</strong><div class="small">${Math.round(h.edges[idx])}–${Math.round(h.edges[idx+1])} nt · ${fmt(d)}</div>`,ev)})
      .on('mouseout',hideTip);
    /* axis baseline */
    g.append('line').attr('x1',0).attr('x2',cellW).attr('y1',cellH-30).attr('y2',cellH-30).attr('stroke','var(--border)');
  });
});

/* ---- HMM profiles ---- */
guard('protein profiles', 'clans-chart', function(){
  const cont = document.getElementById('profiles');
  const items = D.top_profiles.filter(p=>p.name!=='NA' && p.name!=='Unknown').slice(0,16);
  const max = d3.max(items, d=>d.count);
  items.forEach(p=>{
    const row = document.createElement('div'); row.className='profile-row';
    const name = document.createElement('div'); name.className='pname'; name.textContent=p.name;
    const bar = document.createElement('div'); bar.className='pbar';
    const total = Object.values(p.phylum_dist).reduce((a,b)=>a+b,0);
    const ordered = Object.entries(p.phylum_dist).sort((a,b)=>b[1]-a[1]);
    const widthPct = (p.count/max)*100;
    bar.style.width = widthPct+'%';
    ordered.forEach(([k,v])=>{
      const seg=document.createElement('span');
      seg.style.background=phyColor(k);
      seg.style.width=(v/total*100)+'%';
      seg.title=`${k}: ${fmt(v)}`;
      seg.addEventListener('mouseover',ev=>showTip(`<strong>${esc(p.name)}</strong><div class="small">${esc(k)} · ${fmt(v)} hits</div>`,ev));
      seg.addEventListener('mouseout',hideTip);
      bar.appendChild(seg);
    });
    const ct = document.createElement('div'); ct.className='pcount'; ct.textContent=fmt(p.count);
    row.appendChild(name); row.appendChild(bar); row.appendChild(ct);
    cont.appendChild(row);
  });
});

/* ---- Top families per phylum ---- */
guard('top families', 'fam-bars', function(){
  const phyla = Object.keys(D.top_families_per_phylum)
    .filter(p=>!p.startsWith('p.')&&p!=='Unclassified')
    .sort((a,b)=>(D.phylum_count[b]||0)-(D.phylum_count[a]||0));
  const pillsEl=document.getElementById('fam-pills'), out=document.getElementById('fam-bars');

  function render(p){
    out.textContent='';
    const data = D.top_families_per_phylum[p].slice(0,12);
    // The label gutter was a fixed 170px, so below about 220px of width the
    // bar scale ran backwards and d3 emitted negative rect widths. The gutter
    // now scales with the box and the drawing area can never be narrower than
    // a few pixels.
    const W = Math.max(out.clientWidth || 500, 240), H = data.length*30+30;
    const gutter = Math.min(170, Math.max(90, Math.round(W * 0.38)));
    const right = Math.max(gutter + 12, W - 50);
    const svg=d3.select(out).append('svg').attr('viewBox',`0 0 ${W} ${H}`).attr('width','100%').attr('height',H);
    const max=d3.max(data,d=>d.count);
    const x=d3.scaleLinear().domain([0,max]).range([gutter,right]);
    const y=d3.scaleBand().domain(data.map(d=>d.family)).range([10,H-10]).padding(0.3);
    data.forEach(d=>{
      const yy=y(d.family);
      svg.append('text').attr('x',gutter-10).attr('y',yy+y.bandwidth()/2).attr('dy','0.35em').attr('text-anchor','end')
        .style('fill','var(--text)').style('font-size','12px')
        .style('font-family','var(--mono)').text(d.family.length>22?d.family.slice(0,22)+'…':d.family);
      svg.append('rect').attr('x',gutter).attr('y',yy).attr('height',y.bandwidth()).attr('width',Math.max(0,x(d.count)-gutter))
        .attr('fill',phyColor(p)).attr('opacity',0.85).attr('rx',2)
        .on('mouseover',ev=>showTip(`<strong>${esc(d.family)}</strong><div class="small">${fmt(d.count)} viruses</div>`,ev))
        .on('mouseout',hideTip);
      svg.append('text').attr('x',x(d.count)+6).attr('y',yy+y.bandwidth()/2).attr('dy','0.35em')
        .style('fill','var(--muted)').style('font-size','11px').style('font-family','var(--mono)').text(fmt(d.count));
    });
  }

  phyla.forEach((p,i)=>{
    const b=document.createElement('button');b.className='pill'+(i===0?' active':'');b.textContent=p;
    b.onclick=()=>{document.querySelectorAll('#fam-pills .pill').forEach(el=>el.classList.remove('active'));b.classList.add('active');render(p)};
    pillsEl.appendChild(b);
  });
  if(phyla.length) render(phyla[0]);
});

/* ---- Donut helper ---- */
function drawDonut(containerId, data, opts={}){
  const cont=document.getElementById(containerId);
  const W=cont.clientWidth||320, H=300, R=Math.min(W,H)/2-10, r=R*0.62;
  const svg=d3.select(cont).append('svg').attr('viewBox',`${-W/2} ${-H/2} ${W} ${H}`).attr('width','100%').attr('height',H);
  const total=data.reduce((a,b)=>a+b.value,0);
  const pie=d3.pie().value(d=>d.value).sort(null).padAngle(0.012);
  const arc=d3.arc().innerRadius(r).outerRadius(R);
  const palette = opts.palette || ['#7afcff','#ff5e9c','#ffd166','#80ffae','#9d8cff','#ff7a45','#67e8f9','#a78bfa','#fbbf24','#f472b6'];
  svg.selectAll('path').data(pie(data)).join('path')
    .attr('d',arc)
    .attr('fill',(d,i)=>palette[i%palette.length])
    .attr('opacity',0.92)
    .on('mouseover',(ev,d)=>showTip(`<strong>${esc(d.data.name)}</strong><div class="small">${fmt(d.data.value)} (${(d.data.value/total*100).toFixed(1)}%)</div>`,ev))
    .on('mouseout',hideTip);
  const tg = svg.append('g').attr('class','donut-center');
  tg.append('text').attr('class','big').attr('y',-2).text(opts.bigLabel||fmt(total));
  tg.append('text').attr('class','small').attr('y',16).text(opts.smallLabel||'TOTAL');
}

/* Source */
guard('source donut', 'source-donut', function(){
  let raw = Object.entries(D.source).map(([k,v])=>({name:k.length>26?k.slice(0,26)+'…':k,full:k,value:v})).sort((a,b)=>b.value-a.value);
  /* Defensive: cap to top 8 + 'Other' so a schema mishap can never blow up the card */
  let items = raw;
  if (raw.length > 8) {
    const top = raw.slice(0,7);
    const other = raw.slice(7).reduce((a,b)=>a+b.value,0);
    items = [...top, {name:`Other (${raw.length-7})`, full:`Other (${raw.length-7} sources)`, value:other}];
  }
  const palette=['#7afcff','#ff5e9c','#ffd166','#80ffae','#9d8cff','#ff7a45','#67e8f9','#a78bfa'];
  drawDonut('source-donut',items,{bigLabel:fmt(items.reduce((a,b)=>a+b.value,0)),smallLabel:'CONTIGS',palette});
  const cont=document.getElementById('source-donut');
  const leg=document.createElement('div');leg.className='legend';leg.style.marginTop='12px';
  items.forEach((it,i)=>{
    const label = it.full.length > 34 ? it.full.slice(0, 34) + '…' : it.full;
    leg.appendChild(legendItem(palette[i % palette.length], label));
  });
  cont.appendChild(leg);
});

/* Genetic code */
guard('genetic code', 'gcode-donut', function(){
  const items = Object.entries(D.genetic_code).map(([k,v])=>({name:k==='NA'?'Unspecified':'Code '+k,value:v})).sort((a,b)=>b.value-a.value);
  drawDonut('gcode-donut',items,{smallLabel:'CODE TABLE',palette:['#80ffae','#7afcff','#ffd166','#ff6b9d','#9d8cff']});
  const cont=document.getElementById('gcode-donut');
  const leg=document.createElement('div');leg.className='legend';leg.style.marginTop='12px';
  const palette=['#80ffae','#7afcff','#ffd166','#ff6b9d','#9d8cff'];
  items.forEach((it, i) => leg.appendChild(legendItem(palette[i % palette.length], it.name)));
  cont.appendChild(leg);
});

/* Affiliation level bars */
guard('affiliation level', 'afl-bars', function(){
  const cont=document.getElementById('afl-bars');
  const items=Object.entries(D.afflvl).map(([k,v])=>({k,v})).sort((a,b)=>{
    const m=s=>{const m=s.match(/Lvl\s*(\d)/);return m?+m[1]:99};
    return m(a.k)-m(b.k);
  });
  const W=Math.max(cont.clientWidth||320, 240), H=Math.max(220,items.length*38);
  const svg=d3.select(cont).append('svg').attr('viewBox',`0 0 ${W} ${H}`).attr('width','100%').attr('height',H);
  const max=d3.max(items,d=>d.v);
  const x=barScale(max, 60, W, 50);
  const y=d3.scaleBand().domain(items.map(d=>d.k)).range([0,H]).padding(0.3);
  const palette=['#80ffae','#7afcff','#ffd166','#ff7a45','#ff6b9d','#9d8cff'];
  items.forEach((d,i)=>{
    const lbl = (d.k.match(/Lvl\s*(\d)/)||[null,'?'])[1];
    const yy=y(d.k);
    svg.append('text').attr('x',50).attr('y',yy+y.bandwidth()/2).attr('dy','0.35em').attr('text-anchor','end')
      .style('fill','var(--text)').style('font-size','12px').style('font-family','var(--mono)').text('Lvl '+lbl);
    svg.append('rect').attr('x',60).attr('y',yy).attr('height',y.bandwidth()).attr('width',Math.max(0,x(d.v)-60))
      .attr('fill',palette[i%palette.length]).attr('opacity',0.88).attr('rx',2)
      .on('mouseover',ev=>showTip(`<strong>${esc(d.k.split('-')[0].trim())}</strong><div class="small">${fmt(d.v)} viruses<br>${esc(d.k.split('-').slice(1).join('-').trim())}</div>`,ev))
      .on('mouseout',hideTip);
    svg.append('text').attr('x',x(d.v)+6).attr('y',yy+y.bandwidth()/2).attr('dy','0.35em')
      .style('fill','var(--muted)').style('font-size','11px').style('font-family','var(--mono)').text(fmt(d.v));
  });
});

/* ---- Top protein clans ---- */
guard('protein clans', 'clans-chart', function(){
  const cont = document.getElementById('clans-chart');
  if(!cont || !D.top_clans) return;
  const items = D.top_clans.filter(d=>d.name && d.name!=='Unknown' && d.name!=='NA').slice(0,20);
  if(!items.length) return;
  const W = Math.max(cont.clientWidth || 1000, 300), rowH = 30, H = items.length*rowH + 30;
  const svg = d3.select(cont).append('svg').attr('viewBox',`0 0 ${W} ${H}`).attr('width','100%').attr('height',H);
  const max = d3.max(items, d=>d.count);
  const x = barScale(max, 130, W, 90);
  const y = d3.scaleBand().domain(items.map(d=>d.name)).range([10, H-20]).padding(0.28);
  const palette = ['#7afcff','#9d8cff','#ff6b9d','#ffd166','#80ffae','#ff7a45','#67e8f9','#a78bfa','#fbbf24','#f472b6'];
  /* gridlines */
  const ticks = x.ticks(5);
  ticks.forEach(t=>{
    svg.append('line').attr('class','gridline').attr('x1',x(t)).attr('x2',x(t)).attr('y1',8).attr('y2',H-18);
    svg.append('text').attr('x',x(t)).attr('y',H-4).attr('text-anchor','middle')
      .style('fill','var(--dim)').style('font-size','10px').style('font-family','var(--mono)')
      .text(t>=1000?(t/1000).toFixed(0)+'k':t);
  });
  items.forEach((d,i)=>{
    const yy = y(d.name);
    svg.append('text').attr('x',120).attr('y',yy+y.bandwidth()/2).attr('dy','0.35em').attr('text-anchor','end')
      .style('fill','var(--text)').style('font-size','12px').style('font-family','var(--mono)')
      .text(d.name);
    svg.append('rect').attr('x',130).attr('y',yy).attr('height',y.bandwidth()).attr('width',Math.max(0,x(d.count)-130))
      .attr('fill',palette[i%palette.length]).attr('opacity',0.85).attr('rx',3)
      .on('mouseover',ev=>showTip(`<strong>${esc(d.name)}</strong><div class="small">${fmt(d.count)} HMM hits</div>`,ev))
      .on('mouseout',hideTip);
    svg.append('text').attr('x',x(d.count)+6).attr('y',yy+y.bandwidth()/2).attr('dy','0.35em')
      .style('fill','var(--muted)').style('font-size','11px').style('font-family','var(--mono)')
      .text(fmt(d.count));
  });
});

