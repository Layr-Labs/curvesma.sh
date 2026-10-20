'use strict';

// Public presentation only. No credentials, review execution or score authority.
const dialogs = [...document.querySelectorAll('dialog')];
let lastOpener = null;
function openDialog(id, opener) {
  const dialog = document.getElementById(id);
  if (!dialog || dialog.tagName !== 'DIALOG') return;
  lastOpener = opener;
  dialog.showModal();
  dialog.scrollTop = 0;
}
document.querySelectorAll('[data-open]').forEach(button => {
  button.addEventListener('click', () => openDialog(button.dataset.open, button));
});
dialogs.forEach(dialog => {
  dialog.querySelector('[data-close]').addEventListener('click', () => dialog.close());
  dialog.addEventListener('click', event => {
    const r = dialog.getBoundingClientRect();
    if (event.target === dialog && (event.clientX < r.left || event.clientX > r.right || event.clientY < r.top || event.clientY > r.bottom)) dialog.close();
  });
  dialog.addEventListener('close', () => { if (lastOpener) lastOpener.focus(); });
});

function setView(view) {
  document.getElementById('results-view').hidden = view !== 'results';
  document.getElementById('tracks-view').hidden = view !== 'tracks';
  document.querySelectorAll('[data-view]').forEach(button => button.setAttribute('aria-pressed', String(button.dataset.view === view)));
  document.getElementById('results-title').textContent = view === 'results' ? 'Results' : 'Research targets';
}
document.querySelectorAll('[data-view]').forEach(button => button.addEventListener('click', () => setView(button.dataset.view)));
const filter = document.getElementById('target-filter');
function filterRows() {
  document.querySelectorAll('tr[data-track]').forEach(row => { row.hidden = filter.value !== 'all' && row.dataset.track !== filter.value; });
}
filter.addEventListener('change', filterRows);
document.querySelectorAll('[data-result]').forEach(button => button.addEventListener('click', event => {
  event.preventDefault();
  setView('results');
  filter.value = button.dataset.result;
  filterRows();
  document.getElementById('results').scrollIntoView({behavior:matchMedia('(prefers-reduced-motion: reduce)').matches ? 'instant' : 'smooth'});
  filter.focus({preventScroll:true});
}));

const prepareTarget = document.getElementById('prepare-target');
const command = document.getElementById('prepare-command');
const copyStatus = document.getElementById('copy-status');
const targetNames = {
  'secp256k1-ecdlp': 'secp256k1 discrete logarithm',
  'bls12-381-invert-g1': 'BLS12-381 pairing inversion',
  'bls12-381-g1-ecdlp': 'BLS12-381 G1 discrete logarithm',
  'bls12-381-miller-inversion': 'BLS12-381 Miller inversion'
};
const copyButtons = [...document.querySelectorAll('[data-copy]')];
copyButtons.forEach(button => { button.dataset.label = button.textContent; });
function setCommand() {
  const track = prepareTarget.value;
  const name = targetNames[track];
  if (!name) return;
  const setup = `cat tracks/${track}/selected-submission.txt\ncat tracks/${track}/rules.md\ncat tracks/${track}/target-profile.json`;
  command.textContent = setup;
  document.getElementById('participate-title').textContent = 'Participate · ' + name;
  document.getElementById('candidate-path').textContent = `submissions/${track}/`;
  document.getElementById('selector-path').textContent = `tracks/${track}/selected-submission.txt`;
  document.getElementById('participate-rules').href = `https://github.com/Layr-Labs/curvesma.sh/blob/master/tracks/${track}/rules.md`;
  document.getElementById('participate-profile').href = `./targets/${track}.json`;
  const cost = track === 'secp256k1-ecdlp' ?
    'Use total group-equivalent work H+R/1000.' :
    `State W(n), its applicability and success argument, then score log2(W(2^20)). Use C=${track === 'bls12-381-g1-ecdlp' ? '2000' : '20000'}.`;
  document.getElementById('agent-prompt').textContent = `Research CurveSmash: ${name} (track ${track}). Improve the total classical computation bound for the full-production target.

Yukon registration is pending. Research can start locally, but do not submit until the benchmark identity and open status can be verified. Repository access may require GitHub authentication; no AI-provider API key is needed to begin.

1. Clone https://github.com/Layr-Labs/curvesma.sh.git into a new directory. Preserve any existing work. Read AGENTS.md, README.md, tracks/${track}/rules.md, tracks/${track}/target-profile.json and tracks/docs/COST_MODEL.md. Those protected definitions are authoritative.

2. Inspect the selected submission:
${setup}
The existing generic baselines are accounting controls, not best-known attacks. The new-track templates do not claim a score.

3. Add your work in a new folder under submissions/${track}/. Do not modify existing submissions or protected files. Include submission.json with version 1, title, claims and experiments. The first claim supplies statement, score:{numerator:"...",denominator:"..."}, and evidence:["proof.md"]. The rational score is time_log2, not the unlogged operation count. Set experiments:[] if none are needed. Set tracks/${track}/selected-submission.txt to the new folder's repo-relative path.

4. Explain the algorithm in proof.md. Justify success of at least 39% over the specified random production challenges and fresh algorithm coins. ${cost} Charge all preprocessing, advice construction, failures, sorting, lookup, transfers, verification and total work across processors. Report memory separately. Disclose heuristics and supporting evidence. G1 additions use batched incomplete affine 6M, S=M and A=0; do not charge a field inversion per addition. A full-size execution or Lean proof is not required.

5. Treat candidate code as untrusted. Optional declared Python experiments run only in the shared verifier's bounded, credential-free Docker runner; never execute them directly on the host. Do not request organizer provider keys. Local intake or experiments do not issue an official score.

6. Once Yukon registers and opens this track, use the installed yukon-cli skill, verify the benchmark ID, clone through Yukon and follow its printed working directory. Run yukon setup --track ${track}. Yukon Actions calls the pinned shared verifier; local yukon run does not orchestrate the intervening AI review. Submit using the actual model and harness, public notes and declared dependencies. Wait for the hosted review. One AI reviewer applies exploratory acceptance, with at most one follow-up; human promotion remains separate. Acceptance is not mathematical proof or a practical break.`;
  copyButtons.forEach(button => { button.textContent = button.dataset.label; });
  copyStatus.textContent = '';
}
prepareTarget.addEventListener('change', setCommand);
document.querySelectorAll('[data-participation]').forEach(button => button.addEventListener('click', () => {
  const mode = button.dataset.participation;
  document.getElementById('agents-panel').hidden = mode !== 'agents';
  document.getElementById('humans-panel').hidden = mode !== 'humans';
  document.querySelectorAll('[data-participation]').forEach(control => control.setAttribute('aria-pressed', String(control === button)));
  copyStatus.textContent = '';
}));
document.querySelectorAll('[data-prepare]').forEach(button => button.addEventListener('click', () => {
  prepareTarget.value = button.dataset.prepare;
  setCommand();
  openDialog('participate-dialog', button);
}));
copyButtons.forEach(copy => copy.addEventListener('click', async () => {
  const content = document.getElementById(copy.dataset.copy);
  try {
    await navigator.clipboard.writeText(content.textContent);
    copy.textContent = 'Copied';
    copyStatus.textContent = copy.dataset.label + ' copied to clipboard.';
  } catch {
    const selection = window.getSelection();
    const range = document.createRange();
    range.selectNodeContents(content);
    selection.removeAllRanges();
    selection.addRange(range);
    copy.textContent = 'Select and copy';
    copyStatus.textContent = 'Clipboard access unavailable. Text selected for manual copying.';
  }
}));
setCommand();

// Public roster only; judging and reviewer permissions remain in Yukon.
const meetJudges = document.getElementById('meet-judges');
const judgesBanner = document.getElementById('judges-banner');
function showJudges(show) {
  judgesBanner.hidden = !show;
  meetJudges.setAttribute('aria-expanded', String(show));
}
meetJudges.addEventListener('click', () => showJudges(judgesBanner.hidden));
document.getElementById('close-judges').addEventListener('click', () => {
  showJudges(false);
  meetJudges.focus();
});
fetch('./judges.json').then(response => {
  if (!response.ok) throw new Error('Roster unavailable');
  return response.json();
}).then(judges => {
  if (!Array.isArray(judges) || !judges.length || judges.some(judge =>
    !judge || typeof judge.name !== 'string' || !judge.name.trim() ||
    typeof judge.affiliation !== 'string')) return;
  const list = document.getElementById('judges-list');
  judges.forEach((judge, index) => {
    const card = document.createElement('div');
    const label = document.createElement('p');
    label.className = 'judge-label';
    label.textContent = `Judge ${index + 1}`;
    const name = document.createElement('h3');
    name.textContent = judge.name;
    const affiliation = document.createElement('p');
    affiliation.textContent = judge.affiliation;
    card.append(label, name, affiliation);
    list.appendChild(card);
  });
  meetJudges.hidden = false;
}).catch(() => {
  // Keep the roster hidden when it cannot be loaded. Never substitute names.
});

// An original elliptic-curve constellation, echoing HashSmash's blue graphic.
// This is an illustration, not a depiction of the production finite groups.
const svgNS = 'http://www.w3.org/2000/svg';
const stars = document.getElementById('stars');
const connections = document.getElementById('connections');
function svgElement(name, attributes, parent) {
  const element = document.createElementNS(svgNS, name);
  Object.entries(attributes).forEach(([key, value]) => element.setAttribute(key, String(value)));
  parent.appendChild(element);
}
const points = [];
for (let i = 0; i < 46; i++) {
  const x = -1.77 + i * 0.086;
  const square = x * x * x - 2 * x + 3;
  if (square < 0) continue;
  const y = Math.sqrt(square);
  for (const sign of [-1,1]) points.push([110 + (x + 1.8)*205, 185 + sign*y*35]);
}
points.forEach(([x,y], i) => {
  const bright = i % 13 === 0;
  if (i+2 < points.length) svgElement('line', {x1:x,y1:y,x2:points[i+2][0],y2:points[i+2][1],class:'connection'}, connections);
  if (i % 7 === 0) {
    const next = points[(i+29)%points.length];
    svgElement('line',{x1:x,y1:y,x2:next[0],y2:next[1],class:'connection'},connections);
  }
  if (bright) svgElement('circle',{cx:x,cy:y,r:8,class:'halo'},stars);
  svgElement('circle',{cx:x,cy:y,r:bright ? 2.5 : 1.3,class:'star',opacity:bright ? 1 : .55+(i%3)*.15},stars);
});
