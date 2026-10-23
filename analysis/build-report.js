// Builds the shareable report page: template.html + rd.json (from analyze.js) + findings.json.
// usage: node build-report.js rd.json findings.json report.html
const fs = require('fs');
const [rdPath, findingsPath, outPath] = process.argv.slice(2);
const rd = JSON.parse(fs.readFileSync(rdPath, 'utf8'));
const findings = JSON.parse(fs.readFileSync(findingsPath, 'utf8'));
// Colour-relabeling comparison from symmetry.js, if present beside rd.json.
const symPath = require('path').join(require('path').dirname(rdPath), 'symmetry.json');
const sym = fs.existsSync(symPath) ? JSON.parse(fs.readFileSync(symPath, 'utf8')) : null;

// Trim to what the page draws.
const data = {
    conditions: rd.conditions,
    crossRun: rd.crossRun,
    symmetry: sym && {pairs: sym.pairs, rules: sym.rules, canonical: sym.canonical},
    runs: rd.runs.map(r => ({
        name: r.name, seed: r.seed, ticks: r.ticks, phaseSequence: r.phaseSequence,
        lastFounderTick: r.lastFounderTick, spikes: r.spikes,
        phases: r.phases.map(({enc, ...p}) => p),
        finalEnc: r.phases[r.phases.length - 1].enc,
        series: {...r.series, consensusSize: r.consensusSize, consensusAgreement: r.consensusAgreement, meanDistFromConsensus: r.meanDistFromConsensus},
        strip: r.domGenomes,
    })),
};
const tpl = fs.readFileSync(__dirname + '/template.html', 'utf8');
const html = tpl
    .replace('/*__DATA__*/null', JSON.stringify(data))
    .replace('/*__FINDINGS__*/null', JSON.stringify(findings));
fs.writeFileSync(outPath, html);
console.log('wrote', outPath, (html.length / 1024).toFixed(0) + ' KB');
