// Builds the interactive Phase Explorer from deep.js output.
//
// usage: node analysis/explorer.js deep.json explorer.html
//
// The page is self-contained (data embedded, no network): an overview of every run's
// phases, per-phase rule composition and genome similarity, event timing and fill-aligned
// metrics, the parameter sweep, and a per-run view with synced charts and the consensus
// genome strip.

const fs = require('fs');
const path = require('path');
const [inFile, outFile] = process.argv.slice(2);
if (!outFile) { console.error('usage: node explorer.js deep.json explorer.html'); process.exit(1); }
const data = fs.readFileSync(inFile, 'utf8');
JSON.parse(data);   // fail early on a truncated file
const tpl = fs.readFileSync(path.join(__dirname, 'explorer-template.html'), 'utf8');
// Escape "</" so the data cannot close the script element.
fs.writeFileSync(outFile, tpl.replace('/*DATA*/null', () => data.replace(/<\//g, '<\\/')));
console.error(`wrote ${outFile} (${(fs.statSync(outFile).size / 1e6).toFixed(1)} MB)`);
