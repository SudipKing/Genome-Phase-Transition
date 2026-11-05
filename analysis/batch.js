// Runs many seeds / parameter sets in parallel with run.js.
//
// usage: node analysis/batch.js plan.json [--jobs 4]
//
// plan.json:
// {
//   "out": "analysis/runs/batch1",              // one sub-directory per run
//   "args": ["--ticks", "200000", "--settleFull", "15", "--historyTicks", "0"],
//   "sets": [
//     {"name": "base", "seeds": [201, 202], "params": {}},
//     {"name": "mut0.02", "seeds": [301, 302], "params": {"mutationRate": 0.02}}
//   ]
// }
//
// Each run is written to <out>/<set name>-s<seed>; runs whose run.json already exists are
// skipped, so an interrupted batch can be restarted. batch.json in <out> lists every run.

const fs = require('fs');
const path = require('path');
const {spawn} = require('child_process');

const planFile = process.argv[2];
if (!planFile) { console.error('usage: node batch.js plan.json [--jobs N]'); process.exit(1); }
const ji = process.argv.indexOf('--jobs');
const jobs = ji > 0 ? Number(process.argv[ji + 1]) : Math.max(1, require('os').cpus().length);
const plan = JSON.parse(fs.readFileSync(planFile, 'utf8'));
const out = path.resolve(plan.out);
fs.mkdirSync(out, {recursive: true});

const queue = [];
for (const set of plan.sets)
    for (const seed of set.seeds)
        queue.push({set: set.name, seed, params: set.params || {}, dir: path.join(out, `${set.name}-s${seed}`)});
fs.writeFileSync(path.join(out, 'batch.json'), JSON.stringify({plan, runs: queue.map(q => ({...q, dir: path.basename(q.dir)}))}, null, 2));

let next = 0, active = 0, done = 0;
const t0 = Date.now();
function launch() {
    while (active < jobs && next < queue.length) {
        const job = queue[next++];
        if (fs.existsSync(path.join(job.dir, 'run.json'))) { done++; continue; }
        active++;
        const argv = [path.join(__dirname, 'run.js'), '--seed', String(job.seed), '--out', job.dir,
            '--params', JSON.stringify(job.params), ...(plan.args || [])];
        const p = spawn(process.execPath, argv, {stdio: ['ignore', 'pipe', 'inherit']});
        let tail = '';
        p.stdout.on('data', d => { tail = (tail + d).slice(-400); });
        p.on('close', code => {
            active--; done++;
            const lastLine = tail.trim().split('\n').pop();
            console.log(`[${done}/${queue.length} ${((Date.now() - t0) / 60000).toFixed(1)}m] ${job.set} s${job.seed} exit ${code}: ${lastLine}`);
            launch();
        });
    }
}
launch();
