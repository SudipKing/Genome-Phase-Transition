// Local, database-independent storage for simulation run outputs.
// JSON Lines keeps snapshots append-only and allows later streaming into a database.
// By default each line is written as its own gzip member: concatenated members are one
// valid .gz stream, so the file stays append-only and partial runs remain readable,
// while genome strings shrink about 10x.
const fs = require('fs');
const path = require('path');
const zlib = require('zlib');

const PLAIN = 'genome-snapshots.jsonl';
const GZ = PLAIN + '.gz';

function createRun(directory, {gzip = true} = {}) {
    fs.mkdirSync(directory, {recursive: true});
    const snapshotFile = gzip ? GZ : PLAIN;
    const snapshotPath = path.join(directory, snapshotFile);
    fs.writeFileSync(snapshotPath, '');
    // Don't leave a stale file of the other format beside the new one.
    const other = path.join(directory, gzip ? PLAIN : GZ);
    if (fs.existsSync(other)) fs.unlinkSync(other);
    return {directory, snapshotPath, snapshotFile, gzip};
}

function saveSnapshot(run, snapshot) {
    const line = JSON.stringify(snapshot) + String.fromCharCode(10);
    fs.appendFileSync(run.snapshotPath, run.gzip ? zlib.gzipSync(line, {level: 9}) : line);
}

function completeRun(run, metadata) {
    fs.writeFileSync(path.join(run.directory, 'run.json'), JSON.stringify(metadata, null, 2));
}

// All snapshots of a run directory, oldest first, from either format.
function readSnapshots(directory) {
    const gz = path.join(directory, GZ), plain = path.join(directory, PLAIN);
    const text = fs.existsSync(gz) ? zlib.gunzipSync(fs.readFileSync(gz)).toString('utf8')
        : fs.readFileSync(plain, 'utf8');
    return text.split(String.fromCharCode(10)).filter(l => l.trim()).map(l => JSON.parse(l));
}

module.exports = {createRun, saveSnapshot, completeRun, readSnapshots};
