// Local, database-independent storage for simulation run outputs.
// JSON Lines keeps snapshots append-only and allows later streaming into a database.
const fs = require('fs');
const path = require('path');

function createRun(directory) {
    fs.mkdirSync(directory, {recursive: true});
    const snapshotPath = path.join(directory, 'genome-snapshots.jsonl');
    fs.writeFileSync(snapshotPath, '');
    return {directory, snapshotPath};
}

function saveSnapshot(run, snapshot) {
    fs.appendFileSync(run.snapshotPath, JSON.stringify(snapshot) + '\\n');
}

function completeRun(run, metadata) {
    fs.writeFileSync(path.join(run.directory, 'run.json'), JSON.stringify(metadata, null, 2));
}

module.exports = {createRun, saveSnapshot, completeRun};
