// Phase segmentation shared by analyze.js and deep.js.
// A phase is a stretch of samples with the same set of colours above PRESENT of living
// cells, split again when the grid is at least 95% full; blips shorter than MIN_PHASE
// samples merge into a neighbouring phase.

const COLOR_NAMES = 'RGB';
const PRESENT = 0.08;     // a colour "is present" in a phase above this share of living cells
const MIN_PHASE = 5;      // samples; shorter blips merge into the neighbouring phase

function label(s) {
    if (s.pop < 0.01) return 'extinct';
    const colours = s.colors.map((f, i) => f > PRESENT ? COLOR_NAMES[i] : '').join('') || '?';
    // Same colours can be a churning population or a filled grid; those are different phases.
    return colours + (s.pop >= 0.95 ? ' full' : '');
}

function segment(samples) {
    let raw = [];
    for (let i = 0; i < samples.length; i++) {
        const l = label(samples[i]);
        if (raw.length && raw[raw.length - 1].label === l) raw[raw.length - 1].end = i;
        else raw.push({label: l, start: i, end: i});
    }
    // Merge short blips into the previous phase until stable.
    let changed = true;
    while (changed && raw.length > 1) {
        changed = false;
        for (let k = 0; k < raw.length; k++) {
            const p = raw[k];
            if (p.end - p.start + 1 >= MIN_PHASE) continue;
            const into = k > 0 ? k - 1 : k + 1;
            raw[into].start = Math.min(raw[into].start, p.start);
            raw[into].end = Math.max(raw[into].end, p.end);
            raw.splice(k, 1);
            changed = true;
            break;
        }
        // Re-join neighbours that now share a label.
        for (let k = raw.length - 1; k > 0; k--) {
            if (raw[k].label === raw[k - 1].label) { raw[k - 1].end = raw[k].end; raw.splice(k, 1); changed = true; }
        }
    }
    return raw;
}

module.exports = {label, segment, PRESENT, MIN_PHASE, COLOR_NAMES};
