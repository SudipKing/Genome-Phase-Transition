// ─── Organism log ─────────────────────────────────────────────────────────────
// Permanent record of every organism: parent, founder, birth/death tick, and its
// genome stored as the rules that differ from its parent's genome.
//
// Why deltas: across a 10k-tick run, 98% of ~1.8M organisms have a genome no other
// organism shares, so a genome -> [ids] map saves almost nothing. A child differs
// from its parent in ~4 of 252 rules, so storing only those changes is ~30x smaller,
// and any genome is rebuilt by replaying changes from the founder down the tree.
// Founders store their full genome as changes from an empty (all -1) genome.
//
// Columns are growable typed arrays, so millions of organisms stay compact in memory.
// File format "hex-ca-lineage-genomes-v1" (little-endian):
//   "HXLG" | u32 version | u32 headerBytes | header JSON (UTF-8) | pad to 4 | sections
// header.sections lists {name, type, length, offset} for each typed-array section.

class OrganismLog {
    constructor(ncond) {
        this.ncond = ncond;
        this.length = 0;
        this.deltaLength = 0;
        this.enabled = true;       // when false, new organisms are not recorded (ids still advance)
        const cap = 1 << 16;
        this.parent = new Int32Array(cap);
        this.founder = new Int32Array(cap);
        this.birth = new Int32Array(cap);
        this.death = new Int32Array(cap);
        this.deltaStart = new Uint32Array(cap + 1);
        this.deltaCond = new Uint16Array(cap * 4);
        this.deltaVal = new Int8Array(cap * 4);
    }

    _growRows() {
        const cap = this.parent.length * 2;
        for (const k of ['parent', 'founder', 'birth', 'death']) {
            const a = new Int32Array(cap); a.set(this[k]); this[k] = a;
        }
        const ds = new Uint32Array(cap + 1); ds.set(this.deltaStart); this.deltaStart = ds;
    }

    _growDeltas(need) {
        let cap = this.deltaCond.length;
        while (cap < need) cap *= 2;
        const c = new Uint16Array(cap); c.set(this.deltaCond); this.deltaCond = c;
        const v = new Int8Array(cap); v.set(this.deltaVal); this.deltaVal = v;
    }

    // Records an organism whose genome is `genome` (an Int8Array view of ncond rules).
    // `parentGenome` is the parent's genome view, or null for a founder.
    // Returns the id, which must equal the caller's own organism counter.
    add(id, parentId, founderId, tick, genome, parentGenome) {
        if (!this.enabled) return id;
        if (id !== this.length) throw new Error('OrganismLog ids must be sequential (got ' + id + ', expected ' + this.length + ')');
        if (this.length + 1 >= this.parent.length) this._growRows();
        const ncond = this.ncond;
        if (this.deltaLength + ncond > this.deltaCond.length) this._growDeltas(this.deltaLength + ncond);
        let d = this.deltaLength;
        for (let c = 0; c < ncond; c++) {
            const g = genome[c];
            if (parentGenome ? g !== parentGenome[c] : g !== -1) {
                this.deltaCond[d] = c;
                this.deltaVal[d] = g;
                d++;
            }
        }
        this.parent[id] = parentId;
        this.founder[id] = founderId;
        this.birth[id] = tick;
        this.death[id] = -1;
        this.deltaStart[id] = this.deltaLength;
        this.deltaLength = d;
        this.deltaStart[id + 1] = d;
        this.length = id + 1;
        return id;
    }

    // Records the first death only.
    recordDeath(id, tick) {
        if (id >= 0 && id < this.length && this.death[id] === -1) this.death[id] = tick;
    }

    // Serialises to an ArrayBuffer. `extra` is merged into the JSON header.
    toArrayBuffer(extra) {
        const n = this.length, m = this.deltaLength;
        const sections = [
            ['parent', 'Int32Array', this.parent.subarray(0, n)],
            ['founder', 'Int32Array', this.founder.subarray(0, n)],
            ['birth', 'Int32Array', this.birth.subarray(0, n)],
            ['death', 'Int32Array', this.death.subarray(0, n)],
            ['deltaStart', 'Uint32Array', this.deltaStart.subarray(0, n + 1)],
            ['deltaCond', 'Uint16Array', this.deltaCond.subarray(0, m)],
            ['deltaVal', 'Int8Array', this.deltaVal.subarray(0, m)],
        ];
        let offset = 0;
        const meta = sections.map(([name, type, arr]) => {
            const s = {name, type, length: arr.length, offset};
            offset += arr.byteLength;
            offset = (offset + 3) & ~3;        // keep every section 4-byte aligned
            return s;
        });
        const header = Object.assign({format: 'hex-ca-lineage-genomes-v1', ncond: this.ncond,
            organisms: n, deltas: m, sections: meta}, extra || {});
        const headerBytes = new TextEncoder().encode(JSON.stringify(header));
        const start = (12 + headerBytes.length + 3) & ~3;
        const buf = new ArrayBuffer(start + offset);
        const u8 = new Uint8Array(buf);
        u8.set([72, 88, 76, 71], 0);             // "HXLG"
        const dv = new DataView(buf);
        dv.setUint32(4, 1, true);
        dv.setUint32(8, headerBytes.length, true);
        u8.set(headerBytes, 12);
        sections.forEach(([, , arr], i) => u8.set(new Uint8Array(arr.buffer, arr.byteOffset, arr.byteLength), start + meta[i].offset));
        return buf;
    }

    // Parses a buffer written by toArrayBuffer. Returns {header, parent, founder, birth,
    // death, deltaStart, deltaCond, deltaVal, genomeOf(id)} with typed arrays over the buffer.
    static parse(buf) {
        const u8 = new Uint8Array(buf);
        if (u8[0] !== 72 || u8[1] !== 88 || u8[2] !== 76 || u8[3] !== 71) throw new Error('Not a lineage + genomes file (missing HXLG marker).');
        const dv = new DataView(buf);
        const headerLen = dv.getUint32(8, true);
        const header = JSON.parse(new TextDecoder().decode(u8.subarray(12, 12 + headerLen)));
        const start = (12 + headerLen + 3) & ~3;
        const types = {Int32Array, Uint32Array, Uint16Array, Int8Array};
        const out = {header};
        for (const s of header.sections) out[s.name] = new types[s.type](buf, start + s.offset, s.length);
        const ncond = header.ncond;
        // Rebuilds one organism's genome by replaying deltas from its founder down to it.
        out.genomeOf = (id, into) => {
            const g = into || new Int8Array(ncond);
            g.fill(-1);
            const chain = [];
            for (let a = id; a >= 0; a = out.parent[a]) chain.push(a);
            for (let k = chain.length - 1; k >= 0; k--) {
                const a = chain[k];
                for (let d = out.deltaStart[a]; d < out.deltaStart[a + 1]; d++) g[out.deltaCond[d]] = out.deltaVal[d];
            }
            return g;
        };
        return out;
    }
}

if (typeof module !== 'undefined') module.exports = OrganismLog;
