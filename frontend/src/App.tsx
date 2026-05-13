import React, { useState, useRef, useEffect, useCallback, useMemo } from "react";
import "./App.css";

const API_URL = process.env.REACT_APP_API_URL || "http://localhost:8000";

// ---------------------------------------------------------------------------
//  Types
// ---------------------------------------------------------------------------
interface ResidueAnnotation {
  residue: string;
  position: number;
  grn: string;
}

interface AnnotateResponse {
  name: string;
  sequence_length: number;
  annotations: ResidueAnnotation[];
  columns: string[];
  values: string[];
  missing_count: number;
}

interface ReferenceRow {
  name: string;
  values: string[];
}

interface ReferenceTable {
  columns: string[];
  rows: ReferenceRow[];
}

// ---------------------------------------------------------------------------
//  Sections — ordered list of all structural regions
// ---------------------------------------------------------------------------
interface Section {
  key: string;
  label: string;
  color: string;
  match: (grn: string) => boolean;
}

const SECTIONS: Section[] = [
  { key: "N",   label: "N-term", color: "#999",    match: (g) => g.startsWith("n.") },
  { key: "TM1", label: "TM1",    color: "#e63946", match: (g) => grnInHelix(g, 1) },
  { key: "L12", label: "ICL1",   color: "#6b7280", match: (g) => grnInLoop(g, 1, 2) },
  { key: "TM2", label: "TM2",    color: "#f4a261", match: (g) => grnInHelix(g, 2) },
  { key: "L23", label: "ECL1",   color: "#6b7280", match: (g) => grnInLoop(g, 2, 3) },
  { key: "TM3", label: "TM3",    color: "#2a9d8f", match: (g) => grnInHelix(g, 3) },
  { key: "L34", label: "ICL2",   color: "#6b7280", match: (g) => grnInLoop(g, 3, 4) },
  { key: "TM4", label: "TM4",    color: "#264653", match: (g) => grnInHelix(g, 4) },
  { key: "L45", label: "ECL2",   color: "#6b7280", match: (g) => grnInLoop(g, 4, 5) },
  { key: "TM5", label: "TM5",    color: "#e76f51", match: (g) => grnInHelix(g, 5) },
  { key: "L56", label: "ICL3",   color: "#6b7280", match: (g) => grnInLoop(g, 5, 6) },
  { key: "TM6", label: "TM6",    color: "#457b9d", match: (g) => grnInHelix(g, 6) },
  { key: "L67", label: "ECL3",   color: "#6b7280", match: (g) => grnInLoop(g, 6, 7) },
  { key: "TM7", label: "TM7",    color: "#6a0572", match: (g) => grnInHelix(g, 7) },
  { key: "C",   label: "C-term", color: "#999",    match: (g) => g.startsWith("c.") },
];

function grnInHelix(grn: string, h: number): boolean {
  if (grn.startsWith("n.") || grn.startsWith("c.")) return false;
  const parts = grn.split(".");
  return parts[0] === String(h);
}

function grnInLoop(grn: string, a: number, b: number): boolean {
  if (grn.startsWith("n.") || grn.startsWith("c.")) return false;
  const prefix = grn.split(".")[0];
  if (prefix.length !== 2) return false;
  const p = parseInt(prefix, 10);
  return p === a * 10 + b || p === b * 10 + a;
}

function sectionForGrn(grn: string): Section | undefined {
  return SECTIONS.find((s) => s.match(grn));
}

// ---------------------------------------------------------------------------
//  Helpers
// ---------------------------------------------------------------------------
function regionColor(grn: string): string {
  return sectionForGrn(grn)?.color || "#333";
}

function regionLabel(grn: string): string {
  return sectionForGrn(grn)?.label || "?";
}

// Amino acid property groups
const AA_GROUPS: Record<string, { bg: string; fg: string; label: string }> = {
  A: { bg: "#c8dbbe", fg: "#2d4a1e", label: "hydrophobic" },
  V: { bg: "#c8dbbe", fg: "#2d4a1e", label: "hydrophobic" },
  I: { bg: "#c8dbbe", fg: "#2d4a1e", label: "hydrophobic" },
  L: { bg: "#c8dbbe", fg: "#2d4a1e", label: "hydrophobic" },
  M: { bg: "#c8dbbe", fg: "#2d4a1e", label: "hydrophobic" },
  F: { bg: "#e8d5b7", fg: "#5c3d1a", label: "aromatic" },
  W: { bg: "#e8d5b7", fg: "#5c3d1a", label: "aromatic" },
  Y: { bg: "#e8d5b7", fg: "#5c3d1a", label: "aromatic" },
  S: { bg: "#d4e4f7", fg: "#1e3a5f", label: "polar" },
  T: { bg: "#d4e4f7", fg: "#1e3a5f", label: "polar" },
  N: { bg: "#d4e4f7", fg: "#1e3a5f", label: "polar" },
  Q: { bg: "#d4e4f7", fg: "#1e3a5f", label: "polar" },
  K: { bg: "#f7caca", fg: "#8b1a1a", label: "positive" },
  R: { bg: "#f7caca", fg: "#8b1a1a", label: "positive" },
  H: { bg: "#f7caca", fg: "#8b1a1a", label: "positive" },
  D: { bg: "#d0c4e8", fg: "#3b1f6e", label: "negative" },
  E: { bg: "#d0c4e8", fg: "#3b1f6e", label: "negative" },
  G: { bg: "#f0ece0", fg: "#6b6352", label: "special" },
  P: { bg: "#f0ece0", fg: "#6b6352", label: "special" },
  C: { bg: "#fff3bf", fg: "#7a6a00", label: "special" },
};

function residueStyle(val: string): { background: string; color: string } {
  if (!val || val === "-") return { background: "transparent", color: "#ddd" };
  const group = AA_GROUPS[val[0]];
  if (group) return { background: group.bg, color: group.fg };
  return { background: "rgba(0,0,0,0.03)", color: "#333" };
}

function regionBg(grn: string): string {
  const sec = sectionForGrn(grn);
  if (!sec) return "transparent";
  // Parse hex color to rgba with low alpha
  const c = sec.color;
  const r = parseInt(c.slice(1, 3), 16);
  const g = parseInt(c.slice(3, 5), 16);
  const b = parseInt(c.slice(5, 7), 16);
  return `rgba(${r},${g},${b},0.08)`;
}

const EXAMPLE_BR =
  "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLL" +
  "GYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVW" +
  "WAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDV" +
  "SAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATSD";

type Tab = "annotate" | "reference";

// ---------------------------------------------------------------------------
//  Section Legend (shared toggle bar)
// ---------------------------------------------------------------------------
function SectionLegend({
  hidden,
  onToggle,
  columns,
}: {
  hidden: Set<string>;
  onToggle: (key: string) => void;
  columns?: string[]; // if provided, only show sections that have columns
}) {
  const present = useMemo(() => {
    if (!columns) return new Set(SECTIONS.map((s) => s.key));
    const keys = new Set<string>();
    columns.forEach((col) => {
      const sec = sectionForGrn(col);
      if (sec) keys.add(sec.key);
    });
    return keys;
  }, [columns]);

  return (
    <div className="section-legend">
      <span className="section-legend-title">Sections:</span>
      {SECTIONS.filter((s) => present.has(s.key)).map((s) => {
        const off = hidden.has(s.key);
        return (
          <button
            key={s.key}
            className={`section-toggle ${off ? "off" : ""}`}
            style={{
              borderColor: s.color,
              background: off ? "transparent" : s.color,
              color: off ? s.color : "#fff",
            }}
            onClick={() => onToggle(s.key)}
            title={off ? `Show ${s.label}` : `Hide ${s.label}`}
          >
            {s.label}
          </button>
        );
      })}
      <button
        className="section-toggle-all"
        onClick={() => {
          // If any are hidden, show all. Otherwise hide loops+terminals.
          if (hidden.size > 0) {
            // clear all
            hidden.forEach((k) => onToggle(k));
          } else {
            // hide loops + N/C
            SECTIONS.filter(
              (s) => s.key.startsWith("L") || s.key === "N" || s.key === "C"
            ).forEach((s) => {
              if (!hidden.has(s.key)) onToggle(s.key);
            });
          }
        }}
      >
        {hidden.size > 0 ? "Show all" : "TM only"}
      </button>
    </div>
  );
}

// ---------------------------------------------------------------------------
//  App
// ---------------------------------------------------------------------------
function App() {
  const [tab, setTab] = useState<Tab>("annotate");
  const [hiddenSections, setHiddenSections] = useState<Set<string>>(new Set());

  // Annotate state
  const [sequence, setSequence] = useState("");
  const [singleResult, setSingleResult] = useState<AnnotateResponse | null>(null);
  const [batchResults, setBatchResults] = useState<AnnotateResponse[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [progress, setProgress] = useState("");
  const fileRef = useRef<HTMLInputElement>(null);

  // Reference state
  const [refTable, setRefTable] = useState<ReferenceTable | null>(null);
  const [refLoading, setRefLoading] = useState(false);
  const [refAddCount, setRefAddCount] = useState(0);

  useEffect(() => {
    setRefAddCount(0);
  }, [batchResults]);

  const toggleSection = useCallback((key: string) => {
    setHiddenSections((prev) => {
      const next = new Set(prev);
      if (next.has(key)) next.delete(key);
      else next.add(key);
      return next;
    });
  }, []);

  const isVisible = useCallback(
    (grn: string) => {
      if (hiddenSections.size === 0) return true;
      const sec = sectionForGrn(grn);
      return sec ? !hiddenSections.has(sec.key) : true;
    },
    [hiddenSections]
  );

  // ---------------------------------------------------------------------------
  //  Single sequence annotation
  // ---------------------------------------------------------------------------
  const handleSubmit = async () => {
    const clean = sequence.replace(/\s/g, "");
    if (clean.length < 10) {
      setError("Sequence must be at least 10 residues");
      return;
    }
    setLoading(true);
    setError(null);
    setSingleResult(null);
    setBatchResults([]);
    setProgress("Aligning...");

    try {
      const res = await fetch(`${API_URL}/annotate`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ sequence: clean, name: "query" }),
      });
      if (!res.ok) {
        const d = await res.json().catch(() => ({}));
        throw new Error(d.detail || `Server error ${res.status}`);
      }
      const data: AnnotateResponse = await res.json();
      setSingleResult(data);
      setBatchResults([data]);
    } catch (e: any) {
      setError(e.message || "Unknown error");
    } finally {
      setLoading(false);
      setProgress("");
    }
  };

  // ---------------------------------------------------------------------------
  //  FASTA upload
  // ---------------------------------------------------------------------------
  const handleFileUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    setLoading(true);
    setError(null);
    setSingleResult(null);
    setBatchResults([]);
    setProgress(`Uploading ${file.name}...`);

    const formData = new FormData();
    formData.append("file", file);

    try {
      const res = await fetch(`${API_URL}/annotate/fasta`, {
        method: "POST",
        body: formData,
      });
      if (!res.ok) {
        const d = await res.json().catch(() => ({}));
        throw new Error(d.detail || `Server error ${res.status}`);
      }
      const data = await res.json();
      const results: AnnotateResponse[] = data.results;
      setBatchResults(results);
      if (results.length === 1) {
        setSingleResult(results[0]);
      }
    } catch (e: any) {
      setError(e.message || "Unknown error");
    } finally {
      setLoading(false);
      setProgress("");
      if (fileRef.current) fileRef.current.value = "";
    }
  };

  // ---------------------------------------------------------------------------
  //  CSV download
  // ---------------------------------------------------------------------------
  const handleDownloadCSV = async () => {
    if (batchResults.length === 0) return;
    try {
      const res = await fetch(`${API_URL}/export/csv`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ results: batchResults }),
      });
      if (!res.ok) throw new Error("Export failed");
      const blob = await res.blob();
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "grn_annotations.csv";
      a.click();
      URL.revokeObjectURL(url);
    } catch (e: any) {
      setError(e.message);
    }
  };

  // ---------------------------------------------------------------------------
  //  Reference table
  // ---------------------------------------------------------------------------
  const fetchReference = useCallback(async (): Promise<ReferenceTable | null> => {
    try {
      const res = await fetch(`${API_URL}/reference`);
      if (!res.ok) throw new Error("Failed to load reference table");
      return (await res.json()) as ReferenceTable;
    } catch (e: any) {
      setError(e.message);
      return null;
    }
  }, []);

  const loadReference = useCallback(async () => {
    if (refTable) return;
    setRefLoading(true);
    const data = await fetchReference();
    if (data) setRefTable(data);
    setRefLoading(false);
  }, [refTable, fetchReference]);

  const addToReference = useCallback(async () => {
    if (batchResults.length === 0) return;
    let table = refTable;
    if (!table) {
      setRefLoading(true);
      table = await fetchReference();
      setRefLoading(false);
      if (!table) return;
    }

    // Annotation rows already share the canonical column space with the
    // reference table; project each onto the (filtered) reference column list.
    const usedNames = new Set(table.rows.map((r) => r.name));
    const newRows: ReferenceRow[] = [];
    for (const r of batchResults) {
      let name = r.name;
      let counter = 2;
      while (usedNames.has(name)) name = `${r.name}_${counter++}`;
      usedNames.add(name);
      const grnToVal: Record<string, string> = {};
      r.columns.forEach((c, i) => {
        grnToVal[c] = r.values[i];
      });
      const values = table.columns.map((c) => grnToVal[c] ?? "-");
      newRows.push({ name, values });
    }

    setRefTable({ ...table, rows: [...table.rows, ...newRows] });
    setRefAddCount(batchResults.length);
  }, [refTable, batchResults, fetchReference]);

  const handleExportReference = useCallback(() => {
    if (!refTable) return;
    const sortedCols = refTable.columns
      .slice()
      .sort((a, b) => grn2float(a) - grn2float(b));
    const colIdx = refTable.columns.map((c) => sortedCols.indexOf(c));
    const escape = (s: string) =>
      /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
    const lines: string[] = [];
    lines.push(["", ...sortedCols].map(escape).join(","));
    for (const row of refTable.rows) {
      const reordered = new Array(sortedCols.length).fill("-");
      row.values.forEach((v, i) => {
        const target = colIdx[i];
        if (target !== -1) reordered[target] = v;
      });
      lines.push([row.name, ...reordered].map(escape).join(","));
    }
    const blob = new Blob([lines.join("\n")], { type: "text/csv" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "reference_table.csv";
    a.click();
    URL.revokeObjectURL(url);
  }, [refTable]);

  useEffect(() => {
    if (tab === "reference") loadReference();
  }, [tab, loadReference]);

  const loadExample = () => {
    setSequence(EXAMPLE_BR);
    setSingleResult(null);
    setBatchResults([]);
    setError(null);
  };

  // ---------------------------------------------------------------------------
  //  Collect sorted GRN columns for batch table
  // ---------------------------------------------------------------------------
  const batchColumns = useMemo(() => {
    if (batchResults.length <= 1) return [];
    const allGrns = new Set<string>();
    batchResults.forEach((r) => r.annotations.forEach((a) => allGrns.add(a.grn)));
    return Array.from(allGrns).sort((a, b) => grn2float(a) - grn2float(b));
  }, [batchResults]);

  // Filtered columns for tables
  const filteredBatchCols = useMemo(
    () => batchColumns.filter(isVisible),
    [batchColumns, isVisible]
  );
  const filteredRefCols = useMemo(
    () =>
      refTable
        ? refTable.columns
            .filter(isVisible)
            .slice()
            .sort((a, b) => grn2float(a) - grn2float(b))
        : [],
    [refTable, isVisible]
  );

  return (
    <div className="app">
      <header className="header">
        <div className="header-top">
          <h1>GRN Annotator</h1>
          <nav className="tabs">
            <button
              className={`tab ${tab === "annotate" ? "active" : ""}`}
              onClick={() => setTab("annotate")}
            >
              Annotate
            </button>
            <button
              className={`tab ${tab === "reference" ? "active" : ""}`}
              onClick={() => setTab("reference")}
            >
              Reference Table
            </button>
          </nav>
        </div>
        <p className="subtitle">Generic Residue Numbering for Type I Opsins</p>
      </header>

      <main className="main">
        {/* ── ANNOTATE TAB ── */}
        {tab === "annotate" && (
          <>
            <section className="input-section">
              <div className="input-header">
                <label htmlFor="seq-input">Paste sequence</label>
                <div className="input-actions">
                  <button className="link-btn" onClick={loadExample}>
                    Load example (BR)
                  </button>
                  <span className="separator">|</span>
                  <label className="link-btn" htmlFor="fasta-upload">
                    Upload FASTA
                  </label>
                  <input
                    ref={fileRef}
                    id="fasta-upload"
                    type="file"
                    accept=".fasta,.fa,.faa,.txt"
                    style={{ display: "none" }}
                    onChange={handleFileUpload}
                  />
                </div>
              </div>
              <textarea
                id="seq-input"
                className="seq-input"
                placeholder="MLELLPTAVEGVSQAQITGRP..."
                value={sequence}
                onChange={(e) => setSequence(e.target.value)}
                rows={4}
                spellCheck={false}
              />
              <div className="actions">
                <button
                  className="btn-primary"
                  onClick={handleSubmit}
                  disabled={loading || sequence.replace(/\s/g, "").length < 10}
                >
                  {loading ? progress || "Annotating..." : "Annotate"}
                </button>
                {sequence && (
                  <span className="char-count">
                    {sequence.replace(/\s/g, "").length} residues
                  </span>
                )}
                {batchResults.length > 0 && (
                  <button className="btn-secondary" onClick={handleDownloadCSV}>
                    Download CSV
                  </button>
                )}
                {batchResults.length > 0 && (
                  <button className="btn-secondary" onClick={addToReference}>
                    Add to reference table
                  </button>
                )}
                {refAddCount > 0 && (
                  <span className="muted">
                    Added {refAddCount} sequence{refAddCount === 1 ? "" : "s"}.{" "}
                    <button className="link-btn" onClick={() => setTab("reference")}>
                      View
                    </button>
                  </span>
                )}
              </div>
              {error && <div className="error">{error}</div>}
            </section>

            {/* Single-sequence horizontal result */}
            {singleResult && batchResults.length === 1 && (
              <SingleResult
                result={singleResult}
                hidden={hiddenSections}
                onToggle={toggleSection}
                isVisible={isVisible}
              />
            )}

            {/* Multi-sequence table result */}
            {batchResults.length > 1 && (
              <section className="result-section">
                <div className="result-meta">
                  <span>{batchResults.length} sequences annotated</span>
                  <button className="btn-secondary" onClick={handleDownloadCSV}>
                    Download CSV
                  </button>
                  <button className="btn-secondary" onClick={addToReference}>
                    Add to reference table
                  </button>
                </div>
                <SectionLegend
                  hidden={hiddenSections}
                  onToggle={toggleSection}
                  columns={batchColumns}
                />
                <GRNTable
                  columns={filteredBatchCols}
                  rows={batchResults.map((r) => ({
                    name: r.name,
                    grn_map: Object.fromEntries(
                      r.annotations.map((a) => [a.grn, `${a.residue}${a.position}`])
                    ),
                  }))}
                />
              </section>
            )}
          </>
        )}

        {/* ── REFERENCE TAB ── */}
        {tab === "reference" && (
          <section className="result-section">
            {refLoading && <p className="loading-text">Loading reference table...</p>}
            {refTable && (
              <>
                <div className="result-meta">
                  <span>
                    {refTable.rows.length} sequences &times;{" "}
                    {filteredRefCols.length} / {refTable.columns.length} GRN positions
                  </span>
                  <button className="btn-secondary" onClick={handleExportReference}>
                    Download CSV
                  </button>
                </div>
                <SectionLegend
                  hidden={hiddenSections}
                  onToggle={toggleSection}
                  columns={refTable.columns}
                />
                <GRNTable
                  columns={filteredRefCols}
                  rows={refTable.rows.map((r) => ({
                    name: r.name,
                    grn_map: Object.fromEntries(
                      r.values.map((v, i) => [refTable.columns[i], v])
                    ),
                  }))}
                />
              </>
            )}
          </section>
        )}
      </main>

      <footer className="footer">
        <p>GRN assignment based on the Protos framework</p>
      </footer>
    </div>
  );
}

// ---------------------------------------------------------------------------
//  Single-sequence horizontal row
// ---------------------------------------------------------------------------
function SingleResult({
  result,
  hidden,
  onToggle,
  isVisible,
}: {
  result: AnnotateResponse;
  hidden: Set<string>;
  onToggle: (key: string) => void;
  isVisible: (grn: string) => boolean;
}) {
  const filtered = useMemo(
    () => result.annotations.filter((a) => isVisible(a.grn)),
    [result.annotations, isVisible]
  );

  return (
    <section className="result-section">
      <div className="result-meta">
        <span>
          {filtered.length} / {result.sequence_length} shown
        </span>
        {hidden.size > 0 && (
          <span className="muted">
            ({result.annotations.length - filtered.length} hidden)
          </span>
        )}
        {result.missing_count > 0 && (
          <span className="warn">{result.missing_count} missing</span>
        )}
      </div>

      <SectionLegend
        hidden={hidden}
        onToggle={onToggle}
        columns={result.annotations.map((a) => a.grn)}
      />

      <div className="grn-scroll-container">
        <div className="grn-row">
          {filtered.map((a, i) => {
            const color = regionColor(a.grn);
            const showLabel =
              i === 0 ||
              regionLabel(a.grn) !== regionLabel(filtered[i - 1].grn);
            return (
              <div
                key={`${a.position}-${a.grn}`}
                className="grn-cell"
                style={{ borderBottomColor: color }}
                title={`${a.residue}${a.position} = ${a.grn} (${regionLabel(a.grn)})`}
              >
                {showLabel ? (
                  <div className="grn-region" style={{ color }}>
                    {regionLabel(a.grn)}
                  </div>
                ) : (
                  <div className="grn-region">&nbsp;</div>
                )}
                <div className="grn-number">{a.grn}</div>
                <div className="grn-residue" style={{ color }}>
                  {a.residue}
                  <span className="grn-pos">{a.position}</span>
                </div>
              </div>
            );
          })}
        </div>
      </div>
    </section>
  );
}

// ---------------------------------------------------------------------------
//  Shared scrollable GRN table (for reference + batch results)
// ---------------------------------------------------------------------------
function GRNTable({
  columns,
  rows,
}: {
  columns: string[];
  rows: { name: string; grn_map: Record<string, string> }[];
}) {
  return (
    <div className="grn-table-scroll">
      <table className="grn-table">
        <thead>
          <tr>
            <th className="grn-table-name-col sticky-col">Name</th>
            {columns.map((col) => (
              <th
                key={col}
                className="grn-table-header"
                style={{
                  borderBottomColor: regionColor(col),
                  background: regionBg(col),
                }}
              >
                <span className="grn-table-region" style={{ color: regionColor(col) }}>
                  {regionLabel(col)}
                </span>
                <span className="grn-table-colname">{col}</span>
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows.map((row) => (
            <tr key={row.name}>
              <td className="grn-table-name-col sticky-col" title={row.name}>
                {row.name}
              </td>
              {columns.map((col) => {
                const val = row.grn_map[col] || "-";
                const rs = residueStyle(val);
                return (
                  <td
                    key={col}
                    className={`grn-table-cell ${val === "-" ? "gap" : ""}`}
                    style={{ background: rs.background, color: rs.color }}
                    title={
                      val !== "-"
                        ? `${val} (${AA_GROUPS[val[0]]?.label || "unknown"})`
                        : ""
                    }
                  >
                    {val === "-" ? "" : val}
                  </td>
                );
              })}
            </tr>
          ))}
        </tbody>
      </table>

      <div className="aa-legend">
        <span className="aa-legend-title">Residue type:</span>
        {[
          { label: "Hydrophobic", bg: "#c8dbbe", fg: "#2d4a1e" },
          { label: "Aromatic", bg: "#e8d5b7", fg: "#5c3d1a" },
          { label: "Polar", bg: "#d4e4f7", fg: "#1e3a5f" },
          { label: "Positive", bg: "#f7caca", fg: "#8b1a1a" },
          { label: "Negative", bg: "#d0c4e8", fg: "#3b1f6e" },
          { label: "Special", bg: "#f0ece0", fg: "#6b6352" },
        ].map((g) => (
          <span
            key={g.label}
            className="aa-legend-swatch"
            style={{ background: g.bg, color: g.fg }}
          >
            {g.label}
          </span>
        ))}
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
//  GRN float for sorting
//
//  Produces a key such that columns sort in biological N-to-C order:
//    N-tail (n.k)  <  TM1  <  loop 1-2  <  TM2  <  loop 2-3  <  …  <  TM7  <  C-tail (c.k)
//
//  Canonical helix positions "H.pp" keep their natural numeric value.
//  Insertions within a helix "H.ppp" also sort naturally (e.g. 3.521 between 3.52 and 3.53).
//  Inter-helical loops use a two-digit prefix "HxHy" where Hx is the helix the residue
//  is nearer to. For a loop between helix a and helix a+1 (a<b=a+1) we need to map:
//    - ab.pos residues (near Hx=a, ordered ascending with pos toward loop midpoint)
//    - ba.pos residues (near Hx=b, where smaller pos = closer to b = later in sequence)
//  into the open interval (a + last_TM_value, a+1 + first_TM_value) so they sit
//  between helix a and helix b. We pack them into (a+0.9, b) with a small multiplier:
//    ab.pos -> a + 0.90 + pos * 1e-4   (ascending with pos, closer to a first)
//    ba.pos -> b          - pos * 1e-4  (descending with pos, closer to b last)
// ---------------------------------------------------------------------------
function grn2float(grn: string): number {
  if (grn.startsWith("n.")) return -parseInt(grn.slice(2), 10);
  if (grn.startsWith("c.")) return 100 + parseInt(grn.slice(2), 10);
  const dot = grn.indexOf(".");
  if (dot <= 0) return parseFloat(grn) || 0;
  const prefix = grn.slice(0, dot);
  const pos = parseInt(grn.slice(dot + 1), 10) || 0;
  if (prefix.length === 2 && /^\d{2}$/.test(prefix)) {
    const hx = parseInt(prefix[0], 10);
    const hy = parseInt(prefix[1], 10);
    const a = Math.min(hx, hy);
    const b = Math.max(hx, hy);
    if (hx === a) {
      return a + 0.9 + pos * 1e-4;
    }
    return b - pos * 1e-4;
  }
  return parseFloat(grn) || 0;
}

export default App;
