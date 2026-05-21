import React, { useState, useRef, useEffect, useCallback, useMemo } from "react";
import "./App.css";

// 3Dmol (public/3Dmol-min.js) is lazy-loaded the first time the Structure tab
// opens — it attaches to `window.$3Dmol`. This keeps initial page load light
// (≈ 526 KB script + parse cost) for users who never visit the Structure tab.
declare global {
  interface Window {
    $3Dmol: any;
  }
}
let load3DmolPromise: Promise<void> | null = null;
function load3Dmol(): Promise<void> {
  if (load3DmolPromise) return load3DmolPromise;
  load3DmolPromise = new Promise<void>((resolve, reject) => {
    if (typeof window === "undefined") {
      reject(new Error("window not available"));
      return;
    }
    if (window.$3Dmol) {
      resolve();
      return;
    }
    const script = document.createElement("script");
    script.src = "/3Dmol-min.js";
    script.async = true;
    script.onload = () => resolve();
    script.onerror = () => {
      load3DmolPromise = null; // allow retry on next mount
      reject(new Error("Failed to load 3Dmol"));
    };
    document.head.appendChild(script);
  });
  return load3DmolPromise;
}

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

interface CatalogEntry {
  name: string;
  display_name: string | null;
  species: string | null;
  domain: string | null;
  function: string | null;
  function_detail: string | null;
  pdb_id: string | null;
  method: string | null;
  resolution: number | null;
  reference: string | null;
  reference_year: number | null;
  length: number | null;
  uniprot_id: string | null;
  sequence: string | null;
}

interface ReferenceRow {
  name: string;
  values: string[];
  metadata: CatalogEntry | null;
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

// Display label for a structure id: experimental PDB ids pass through; predicted
// model ids replace `_model_0` with `_boltz_model` to make their provenance clear.
function displayStructureId(structureId: string): string {
  return structureId.replace(/_model_0$/, "_boltz_model");
}

function rowTooltip(r: ReferenceRow): string {
  const m = r.metadata;
  if (!m) return r.name;
  const parts = [m.name];
  if (m.species) parts.push(m.species);
  if (m.function) parts.push(m.function);
  if (m.method && m.resolution) parts.push(`${m.method} ${m.resolution}Å`);
  else if (m.method) parts.push(m.method);
  if (m.reference) parts.push(m.reference);
  return parts.join(" · ");
}

// Display name for a reference row: prefer the catalog descriptive name when
// the join succeeded (falls back to the short catalog name, then to the raw
// row id with the `_model_N` suffix stripped).
function refRowDisplayName(r: ReferenceRow): string {
  if (r.metadata?.display_name) return r.metadata.display_name;
  if (r.metadata?.name) return r.metadata.name;
  let base = r.name;
  for (const suffix of ["_model_0", "_model_1", "_model_2", "_model_3"]) {
    if (base.endsWith(suffix)) {
      base = base.slice(0, -suffix.length);
      break;
    }
  }
  return base.replace(/_+$/, "");
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

type Tab = "annotate" | "reference" | "catalog" | "grn" | "structure";

interface StructureManifestEntry {
  structure_id: string;
  structure_type: string;
  n_atoms: number;
  n_grn_residues: number;
  name: string;
  display_name: string | null;
  species: string | null;
  function: string | null;
  pdb_id: string | null;
}

interface StructureGRNResidue {
  resi: number;
  res1: string;
}

interface StructureResponse {
  structure_id: string;
  structure_type: string;
  pdb: string;
  residue_grn: Record<string, string>;
  grn_residue: Record<string, StructureGRNResidue>;
  metadata: CatalogEntry | null;
}

interface GRNPositionEntry {
  grn: string;
  residue: string | null;
  label: string;
  description: string;
}

interface GRNCategoryEntry {
  key: string;
  name: string;
  color: string;
  summary: string;
  positions: GRNPositionEntry[];
  relevant_for: string[];
}

const FUNCTION_TONES: Record<string, string> = {
  "Proton Pump": "#c81e1e",
  "Sodium Pump": "#d97706",
  "Chloride Pump": "#0e7490",
  "Anion Channel": "#1d4ed8",
  "Cation Channel": "#7c3aed",
  "Sensor / Regulatory": "#057a55",
  "Unknown": "#8a8780",
};

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
  const [tab, setTab] = useState<Tab>("catalog");
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

  // Catalog state
  const [catalog, setCatalog] = useState<CatalogEntry[] | null>(null);
  const [catalogLoading, setCatalogLoading] = useState(false);
  const [catalogQuery, setCatalogQuery] = useState("");
  const [catalogFunction, setCatalogFunction] = useState<string | null>(null);

  // GRN function-map state
  const [grnFunctions, setGrnFunctions] = useState<GRNCategoryEntry[] | null>(null);
  const [grnLoading, setGrnLoading] = useState(false);
  const [grnActiveKey, setGrnActiveKey] = useState<string | null>(null);
  const [grnQuery, setGrnQuery] = useState("");

  // About / contact modal
  const [aboutOpen, setAboutOpen] = useState(false);
  useEffect(() => {
    if (!aboutOpen) return;
    const onKey = (e: KeyboardEvent) => {
      if (e.key === "Escape") setAboutOpen(false);
    };
    document.addEventListener("keydown", onKey);
    return () => document.removeEventListener("keydown", onKey);
  }, [aboutOpen]);

  // Structure tab state
  const [structureManifest, setStructureManifest] = useState<StructureManifestEntry[] | null>(null);
  const [structureManifestLoading, setStructureManifestLoading] = useState(false);
  const [structureData, setStructureData] = useState<Record<string, StructureResponse>>({});
  const [selectedStructures, setSelectedStructures] = useState<string[]>(["1c3w"]);
  const [structureLoadingIds, setStructureLoadingIds] = useState<Set<string>>(new Set());
  const [selectedGrns, setSelectedGrns] = useState<Set<string>>(new Set());

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
    // Empty textarea → fall back to the placeholder example (HsBR).
    const usingExample = clean.length === 0;
    const seqToUse = usingExample ? EXAMPLE_BR.replace(/\s/g, "") : clean;
    if (seqToUse.length < 10) {
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
        body: JSON.stringify({
          sequence: seqToUse,
          name: usingExample ? "HsBR" : "query",
        }),
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
      r.annotations.forEach((a) => {
        grnToVal[a.grn] = `${a.residue}${a.position}`;
      });
      const values = table.columns.map((c) => grnToVal[c] ?? "-");
      newRows.push({ name, values, metadata: null });
    }

    setRefTable({ ...table, rows: [...table.rows, ...newRows] });
    setRefAddCount(batchResults.length);
    setTab("reference");
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

  // ---------------------------------------------------------------------------
  //  Catalog (Zenodo)
  // ---------------------------------------------------------------------------
  const loadCatalog = useCallback(async () => {
    if (catalog) return;
    setCatalogLoading(true);
    try {
      const res = await fetch(`${API_URL}/catalog`);
      if (!res.ok) throw new Error("Failed to load catalog");
      const data = await res.json();
      setCatalog(data.entries as CatalogEntry[]);
    } catch (e: any) {
      setError(e.message);
    } finally {
      setCatalogLoading(false);
    }
  }, [catalog]);

  useEffect(() => {
    if (tab === "catalog") loadCatalog();
  }, [tab, loadCatalog]);

  // ---------------------------------------------------------------------------
  //  GRN function map
  // ---------------------------------------------------------------------------
  const loadGrnFunctions = useCallback(async () => {
    if (grnFunctions) return;
    setGrnLoading(true);
    try {
      const res = await fetch(`${API_URL}/grn-functions`);
      if (!res.ok) throw new Error("Failed to load GRN function map");
      const data = await res.json();
      setGrnFunctions(data.categories as GRNCategoryEntry[]);
    } catch (e: any) {
      setError(e.message);
    } finally {
      setGrnLoading(false);
    }
  }, [grnFunctions]);

  useEffect(() => {
    if (tab === "grn" || tab === "structure") loadGrnFunctions();
  }, [tab, loadGrnFunctions]);

  // ---------------------------------------------------------------------------
  //  Structure tab — manifest + per-structure fetch
  // ---------------------------------------------------------------------------
  const loadStructureManifest = useCallback(async () => {
    if (structureManifest) return;
    setStructureManifestLoading(true);
    try {
      const res = await fetch(`${API_URL}/structures`);
      if (!res.ok) throw new Error("Failed to load structure manifest");
      const data = await res.json();
      setStructureManifest(data.entries as StructureManifestEntry[]);
    } catch (e: any) {
      setError(e.message);
    } finally {
      setStructureManifestLoading(false);
    }
  }, [structureManifest]);

  const loadStructure = useCallback(
    async (structure_id: string) => {
      if (structureData[structure_id]) return;
      setStructureLoadingIds((s) => new Set(s).add(structure_id));
      try {
        const res = await fetch(`${API_URL}/structures/${encodeURIComponent(structure_id)}`);
        if (!res.ok) throw new Error(`Failed to load ${structure_id}`);
        const data: StructureResponse = await res.json();
        setStructureData((prev) => ({ ...prev, [structure_id]: data }));
      } catch (e: any) {
        setError(e.message);
      } finally {
        setStructureLoadingIds((s) => {
          const next = new Set(s);
          next.delete(structure_id);
          return next;
        });
      }
    },
    [structureData]
  );

  useEffect(() => {
    if (tab === "structure" || tab === "catalog" || tab === "reference") loadStructureManifest();
  }, [tab, loadStructureManifest]);

  // Navigate from Catalog → Structure tab. Resolves a catalog entry to the
  // matching structure_id (lowercased PDB for experimental, "<name>_model_0"
  // for predicted). Returns null if no structure is bundled for the entry.
  const resolveStructureId = useCallback(
    (entry: CatalogEntry): string | null => {
      if (!structureManifest) return null;
      const ids = new Set(structureManifest.map((m) => m.structure_id));
      if (entry.pdb_id) {
        const sid = entry.pdb_id.toLowerCase();
        if (ids.has(sid)) return sid;
      }
      const predicted = `${entry.name}_model_0`;
      if (ids.has(predicted)) return predicted;
      // some names carry trailing underscores in the manifest (e.g. ZipACR__model_0)
      const predicted2 = `${entry.name}__model_0`;
      if (ids.has(predicted2)) return predicted2;
      return null;
    },
    [structureManifest]
  );

  const openStructureFromCatalog = useCallback(
    (entry: CatalogEntry) => {
      const sid = resolveStructureId(entry);
      if (!sid) {
        setError(`No bundled structure for ${entry.name}.`);
        return;
      }
      setSelectedStructures([sid]);
      setSelectedGrns(new Set());
      setTab("structure");
    },
    [resolveStructureId]
  );

  // Reference rows often lowercase predicted-model names; resolve them back to
  // the manifest's original-case structure_id via a case-insensitive lookup.
  const resolveStructureIdByRaw = useCallback(
    (rawName: string): string | null => {
      if (!structureManifest) return null;
      const lower = rawName.toLowerCase();
      for (const m of structureManifest) {
        if (m.structure_id.toLowerCase() === lower) return m.structure_id;
      }
      return null;
    },
    [structureManifest]
  );

  const openStructureFromRefRow = useCallback(
    (rawName: string) => {
      const sid = resolveStructureIdByRaw(rawName);
      if (!sid) {
        setError(`No bundled structure for ${rawName}.`);
        return;
      }
      setSelectedStructures([sid]);
      setSelectedGrns(new Set());
      setTab("structure");
    },
    [resolveStructureIdByRaw]
  );

  useEffect(() => {
    if (tab !== "structure") return;
    selectedStructures.forEach((id) => {
      if (!structureData[id] && !structureLoadingIds.has(id)) {
        loadStructure(id);
      }
    });
  }, [tab, selectedStructures, structureData, structureLoadingIds, loadStructure]);

  const toggleStructure = useCallback((id: string) => {
    setSelectedStructures((prev) =>
      prev.includes(id) ? prev.filter((x) => x !== id) : [...prev, id]
    );
  }, []);

  const toggleGrn = useCallback((grn: string, additive: boolean) => {
    setSelectedGrns((prev) => {
      const next = new Set(additive ? prev : []);
      if (prev.has(grn) && additive) next.delete(grn);
      else next.add(grn);
      return next;
    });
  }, []);

  const clearGrnSelection = useCallback(() => setSelectedGrns(new Set()), []);

  const grnVisibleCategories = useMemo(() => {
    if (!grnFunctions) return [];
    const q = grnQuery.trim().toLowerCase();
    let list = grnActiveKey
      ? grnFunctions.filter((c) => c.key === grnActiveKey)
      : grnFunctions;
    if (q) {
      list = list
        .map((c) => ({
          ...c,
          positions: c.positions.filter(
            (p) =>
              p.grn.toLowerCase().includes(q) ||
              p.label.toLowerCase().includes(q) ||
              p.description.toLowerCase().includes(q) ||
              (p.residue || "").toLowerCase().includes(q)
          ),
        }))
        .filter(
          (c) =>
            c.positions.length > 0 ||
            c.name.toLowerCase().includes(q) ||
            c.summary.toLowerCase().includes(q)
        );
    }
    return list;
  }, [grnFunctions, grnActiveKey, grnQuery]);

  const filteredCatalog = useMemo(() => {
    if (!catalog) return [];
    const q = catalogQuery.trim().toLowerCase();
    return catalog.filter((e) => {
      if (catalogFunction && e.function !== catalogFunction) return false;
      if (!q) return true;
      return (
        e.name.toLowerCase().includes(q) ||
        (e.display_name || "").toLowerCase().includes(q) ||
        (e.species || "").toLowerCase().includes(q) ||
        (e.pdb_id || "").toLowerCase().includes(q) ||
        (e.uniprot_id || "").toLowerCase().includes(q) ||
        (e.function || "").toLowerCase().includes(q)
      );
    });
  }, [catalog, catalogQuery, catalogFunction]);

  const catalogFunctions = useMemo(() => {
    if (!catalog) return [];
    const fns = new Set<string>();
    catalog.forEach((e) => e.function && fns.add(e.function));
    return Array.from(fns).sort();
  }, [catalog]);

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
          <nav className="tabs">
            <button
              className={`tab ${tab === "catalog" ? "active" : ""}`}
              onClick={() => setTab("catalog")}
            >
              Catalog
            </button>
            <button
              className={`tab ${tab === "reference" ? "active" : ""}`}
              onClick={() => setTab("reference")}
            >
              Reference GRNs
            </button>
            <button
              className={`tab ${tab === "grn" ? "active" : ""}`}
              onClick={() => setTab("grn")}
            >
              GRN Overview
            </button>
            <button
              className={`tab ${tab === "structure" ? "active" : ""}`}
              onClick={() => setTab("structure")}
            >
              Structure
            </button>
            <button
              className={`tab ${tab === "annotate" ? "active" : ""}`}
              onClick={() => setTab("annotate")}
            >
              Annotate
            </button>
          </nav>
        </div>
      </header>

      <main className="main">
        {/* ── ANNOTATE TAB ── */}
        {tab === "annotate" && (
          <>
            <section className="input-section">
              <div className="input-actions">
                <label className="link-btn" htmlFor="fasta-upload">
                  Upload FASTA
                </label>
                <span className="info-wrapper">
                  <span
                    className="info-sup"
                    tabIndex={0}
                    role="button"
                    aria-label="FASTA upload help"
                  >
                    i
                  </span>
                  <div className="info-popover" role="tooltip">
                    <p className="info-popover-text">
                      A FASTA file may contain <strong>one or more</strong>{" "}
                      sequences. Each is annotated independently against the MO
                      numbering reference. Accepts{" "}
                      <code>.fasta / .fa / .faa / .txt</code> (UTF-8). Up to
                      500 sequences per file.
                    </p>
                    <pre className="info-popover-code">
{`>sp|P02945|BACH_HALSA Bacteriorhodopsin OS=Halobacterium salinarum
MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKK...
>sp|P71411|BOPS_HALSO Bacteriorhodopsin-like protein OS=Haloarcula sp.
MLELQPTIAEHSIELQAEVIGRPAWIWLALGTLMGLGTLYFLVKGAGVADPQTK...
>sp|Q5NR01|CHOP1_CHLRE Channelrhodopsin-1 OS=Chlamydomonas reinhardtii
MSRRPWLLALALAVALAAGSAGASTGSDATVPVATQDGPDYVFHRAHERMLFQT...`}
                    </pre>
                  </div>
                </span>
                <input
                  ref={fileRef}
                  id="fasta-upload"
                  type="file"
                  accept=".fasta,.fa,.faa,.txt"
                  style={{ display: "none" }}
                  onChange={handleFileUpload}
                />
              </div>
              <textarea
                id="seq-input"
                className="seq-input"
                placeholder={EXAMPLE_BR}
                value={sequence}
                onChange={(e) => setSequence(e.target.value)}
                rows={4}
                spellCheck={false}
              />
              <div className="actions">
                <button
                  className="btn-primary"
                  onClick={handleSubmit}
                  disabled={loading}
                  title={
                    sequence.replace(/\s/g, "").length === 0
                      ? "No sequence entered — clicking annotates the displayed example (HsBR)"
                      : undefined
                  }
                >
                  {loading
                    ? progress || "Annotating..."
                    : sequence.replace(/\s/g, "").length === 0
                    ? "Annotate example"
                    : "Annotate"}
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
                    Add to reference
                  </button>
                )}
                {refAddCount > 0 && (
                  <span className="muted">
                    +{refAddCount}
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
                    Add to reference
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
                  rows={refTable.rows.map((r) => {
                    const hasStructure = resolveStructureIdByRaw(r.name) !== null;
                    return {
                      name: refRowDisplayName(r),
                      grn_map: Object.fromEntries(
                        r.values.map((v, i) => [refTable.columns[i], v])
                      ),
                      nameTitle: rowTooltip(r),
                      onNameClick: hasStructure
                        ? () => openStructureFromRefRow(r.name)
                        : undefined,
                    };
                  })}
                />
              </>
            )}
          </section>
        )}

        {/* ── CATALOG TAB ── */}
        {tab === "catalog" && (
          <section className="result-section">
            {catalogLoading && <p className="loading-text">Loading catalog…</p>}
            {catalog && (
              <CatalogView
                entries={filteredCatalog}
                total={catalog.length}
                query={catalogQuery}
                onQuery={setCatalogQuery}
                fn={catalogFunction}
                onFn={setCatalogFunction}
                functions={catalogFunctions}
                onOpenStructure={openStructureFromCatalog}
                resolveStructureId={resolveStructureId}
              />
            )}
          </section>
        )}

        {/* ── GRN FUNCTION-MAP TAB ── */}
        {tab === "grn" && (
          <section className="result-section">
            {grnLoading && <p className="loading-text">Loading GRN function map…</p>}
            {grnFunctions && (
              <GRNFunctionsView
                categories={grnFunctions}
                visible={grnVisibleCategories}
                activeKey={grnActiveKey}
                onActiveKey={setGrnActiveKey}
                query={grnQuery}
                onQuery={setGrnQuery}
              />
            )}
          </section>
        )}

        {/* ── STRUCTURE TAB ── */}
        {tab === "structure" && (
          <section className="result-section">
            {structureManifestLoading && (
              <p className="loading-text">Loading structure index…</p>
            )}
            {structureManifest && (
              <StructureView
                manifest={structureManifest}
                selected={selectedStructures}
                onToggle={toggleStructure}
                data={structureData}
                loading={structureLoadingIds}
                selectedGrns={selectedGrns}
                onToggleGrn={toggleGrn}
                onClearGrns={clearGrnSelection}
                onSetSelectedGrns={setSelectedGrns}
                grnFunctions={grnFunctions}
              />
            )}
          </section>
        )}
      </main>

      <footer className="footer">
        <p>
          Driven by the{" "}
          <a
            href="https://github.com/flurinh/protos"
            target="_blank"
            rel="noopener noreferrer"
            className="footer-link"
          >
            ProtOS framework
          </a>
          <span className="footer-sep">·</span>
          <button
            type="button"
            className="footer-link footer-btn"
            onClick={() => setAboutOpen(true)}
          >
            About / Contact
          </button>
          <span className="footer-sep">·</span>
          <span className="footer-cite">Manuscript under review</span>
          <span className="footer-sep">·</span>
          <a
            href="https://creativecommons.org/licenses/by/4.0/"
            target="_blank"
            rel="noopener noreferrer"
            className="footer-link"
          >
            CC-BY-4.0
          </a>
        </p>
      </footer>

      {aboutOpen && (
        <AboutModal onClose={() => setAboutOpen(false)} />
      )}
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
                  {a.residue}{a.position}
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
  rows: {
    name: string;
    grn_map: Record<string, string>;
    nameTitle?: string;
    onNameClick?: () => void;
  }[];
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
              <td className="grn-table-name-col sticky-col" title={row.nameTitle || row.name}>
                {row.onNameClick ? (
                  <button
                    type="button"
                    className="grn-table-name-link"
                    onClick={row.onNameClick}
                    title={`Open ${row.name} in the structure viewer`}
                  >
                    {row.name}
                  </button>
                ) : (
                  row.name
                )}
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
//  Catalog view
// ---------------------------------------------------------------------------
function CatalogView({
  entries,
  total,
  query,
  onQuery,
  fn,
  onFn,
  functions,
  onOpenStructure,
  resolveStructureId,
}: {
  entries: CatalogEntry[];
  total: number;
  query: string;
  onQuery: (q: string) => void;
  fn: string | null;
  onFn: (f: string | null) => void;
  functions: string[];
  onOpenStructure: (entry: CatalogEntry) => void;
  resolveStructureId: (entry: CatalogEntry) => string | null;
}) {
  return (
    <>
      <div className="catalog-controls">
        <input
          type="search"
          className="catalog-search"
          placeholder="Search name, species, PDB, UniProt…"
          value={query}
          onChange={(e) => onQuery(e.target.value)}
        />
        <div className="catalog-fn-filter">
          <button
            className={`catalog-fn-chip ${fn === null ? "active" : ""}`}
            onClick={() => onFn(null)}
          >
            All
          </button>
          {functions.map((f) => {
            const c = FUNCTION_TONES[f] || "#8a8780";
            const active = fn === f;
            return (
              <button
                key={f}
                className={`catalog-fn-chip ${active ? "active" : ""}`}
                style={{
                  borderColor: c,
                  background: active ? c : "transparent",
                  color: active ? "#fff" : c,
                }}
                onClick={() => onFn(active ? null : f)}
              >
                {f}
              </button>
            );
          })}
        </div>
        <span className="catalog-count">
          {entries.length} / {total}
        </span>
      </div>

      <div className="catalog-table-scroll">
        <table className="catalog-table">
          <thead>
            <tr>
              <th>Name</th>
              <th>Species</th>
              <th>Function</th>
              <th className="num">Len</th>
              <th>PDB</th>
              <th>Method</th>
              <th className="num">Res (Å)</th>
              <th>UniProt</th>
              <th>Reference</th>
            </tr>
          </thead>
          <tbody>
            {entries.map((e) => {
              const c = e.function ? FUNCTION_TONES[e.function] || "#8a8780" : "#8a8780";
              const hasStructure = resolveStructureId(e) !== null;
              return (
                <tr key={`${e.name}-${e.pdb_id || "x"}`}>
                  <td className="catalog-name" title={e.name}>
                    {hasStructure ? (
                      <button
                        type="button"
                        className="catalog-name-link"
                        onClick={() => onOpenStructure(e)}
                        title={`Open ${e.display_name || e.name} in the structure viewer`}
                      >
                        {e.display_name || e.name}
                      </button>
                    ) : (
                      e.display_name || e.name
                    )}
                  </td>
                  <td className="catalog-species">{e.species || "—"}</td>
                  <td>
                    {e.function && (
                      <span className="catalog-fn-pill" style={{ borderColor: c, color: c }}>
                        {e.function}
                      </span>
                    )}
                  </td>
                  <td className="num">{e.length ?? "—"}</td>
                  <td className="catalog-pdb">
                    {e.pdb_id ? (
                      <a
                        href={`https://www.rcsb.org/structure/${e.pdb_id}`}
                        target="_blank"
                        rel="noopener noreferrer"
                      >
                        {e.pdb_id}
                      </a>
                    ) : (
                      "—"
                    )}
                  </td>
                  <td>{e.method || "—"}</td>
                  <td className="num">{e.resolution ?? "—"}</td>
                  <td className="catalog-uniprot">
                    {e.uniprot_id ? (
                      <a
                        href={`https://www.uniprot.org/uniprotkb/${e.uniprot_id}`}
                        target="_blank"
                        rel="noopener noreferrer"
                      >
                        {e.uniprot_id}
                      </a>
                    ) : (
                      "—"
                    )}
                  </td>
                  <td className="catalog-ref" title={e.reference || ""}>
                    {e.reference ? (
                      <a
                        href={`https://scholar.google.com/scholar?q=${encodeURIComponent(e.reference)}`}
                        target="_blank"
                        rel="noopener noreferrer"
                        className="catalog-ref-link"
                      >
                        {e.reference}
                      </a>
                    ) : (
                      "—"
                    )}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </>
  );
}

// ---------------------------------------------------------------------------
//  About / Contact modal
// ---------------------------------------------------------------------------
function AboutModal({ onClose }: { onClose: () => void }) {
  const contacts = [
    {
      name: "Xavier Deupi",
      email: "xavier.deupi@psi.ch",
      topic: "Analytical questions",
    },
    {
      name: "Hideaki E. Kato",
      email: "c-hekato@g.ecc.u-tokyo.ac.jp",
      topic: "Functional properties, GRN assignments",
    },
    {
      name: "Flurin Hidber",
      email: "flurin.hidber@psi.ch",
      topic: "Bugs, feature requests",
    },
  ];
  return (
    <div
      className="modal-backdrop"
      onClick={(e) => {
        if (e.target === e.currentTarget) onClose();
      }}
    >
      <div className="modal-card" role="dialog" aria-modal="true" aria-label="About">
        <header className="modal-header">
          <span className="modal-title">About this resource</span>
          <button
            type="button"
            className="modal-close"
            onClick={onClose}
            aria-label="Close"
          >
            ×
          </button>
        </header>
        <section className="modal-section">
          <p className="modal-cite-label">How MO numbering works</p>
          <p className="modal-cite-body">
            The MO (Microbial Opsin) numbering system is a structure-guided
            framework anchored to the retinal-binding pocket of microbial
            rhodopsins. Each transmembrane residue receives a generic identifier
            of the form <code className="modal-inline-code">helix.position</code>
            {" "}— for example, <code className="modal-inline-code">3.45</code>{" "}
            denotes a residue on TM3.
          </p>
          <ul className="modal-bullet-list">
            <li>
              <strong>Anchor positions (X.50)</strong> are the closest residues
              to retinal in each helix across 129 reference structures, with
              sequence conservation as a tiebreaker. Numbers decrease toward the
              N-terminus and increase toward the C-terminus.
            </li>
            <li>
              <strong>The Schiff-base lysine</strong> on TM7 — covalently bound
              to all-trans retinal — defines{" "}
              <code className="modal-inline-code">7.50</code> and grounds the
              entire system. <code className="modal-inline-code">D85<sup>3.45</sup></code>{" "}
              in HsBR means: Asp85 in HsBR, at generic position 3.45.
            </li>
            <li>
              <strong>Loops</strong> use a two-digit helix prefix followed by a
              distance:{" "}
              <code className="modal-inline-code">xy.dist</code>, where{" "}
              <code className="modal-inline-code">x</code> is the helix the
              residue is closer to and{" "}
              <code className="modal-inline-code">y</code> is the more distant
              flanking helix.{" "}
              <code className="modal-inline-code">dist</code> counts residues
              from helix <code className="modal-inline-code">x</code> into the
              loop, zero-padded to three digits. So a residue in the TM1–TM2
              loop sitting close to TM1 reads e.g.{" "}
              <code className="modal-inline-code">12.004</code>, while one
              closer to TM2 reads e.g.{" "}
              <code className="modal-inline-code">21.003</code>.
            </li>
            <li>
              <strong>Gaps and insertions</strong> are handled explicitly: a
              missing residue skips the position (e.g.{" "}
              <code className="modal-inline-code">6.48</code> in ChRmine), and a
              bulged insertion receives a fractional identifier (e.g.{" "}
              <code className="modal-inline-code">5.451</code> in TaraRRB-R1).
            </li>
          </ul>
          <p className="modal-cite-body modal-cite-aside">
            For motif-level interpretation per protein family — DTD/FSE/NDQ/NTQ
            motifs, the channelrhodopsin TM2 glutamates and DC gate, sensory
            signaling residues, etc. — see the <em>GRN</em> tab.
          </p>
        </section>
        <section className="modal-section">
          <p className="modal-cite-label">Citation</p>
          <p className="modal-cite-body">
            If you use this resource in your work, please cite the MO numbering
            system paper (Hidber, Tajima, Kishi et al.,{" "}
            <em>manuscript under review</em>). A DOI will be added here once it is
            available.
          </p>
        </section>
        <section className="modal-section">
          <p className="modal-cite-label">Contact</p>
          <ul className="modal-contact-list">
            {contacts.map((c) => (
              <li key={c.email} className="modal-contact-item">
                <div className="modal-contact-name">{c.name}</div>
                <a
                  className="modal-contact-email"
                  href={`mailto:${c.email}`}
                >
                  {c.email}
                </a>
                <div className="modal-contact-topic">{c.topic}</div>
              </li>
            ))}
          </ul>
        </section>
        <section className="modal-section">
          <p className="modal-cite-label">Data</p>
          <p className="modal-cite-body">
            Catalog and aligned-structure data: Zenodo records{" "}
            <a
              href="https://doi.org/10.5281/zenodo.18147121"
              target="_blank"
              rel="noopener noreferrer"
            >
              10.5281/zenodo.18147121
            </a>{" "}
            and{" "}
            <a
              href="https://doi.org/10.5281/zenodo.20328414"
              target="_blank"
              rel="noopener noreferrer"
            >
              10.5281/zenodo.20328414
            </a>{" "}
            (CC-BY-4.0). Backend uses the{" "}
            <a
              href="https://github.com/flurinh/protos"
              target="_blank"
              rel="noopener noreferrer"
            >
              ProtOS
            </a>{" "}
            framework.
          </p>
        </section>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
//  Structure view (3Dmol viewer + multi-select + GRN selection bar)
// ---------------------------------------------------------------------------
const STRUCT_OVERLAY_COLORS = [
  "#0a0a0a",
  "#c81e1e",
  "#1d4ed8",
  "#057a55",
  "#d97706",
  "#7c3aed",
  "#0e7490",
  "#a16207",
];

function StructureView({
  manifest,
  selected,
  onToggle,
  data,
  loading,
  selectedGrns,
  onToggleGrn,
  onClearGrns,
  onSetSelectedGrns,
  grnFunctions,
}: {
  manifest: StructureManifestEntry[];
  selected: string[];
  onToggle: (id: string) => void;
  data: Record<string, StructureResponse>;
  loading: Set<string>;
  selectedGrns: Set<string>;
  onToggleGrn: (grn: string, additive: boolean) => void;
  onClearGrns: () => void;
  onSetSelectedGrns: React.Dispatch<React.SetStateAction<Set<string>>>;
  grnFunctions: GRNCategoryEntry[] | null;
}) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<any>(null);
  const [search, setSearch] = useState("");
  const [dropdownOpen, setDropdownOpen] = useState(false);
  const dropdownRef = useRef<HTMLDivElement>(null);

  // Initialize the 3Dmol viewer once (lazy-loaded on Structure-tab entry).
  const [viewerReady, setViewerReady] = useState(false);
  useEffect(() => {
    let cancelled = false;
    let handleResize: (() => void) | null = null;
    (async () => {
      if (!containerRef.current || viewerRef.current) return;
      try {
        await load3Dmol();
      } catch (e: any) {
        console.error("3Dmol load failed:", e);
        return;
      }
      if (cancelled || !containerRef.current) return;
      const lib = window.$3Dmol;
      const v = lib.createViewer(containerRef.current, {
        backgroundColor: "#ffffff",
        antialias: true,
      });
      viewerRef.current = v;
      handleResize = () => v.resize();
      window.addEventListener("resize", handleResize);
      setViewerReady(true);
    })();
    return () => {
      cancelled = true;
      if (handleResize) window.removeEventListener("resize", handleResize);
    };
  }, []);

  // Close dropdown when clicking outside.
  useEffect(() => {
    if (!dropdownOpen) return;
    const onClick = (e: MouseEvent) => {
      if (dropdownRef.current && !dropdownRef.current.contains(e.target as Node)) {
        setDropdownOpen(false);
      }
    };
    document.addEventListener("mousedown", onClick);
    return () => document.removeEventListener("mousedown", onClick);
  }, [dropdownOpen]);

  // Rebuild scene whenever selection / data / highlights change.
  useEffect(() => {
    const v = viewerRef.current;
    if (!v) return;
    v.removeAllModels();
    v.removeAllLabels();

    // Preserve the user's selection order so model indices are predictable.
    const loadedIds = selected.filter((id) => data[id]);
    loadedIds.forEach((id, idx) => {
      const struct = data[id];
      const overlayColor = STRUCT_OVERLAY_COLORS[idx % STRUCT_OVERLAY_COLORS.length];
      v.addModel(struct.pdb, "pdb");

      // Default cartoon: per-structure tint when overlaying >1, else region color.
      if (loadedIds.length > 1) {
        v.setStyle({ model: idx }, { cartoon: { color: overlayColor, opacity: 0.85 } });
      } else {
        // Single structure: color cartoon by TM region using each residue's GRN.
        v.setStyle({ model: idx }, { cartoon: { color: "#cccccc" } });
        Object.entries(struct.residue_grn).forEach(([resi, grn]) => {
          v.setStyle(
            { model: idx, resi: parseInt(resi, 10) },
            { cartoon: { color: regionColor(grn) } }
          );
        });
      }

      // Highlight selected GRNs as sticks across every loaded structure,
      // keeping each residue's existing color (region color in single-structure
      // mode, per-structure overlay color in multi-structure mode).
      if (selectedGrns.size > 0) {
        Object.entries(struct.residue_grn)
          .filter(([, g]) => selectedGrns.has(g))
          .forEach(([resi, grn]) => {
            const color = loadedIds.length > 1 ? overlayColor : regionColor(grn);
            v.setStyle(
              { model: idx, resi: parseInt(resi, 10) },
              {
                cartoon: { color },
                stick: { color, radius: 0.2 },
              }
            );
          });
      }

      // Hide retinal HETATMs for now (keeps the viewer focused on the cartoon).
      v.setStyle({ model: idx, hetflag: true }, { stick: { color: "#a16207", radius: 0.12 } });
    });

    // Hover labels using residue → GRN lookup, indexed by model.
    const idForModel = (m: any): string | null => {
      const idx = typeof m === "number" ? m : m;
      return loadedIds[idx] || null;
    };
    v.setHoverable(
      {},
      true,
      (atom: any) => {
        if (!atom || atom.label) return;
        const sid = idForModel(atom.model);
        if (!sid) return;
        const struct = data[sid];
        const grn = struct?.residue_grn[String(atom.resi)];
        const text = `${atom.resn}${atom.resi}${grn ? "  ·  " + grn : ""}${
          loadedIds.length > 1 ? "  ·  " + sid : ""
        }`;
        atom.label = v.addLabel(text, {
          position: { x: atom.x, y: atom.y, z: atom.z },
          backgroundColor: "#0a0a0a",
          fontColor: "#ffffff",
          fontSize: 11,
          padding: 4,
          borderThickness: 0,
        });
      },
      (atom: any) => {
        if (atom && atom.label) {
          v.removeLabel(atom.label);
          delete atom.label;
        }
      }
    );

    // Click any residue → toggle that residue's GRN in the selection set.
    // Same effect as clicking a cell in the GRN bar (additive toggle).
    v.setClickable({}, true, (atom: any) => {
      const sid = idForModel(atom.model);
      if (!sid) return;
      const struct = data[sid];
      const grn = struct?.residue_grn[String(atom.resi)];
      if (grn) onToggleGrn(grn, true);
    });

    v.zoomTo();
    v.render();
  }, [selected, data, selectedGrns, onToggleGrn, viewerReady]);

  const filteredManifest = useMemo(() => {
    const q = search.trim().toLowerCase();
    if (!q) return manifest;
    return manifest.filter((e) =>
      [e.structure_id, e.name, e.display_name, e.species, e.function, e.pdb_id]
        .filter(Boolean)
        .some((s) => (s as string).toLowerCase().includes(q))
    );
  }, [manifest, search]);

  // Build the GRN bar: union of every loaded structure's annotated GRNs, sorted N→C.
  const grnBar = useMemo(() => {
    const present = new Set<string>();
    selected.forEach((id) => {
      const d = data[id];
      if (!d) return;
      Object.values(d.residue_grn).forEach((g) => present.add(g));
    });
    return Array.from(present).sort((a, b) => grn2float(a) - grn2float(b));
  }, [selected, data]);

  const manifestById = useMemo(() => {
    const m = new Map<string, StructureManifestEntry>();
    manifest.forEach((e) => m.set(e.structure_id, e));
    return m;
  }, [manifest]);

  // Functions of currently loaded structures — used to filter GRN descriptions.
  const loadedFunctions = useMemo(() => {
    const fns = new Set<string>();
    selected.forEach((id) => {
      const m = manifestById.get(id);
      if (m?.function) fns.add(m.function);
    });
    return fns;
  }, [selected, manifestById]);

  const loadedFunctionsList = useMemo(
    () => Array.from(loadedFunctions).sort(),
    [loadedFunctions]
  );

  // For 2+ loaded structures: map each structure id to its assigned overlay
  // colour (same cycle the viewer effect uses). Drives the in-viewport legend.
  const overlayLegend = useMemo(() => {
    const loadedIds = selected.filter((id) => data[id]);
    if (loadedIds.length < 2) return null;
    return loadedIds.map((id, idx) => ({
      id,
      color: STRUCT_OVERLAY_COLORS[idx % STRUCT_OVERLAY_COLORS.length],
      label: displayStructureId(id),
      name: manifestById.get(id)?.name || id,
    }));
  }, [selected, data, manifestById]);

  // GRNs that exist in at least one loaded structure AND have at least one curated
  // position relevant to the loaded function(s). Powers the quick-pick dropdown.
  const functionalResidues = useMemo<
    { grn: string; label: string; categoryName: string; categoryColor: string }[]
  >(() => {
    if (!grnFunctions) return [];
    const presentInStructures = new Set<string>();
    selected.forEach((id) => {
      const d = data[id];
      if (!d) return;
      Object.values(d.residue_grn).forEach((g) => presentInStructures.add(g));
    });
    const isRelevant = (cat: GRNCategoryEntry) => {
      if (!cat.relevant_for || cat.relevant_for.length === 0) return true;
      if (cat.relevant_for.includes("*")) return true;
      if (loadedFunctions.size === 0) return true;
      return cat.relevant_for.some((f) => loadedFunctions.has(f));
    };
    const seen = new Map<
      string,
      { grn: string; label: string; categoryName: string; categoryColor: string }
    >();
    for (const cat of grnFunctions) {
      if (!isRelevant(cat)) continue;
      for (const pos of cat.positions) {
        if (!presentInStructures.has(pos.grn)) continue;
        if (!seen.has(pos.grn)) {
          seen.set(pos.grn, {
            grn: pos.grn,
            label: pos.label,
            categoryName: cat.name,
            categoryColor: cat.color,
          });
        }
      }
    }
    return Array.from(seen.values()).sort((a, b) => grn2float(a.grn) - grn2float(b.grn));
  }, [grnFunctions, selected, data, loadedFunctions]);

  // For each selected GRN: list of (category, position) entries from the function
  // map, restricted to categories applicable to any loaded structure's function.
  type PanelEntry = { grn: string; category: GRNCategoryEntry; position: GRNPositionEntry };
  const panelEntries = useMemo<PanelEntry[]>(() => {
    if (!grnFunctions || selectedGrns.size === 0) return [];
    const isRelevant = (cat: GRNCategoryEntry) => {
      if (!cat.relevant_for || cat.relevant_for.length === 0) return true;
      if (cat.relevant_for.includes("*")) return true;
      if (loadedFunctions.size === 0) return true; // no structures loaded → show everything
      return cat.relevant_for.some((f) => loadedFunctions.has(f));
    };
    const sortedGrns = Array.from(selectedGrns).sort((a, b) => grn2float(a) - grn2float(b));
    const out: PanelEntry[] = [];
    for (const grn of sortedGrns) {
      for (const cat of grnFunctions) {
        if (!isRelevant(cat)) continue;
        for (const pos of cat.positions) {
          if (pos.grn === grn) out.push({ grn, category: cat, position: pos });
        }
      }
    }
    return out;
  }, [grnFunctions, selectedGrns, loadedFunctions]);

  return (
    <>
      <div className="struct-controls">
        <div className="struct-selector" ref={dropdownRef}>
          <div className="struct-pills">
            {selected.map((id) => {
              const entry = manifestById.get(id);
              const isLoading = loading.has(id);
              const fnColor = entry?.function
                ? FUNCTION_TONES[entry.function] || "#8a8780"
                : "#8a8780";
              return (
                <span
                  key={id}
                  className="struct-pill"
                  title={entry?.display_name || entry?.name || id}
                >
                  <span className="struct-pill-dot" style={{ background: fnColor }} />
                  <span className="struct-pill-label" style={{ color: fnColor }}>
                    {displayStructureId(id)}
                  </span>
                  {isLoading && <span className="struct-pill-loading">…</span>}
                  <button
                    type="button"
                    className="struct-pill-x"
                    onClick={() => onToggle(id)}
                    aria-label={`Remove ${id}`}
                  >
                    ×
                  </button>
                </span>
              );
            })}
            <input
              className="struct-pill-input"
              placeholder={selected.length === 0 ? "Add structures…" : "+ add"}
              value={search}
              onFocus={() => setDropdownOpen(true)}
              onChange={(e) => {
                setSearch(e.target.value);
                setDropdownOpen(true);
              }}
            />
          </div>
          {dropdownOpen && (
            <div className="struct-dropdown">
              {filteredManifest.length === 0 && (
                <div className="struct-dropdown-empty">No matches.</div>
              )}
              {filteredManifest.map((e) => {
                const isSelected = selected.includes(e.structure_id);
                const fnColor = e.function ? FUNCTION_TONES[e.function] || "#8a8780" : "#8a8780";
                return (
                  <button
                    key={e.structure_id}
                    type="button"
                    className={`struct-option ${isSelected ? "selected" : ""}`}
                    onClick={() => onToggle(e.structure_id)}
                    title={e.display_name || e.name}
                  >
                    <span className="struct-option-check">{isSelected ? "✓" : ""}</span>
                    <span className="struct-option-dot" style={{ background: fnColor }} />
                    <span className="struct-option-name">{displayStructureId(e.structure_id)}</span>
                    <span className="struct-option-sid">{e.name}</span>
                    <span className="struct-option-type">
                      {e.structure_type === "experimental" ? "exp" : "pred"}
                    </span>
                    <span className="struct-option-species">{e.species || ""}</span>
                  </button>
                );
              })}
            </div>
          )}
        </div>
        <span className="struct-counter">
          {selected.length} loaded · {manifest.length} total
        </span>
      </div>

      {functionalResidues.length > 0 && (
        <div className="struct-quickpick">
          <label htmlFor="struct-quickpick-select" className="struct-quickpick-label">
            Functional residues
          </label>
          <select
            id="struct-quickpick-select"
            className="struct-quickpick-select"
            value=""
            onChange={(e) => {
              const v = e.target.value;
              if (v) onToggleGrn(v, true);
              e.target.value = "";
            }}
          >
            <option value="">
              {functionalResidues.length} available · jump to…
            </option>
            {functionalResidues.map((r) => (
              <option key={r.grn} value={r.grn}>
                {r.grn} — {r.label} ({r.categoryName})
              </option>
            ))}
          </select>
        </div>
      )}

      {loadedFunctionsList.length > 0 && (
        <div className="struct-legend">
          <span className="struct-legend-label">Families loaded</span>
          {loadedFunctionsList.map((fn) => {
            const c = FUNCTION_TONES[fn] || "#8a8780";
            return (
              <span
                key={fn}
                className="struct-legend-chip"
                style={{ borderColor: c, color: c }}
              >
                {fn}
              </span>
            );
          })}
        </div>
      )}

      <div className={`struct-stage ${panelEntries.length > 0 ? "with-panel" : ""}`}>
        <div className="struct-viewer-wrap">
          <div ref={containerRef} className="struct-viewer" />
          {overlayLegend && (
            <div className="struct-viewer-legend" role="legend">
              {overlayLegend.map((o) => (
                <div key={o.id} className="struct-viewer-legend-item" title={o.name}>
                  <span
                    className="struct-viewer-legend-swatch"
                    style={{ background: o.color }}
                  />
                  <span className="struct-viewer-legend-id">{o.label}</span>
                </div>
              ))}
            </div>
          )}
        </div>
        {panelEntries.length > 0 && (
          <aside className="struct-panel">
            {(() => {
              // Group panelEntries by GRN — each group becomes its own card.
              const groups: { grn: string; residue?: string | null; items: PanelEntry[] }[] = [];
              for (const e of panelEntries) {
                const last = groups[groups.length - 1];
                if (last && last.grn === e.grn) {
                  last.items.push(e);
                  if (!last.residue && e.position.residue) last.residue = e.position.residue;
                } else {
                  groups.push({ grn: e.grn, residue: e.position.residue, items: [e] });
                }
              }
              return groups.map((g) => (
                <article key={g.grn} className="struct-card">
                  <header className="struct-card-header">
                    <span className="struct-card-grn">{g.grn}</span>
                    {g.residue && <span className="struct-card-residue">{g.residue}</span>}
                  </header>
                  {g.items.map((it, i) => (
                    <div key={`${g.grn}-${it.category.key}-${i}`} className="struct-card-entry">
                      <div className="struct-card-category">
                        <span
                          className="struct-card-category-dot"
                          style={{ background: it.category.color }}
                        />
                        {it.category.name}
                      </div>
                      <div className="struct-card-label">{it.position.label}</div>
                      <p className="struct-card-desc">{it.position.description}</p>
                    </div>
                  ))}
                </article>
              ));
            })()}
          </aside>
        )}
      </div>

      <GRNSelectionBar
        grnBar={grnBar}
        selected={selected}
        data={data}
        manifestById={manifestById}
        selectedGrns={selectedGrns}
        onToggleGrn={onToggleGrn}
        onClearGrns={onClearGrns}
        onSetSelectedGrns={onSetSelectedGrns}
        grnFunctions={grnFunctions}
        loadedFunctions={loadedFunctions}
      />
    </>
  );
}

// ---------------------------------------------------------------------------
//  GRN selection bar: GRN cells (drag-select + tooltips) + per-structure rows
// ---------------------------------------------------------------------------
function GRNSelectionBar({
  grnBar,
  selected,
  data,
  manifestById,
  selectedGrns,
  onToggleGrn,
  onClearGrns,
  onSetSelectedGrns,
  grnFunctions,
  loadedFunctions,
}: {
  grnBar: string[];
  selected: string[];
  data: Record<string, StructureResponse>;
  manifestById: Map<string, StructureManifestEntry>;
  selectedGrns: Set<string>;
  onToggleGrn: (grn: string, additive: boolean) => void;
  onClearGrns: () => void;
  onSetSelectedGrns: React.Dispatch<React.SetStateAction<Set<string>>>;
  grnFunctions: GRNCategoryEntry[] | null;
  loadedFunctions: Set<string>;
}) {
  const curatedLabel = useMemo<Record<string, string>>(() => {
    const out: Record<string, string> = {};
    if (!grnFunctions) return out;
    const isRelevant = (cat: GRNCategoryEntry) => {
      if (!cat.relevant_for || cat.relevant_for.length === 0) return true;
      if (cat.relevant_for.includes("*")) return true;
      if (loadedFunctions.size === 0) return true;
      return cat.relevant_for.some((f) => loadedFunctions.has(f));
    };
    for (const cat of grnFunctions) {
      if (!isRelevant(cat)) continue;
      for (const pos of cat.positions) {
        if (out[pos.grn]) continue;
        out[pos.grn] = `${pos.label} · ${cat.name}`;
      }
    }
    return out;
  }, [grnFunctions, loadedFunctions]);

  const [dragStartIdx, setDragStartIdx] = useState<number | null>(null);
  const [dragCurrentIdx, setDragCurrentIdx] = useState<number | null>(null);
  const [dragAdditive, setDragAdditive] = useState(false);
  const dragMovedRef = useRef(false);

  const onCellMouseDown = (idx: number, e: React.MouseEvent) => {
    setDragStartIdx(idx);
    setDragCurrentIdx(idx);
    setDragAdditive(e.shiftKey || e.metaKey || e.ctrlKey);
    dragMovedRef.current = false;
  };
  const onCellMouseEnter = (idx: number) => {
    if (dragStartIdx == null) return;
    setDragCurrentIdx(idx);
    if (idx !== dragStartIdx) dragMovedRef.current = true;
  };

  useEffect(() => {
    if (dragStartIdx == null) return;
    const onUp = () => {
      const start = dragStartIdx;
      const end = dragCurrentIdx ?? dragStartIdx;
      const lo = Math.min(start, end);
      const hi = Math.max(start, end);
      const moved = dragMovedRef.current;
      if (!moved && lo === hi) {
        onToggleGrn(grnBar[lo], dragAdditive);
      } else {
        const rangeGrns = grnBar.slice(lo, hi + 1);
        onSetSelectedGrns((prev) => {
          if (!dragAdditive) return new Set(rangeGrns);
          const next = new Set<string>(prev);
          rangeGrns.forEach((g) => next.add(g));
          return next;
        });
      }
      setDragStartIdx(null);
      setDragCurrentIdx(null);
      dragMovedRef.current = false;
    };
    document.addEventListener("mouseup", onUp);
    return () => document.removeEventListener("mouseup", onUp);
  }, [dragStartIdx, dragCurrentIdx, dragAdditive, grnBar, onToggleGrn, onSetSelectedGrns]);

  const dragRangeSet = useMemo(() => {
    if (dragStartIdx == null || dragCurrentIdx == null) return null;
    const lo = Math.min(dragStartIdx, dragCurrentIdx);
    const hi = Math.max(dragStartIdx, dragCurrentIdx);
    const s = new Set<number>();
    for (let i = lo; i <= hi; i++) s.add(i);
    return s;
  }, [dragStartIdx, dragCurrentIdx]);

  const loadedRows = selected.filter((id) => data[id]);

  return (
    <div className="struct-grnbar">
      <div className="struct-grnbar-meta">
        <span>
          {grnBar.length} GRN positions · {selectedGrns.size} selected
        </span>
        {selectedGrns.size > 0 && (
          <button className="link-btn" onClick={onClearGrns}>
            Clear
          </button>
        )}
      </div>
      <div className="struct-grnbar-scroll">
        <table className="struct-grnbar-table">
          <tbody>
            {/* GRN header row */}
            <tr className="struct-grnbar-row-grn">
              <th className="struct-grnbar-rowname" scope="row">
                <span className="struct-grnbar-rowname-label">GRN</span>
              </th>
              {grnBar.map((g, i) => {
                const color = regionColor(g);
                const isSelected = selectedGrns.has(g);
                const inDrag = dragRangeSet?.has(i) ?? false;
                const showLabel = i === 0 || regionLabel(g) !== regionLabel(grnBar[i - 1]);
                const curated = curatedLabel[g];
                const title = curated
                  ? `${g} · ${regionLabel(g)}\n${curated}`
                  : `${g} · ${regionLabel(g)}`;
                return (
                  <td
                    key={g}
                    className={`struct-grnbar-cell${isSelected ? " selected" : ""}${
                      inDrag && !isSelected ? " dragging" : ""
                    }${curated ? " has-curated" : ""}`}
                    style={{ borderBottomColor: color }}
                    onMouseDown={(e) => onCellMouseDown(i, e)}
                    onMouseEnter={() => onCellMouseEnter(i)}
                    title={title}
                  >
                    <span className="struct-grnbar-region" style={{ color }}>
                      {showLabel ? regionLabel(g) : " "}
                    </span>
                    <span className="struct-grnbar-grn">{g}</span>
                  </td>
                );
              })}
            </tr>
            {/* Per-structure residue rows */}
            {loadedRows.map((id) => {
              const struct = data[id];
              const entry = manifestById.get(id);
              const shortName = entry?.name || id;
              const fnColor = entry?.function
                ? FUNCTION_TONES[entry.function] || "var(--ink)"
                : "var(--ink)";
              return (
                <tr key={id} className="struct-grnbar-row-res">
                  <th className="struct-grnbar-rowname" scope="row" title={entry?.display_name || shortName}>
                    <span className="struct-grnbar-rowname-label" style={{ color: fnColor }}>
                      {shortName}
                    </span>
                  </th>
                  {grnBar.map((g, i) => {
                    const hit = struct.grn_residue[g];
                    const isSelected = selectedGrns.has(g);
                    const inDrag = dragRangeSet?.has(i) ?? false;
                    const baseClass = `struct-resrow-cell${isSelected ? " selected" : ""}${
                      inDrag && !isSelected ? " dragging" : ""
                    }`;
                    if (!hit) {
                      return (
                        <td
                          key={g}
                          className={`${baseClass} gap`}
                          onMouseDown={(e) => onCellMouseDown(i, e)}
                          onMouseEnter={() => onCellMouseEnter(i)}
                          title={`${shortName} · — · ${g}`}
                        >
                          —
                        </td>
                      );
                    }
                    const style = residueStyle(hit.res1);
                    return (
                      <td
                        key={g}
                        className={baseClass}
                        style={{ background: style.background, color: style.color }}
                        onMouseDown={(e) => onCellMouseDown(i, e)}
                        onMouseEnter={() => onCellMouseEnter(i)}
                        title={`${shortName} · ${hit.res1}${hit.resi} · ${g}`}
                      >
                        <span className="struct-resrow-res">{hit.res1}</span>
                        <span className="struct-resrow-num">{hit.resi}</span>
                      </td>
                    );
                  })}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
//  GRN function-map view
// ---------------------------------------------------------------------------
function GRNFunctionsView({
  categories,
  visible,
  activeKey,
  onActiveKey,
  query,
  onQuery,
}: {
  categories: GRNCategoryEntry[];
  visible: GRNCategoryEntry[];
  activeKey: string | null;
  onActiveKey: (k: string | null) => void;
  query: string;
  onQuery: (q: string) => void;
}) {
  return (
    <>
      <div className="grnfn-controls">
        <input
          type="search"
          className="catalog-search"
          placeholder="Search GRN, residue, term…"
          value={query}
          onChange={(e) => onQuery(e.target.value)}
        />
        <div className="catalog-fn-filter">
          <button
            className={`catalog-fn-chip ${activeKey === null ? "active" : ""}`}
            onClick={() => onActiveKey(null)}
          >
            All
          </button>
          {categories.map((c) => {
            const active = activeKey === c.key;
            return (
              <button
                key={c.key}
                className={`catalog-fn-chip ${active ? "active" : ""}`}
                style={{
                  borderColor: c.color,
                  background: active ? c.color : "transparent",
                  color: active ? "#fff" : c.color,
                }}
                onClick={() => onActiveKey(active ? null : c.key)}
              >
                {c.name}
              </button>
            );
          })}
        </div>
      </div>

      <div className="grnfn-grid">
        {visible.map((c) => (
          <article key={c.key} className="grnfn-card">
            <header className="grnfn-card-header">
              <span className="grnfn-rule" style={{ background: c.color }} />
              <h2 className="grnfn-card-title">{c.name}</h2>
            </header>
            <p className="grnfn-summary">{c.summary}</p>
            <ul className="grnfn-positions">
              {c.positions.map((p) => (
                <li key={`${c.key}-${p.grn}-${p.label}`} className="grnfn-position">
                  <div className="grnfn-position-key">
                    <span className="grnfn-grn" style={{ color: c.color }}>
                      {p.grn}
                    </span>
                    {p.residue && <span className="grnfn-residue">{p.residue}</span>}
                  </div>
                  <div className="grnfn-position-body">
                    <span className="grnfn-label">{p.label}</span>
                    <span className="grnfn-description">{p.description}</span>
                  </div>
                </li>
              ))}
            </ul>
          </article>
        ))}
        {visible.length === 0 && (
          <p className="muted">No matches.</p>
        )}
      </div>
    </>
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
