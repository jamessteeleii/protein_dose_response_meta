import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const [inputPath, outputPath] = process.argv.slice(2);

if (!inputPath || !outputPath) {
  throw new Error(
    "Usage: node export_strategy_audit_csv.mjs <audit.xlsx> <audit.csv>",
  );
}

const blob = await FileBlob.load(inputPath);
const workbook = await SpreadsheetFile.importXlsx(blob);
const sheet = workbook.worksheets.getItem("Study audit");
const matrix = sheet.getUsedRange(true).values;

if (!Array.isArray(matrix) || matrix.length !== 110) {
  throw new Error("Expected one header row and 109 study-audit rows.");
}

const headers = matrix[0].map((value) => String(value ?? "").trim());
const rows = matrix.slice(1).map((values) =>
  Object.fromEntries(headers.map((header, index) => [header, values[index]])),
);

const requiredHeaders = [
  "Study No",
  "PMID",
  "DOI",
  "Authors",
  "Year",
  "Title",
  "Arms",
  "Arm mix",
  "Primary strategy",
  "Strategy detail",
  "Energy relation",
  "Displaced energy",
  "Confidence",
  "Review status",
  "Coding reason",
  "Reported intervention",
  "PDF evidence excerpt",
  "Page",
  "Identifier match",
];

const missingHeaders = requiredHeaders.filter((header) => !headers.includes(header));
if (missingHeaders.length > 0) {
  throw new Error(`Missing audit columns: ${missingHeaders.join(", ")}`);
}

const outputHeaders = [
  "study_no",
  "pmid",
  "doi",
  "authors",
  "year",
  "title",
  "arm_count",
  "arm_mix",
  "primary_strategy",
  "strategy_detail",
  "energy_relation",
  "displaced_energy",
  "confidence",
  "coding_reason",
  "reported_intervention",
  "evidence_excerpt",
  "page",
];

const normalize = (value) => String(value ?? "").trim();
const isResolved = (strategy) => ["Addition", "Substitution"].includes(strategy);

const outputRows = rows.map((row) => {
  const strategy = normalize(row["Primary strategy"]);
  const reviewStatus = normalize(row["Review status"]);
  const confidence = normalize(row.Confidence);
  if (reviewStatus !== "Confirmed") {
    throw new Error(
      `Study ${normalize(row["Study No"])} is not confirmed and cannot be exported to the final audit.`,
    );
  }

  return {
    study_no: row["Study No"],
    pmid: row.PMID,
    doi: row.DOI,
    authors: row.Authors,
    year: row.Year,
    title: row.Title,
    arm_count: row.Arms,
    arm_mix: row["Arm mix"],
    primary_strategy: strategy,
    strategy_detail: row["Strategy detail"],
    energy_relation: row["Energy relation"],
    displaced_energy: row["Displaced energy"],
    confidence,
    coding_reason: row["Coding reason"],
    reported_intervention: row["Reported intervention"],
    evidence_excerpt: row["PDF evidence excerpt"],
    page: row.Page,
  };
});

const clearlyMatched = new Set([
  "Matched",
  "Matched/controlled",
  "Approximately matched",
  "Matched by powder dose",
  "Probably matched",
]);
const clearlyNotMatched = new Set(["Not matched", "Not constrained"]);
for (const row of outputRows) {
  if (
    (row.primary_strategy === "Addition" &&
      clearlyMatched.has(normalize(row.energy_relation))) ||
    (row.primary_strategy === "Substitution" &&
      clearlyNotMatched.has(normalize(row.energy_relation)))
  ) {
    throw new Error(
      `Study ${normalize(row.study_no)} has a strategy/energy-relation contradiction.`,
    );
  }
}

const studyNumbers = outputRows.map((row) => Number(row.study_no));
if (
  studyNumbers.some((value) => !Number.isInteger(value) || value < 1) ||
  new Set(studyNumbers).size !== outputRows.length
) {
  throw new Error("Study numbers must be unique positive integers.");
}

const csvEscape = (value) => {
  if (value === null || value === undefined) return "";
  if (typeof value === "boolean") return value ? "TRUE" : "FALSE";
  const text = String(value).replace(/\r\n/g, "\n").replace(/\r/g, "\n");
  if (/[",\n]/.test(text)) return `"${text.replace(/"/g, '""')}"`;
  return text;
};

const csv = [
  outputHeaders.join(","),
  ...outputRows.map((row) =>
    outputHeaders.map((header) => csvEscape(row[header])).join(","),
  ),
].join("\n") + "\n";

await fs.mkdir(path.dirname(outputPath), { recursive: true });
await fs.writeFile(outputPath, csv, "utf8");

// Round-trip the machine input through the spreadsheet parser. This catches
// quoting/newline errors in evidence excerpts before the CSV reaches targets.
const csvWorkbook = await Workbook.fromCSV(csv, { sheetName: "Strategy audit" });
const csvMatrix = csvWorkbook.worksheets
  .getItem("Strategy audit")
  .getUsedRange(true).values;
if (
  !Array.isArray(csvMatrix) ||
  csvMatrix.length !== outputRows.length + 1 ||
  csvMatrix[0].length !== outputHeaders.length
) {
  throw new Error("The exported strategy CSV failed its round-trip audit.");
}

const countBy = (key) =>
  Object.fromEntries(
    [...new Set(outputRows.map((row) => String(row[key])))].sort().map((value) => [
      value,
      outputRows.filter((row) => String(row[key]) === value).length,
    ]),
  );

process.stdout.write(
  `${JSON.stringify(
    {
      rows: outputRows.length,
      strategies: countBy("primary_strategy"),
      principal_eligible: outputRows.filter((row) =>
        isResolved(row.primary_strategy),
      ).length,
      strict_eligible: outputRows.filter((row) =>
        isResolved(row.primary_strategy) && row.confidence === "High",
      ).length,
      output: path.resolve(outputPath),
    },
    null,
    2,
  )}\n`,
);
