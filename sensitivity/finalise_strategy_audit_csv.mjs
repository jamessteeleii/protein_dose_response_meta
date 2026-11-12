import fs from "node:fs/promises";
import path from "node:path";

const { Workbook } = await import(
  process.env.OAI_ARTIFACT_TOOL_MODULE ?? "@oai/artifact-tool"
);

const inputPath = process.argv[2] ??
  "data/protein_manipulation_strategy_audit.csv";
const outputPath = process.argv[3] ?? inputPath;

const finalHeaders = [
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

const csvText = await fs.readFile(inputPath, "utf8");
const workbook = await Workbook.fromCSV(csvText, { sheetName: "Strategy audit" });
const matrix = workbook.worksheets.getItem("Strategy audit").getUsedRange(true).values;

if (!Array.isArray(matrix) || matrix.length !== 110) {
  throw new Error("Expected one header row and 109 study-audit rows.");
}

const sourceHeaders = matrix[0].map((value) => String(value ?? "").trim());
const missing = finalHeaders.filter((header) => !sourceHeaders.includes(header));
if (missing.length > 0) {
  throw new Error(`Missing final audit columns: ${missing.join(", ")}`);
}

const sourceRows = matrix.slice(1).map((values) =>
  Object.fromEntries(
    sourceHeaders.map((header, index) => [header, values[index]]),
  ),
);
const normalize = (value) => String(value ?? "").trim();
const strategies = new Set([
  "Addition",
  "Substitution",
  "Mixed",
  "Protein-only comparison",
  "Unclear",
]);
const confidenceLevels = new Set(["High", "Moderate", "Low"]);
const clearlyMatched = new Set([
  "Matched",
  "Matched/controlled",
  "Approximately matched",
  "Matched by powder dose",
  "Probably matched",
]);
const clearlyNotMatched = new Set(["Not matched", "Not constrained"]);

const studyNumbers = sourceRows.map((row) => normalize(row.study_no));
if (
  studyNumbers.some((value) => !/^\d+$/.test(value) || Number(value) < 1) ||
  new Set(studyNumbers).size !== sourceRows.length
) {
  throw new Error("Study numbers must be unique positive integers.");
}

if (
  sourceHeaders.includes("review_status") &&
  sourceRows.some((row) => normalize(row.review_status) !== "Confirmed")
) {
  throw new Error("Every study must be Confirmed before finalisation.");
}

for (const row of sourceRows) {
  const study = normalize(row.study_no);
  if (!strategies.has(normalize(row.primary_strategy))) {
    throw new Error(`Study ${study} has an invalid primary strategy.`);
  }
  if (!confidenceLevels.has(normalize(row.confidence))) {
    throw new Error(`Study ${study} has an invalid confidence level.`);
  }
  const strategy = normalize(row.primary_strategy);
  const energyRelation = normalize(row.energy_relation);
  if (
    (strategy === "Addition" && clearlyMatched.has(energyRelation)) ||
    (strategy === "Substitution" && clearlyNotMatched.has(energyRelation))
  ) {
    throw new Error(
      `Study ${study} has a strategy/energy-relation contradiction.`,
    );
  }
  for (const field of [
    "authors",
    "year",
    "title",
    "arm_count",
    "arm_mix",
    "strategy_detail",
    "energy_relation",
    "displaced_energy",
    "coding_reason",
    "reported_intervention",
    "evidence_excerpt",
  ]) {
    if (normalize(row[field]) === "") {
      throw new Error(`Study ${study} is missing ${field}.`);
    }
  }
}

const outputRows = sourceRows.map((row) =>
  Object.fromEntries(finalHeaders.map((header) => [header, row[header]])),
);

const csvEscape = (value) => {
  const text = String(value ?? "").replace(/\r\n/g, "\n").replace(/\r/g, "\n");
  if (/[",\n]/.test(text)) return `"${text.replace(/"/g, '""')}"`;
  return text;
};
const outputCsv = [
  finalHeaders.join(","),
  ...outputRows.map((row) =>
    finalHeaders.map((header) => csvEscape(row[header])).join(","),
  ),
].join("\n") + "\n";

// Validate the lean file through the spreadsheet parser before replacing the
// user-reviewed source.
const roundTrip = await Workbook.fromCSV(outputCsv, { sheetName: "Strategy audit" });
const roundTripMatrix = roundTrip.worksheets
  .getItem("Strategy audit")
  .getUsedRange(true).values;
if (
  roundTripMatrix.length !== 110 ||
  roundTripMatrix[0].length !== finalHeaders.length ||
  roundTripMatrix[0].some((value, index) => value !== finalHeaders[index])
) {
  throw new Error("The final strategy CSV failed its round-trip audit.");
}

await fs.mkdir(path.dirname(outputPath), { recursive: true });
const temporaryPath = `${outputPath}.tmp`;
await fs.writeFile(temporaryPath, outputCsv, "utf8");
await fs.rename(temporaryPath, outputPath);

const countBy = (key) => Object.fromEntries(
  [...new Set(outputRows.map((row) => normalize(row[key])))].sort().map((value) => [
    value,
    outputRows.filter((row) => normalize(row[key]) === value).length,
  ]),
);
const resolved = (row) =>
  ["Addition", "Substitution"].includes(normalize(row.primary_strategy));

process.stdout.write(`${JSON.stringify({
  rows: outputRows.length,
  columns: finalHeaders.length,
  strategies: countBy("primary_strategy"),
  confidence: countBy("confidence"),
  principal_eligible: outputRows.filter(resolved).length,
  strict_eligible: outputRows.filter(
    (row) => resolved(row) && normalize(row.confidence) === "High",
  ).length,
  output: path.resolve(outputPath),
}, null, 2)}\n`);
