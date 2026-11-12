/** Run the same recommendation and optional pilot update as the browser. */
import fs from "node:fs";
import { evaluate } from "../docs/model-core.js";

if (process.argv.length > 4 || process.argv.includes("--help")) {
  console.log("node analysis/predict.mjs [settings.json [pilot-summary.json]]");
  process.exit(process.argv.includes("--help") ? 0 : 2);
}

const read = name => JSON.parse(fs.readFileSync(name, "utf8"));
const dyes = read(new URL("../docs/data/dyes.json", import.meta.url));
const config = read(new URL("../docs/data/model.json", import.meta.url));
const settings = process.argv[2] ? read(process.argv[2]) : {};
const pilot = process.argv[3] ? read(process.argv[3]) : null;
const decision = evaluate(dyes, config, settings, pilot);
console.log(JSON.stringify(decision, null, 2));
if (decision.status === "invalid") process.exitCode = 2;
