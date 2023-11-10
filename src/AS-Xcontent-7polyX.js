// @ts-check

const fs = require("fs");

const { ArgvParser, } = require("./ArgvParser.js");
const { readFasta, } = require("./fasta_util.js");
const { findPolyXX, } = require("./find_polyX.js");
const { find_SCD_v4, MAX_DISTANCE, } = require("./find_SCD.js");


class ASD {
	/** @param {string | string[]} poly */
	get_name(poly) { return ""; }

	/**
	 * @param {number} num_AA
	 * @param {number} polyQ_region_length
	 */
	 cond(num_AA, polyQ_region_length) { return false; }
}

/** @type {ASD[]} */
const polyAA_subTypes = [
	{ get_name(poly) { return     `4${Array.isArray(poly) ? poly.join("⁄") : poly}(4⁄4)`;     }, cond: (x, y) => (  x  ) >= 4   && y == 4, },
	{ get_name(poly) { return     `5${Array.isArray(poly) ? poly.join("⁄") : poly}(4⁄5-5⁄5)`; }, cond: (x, y) => (  x  ) >= 4   && y == 5, },
	{ get_name(poly) { return     `6${Array.isArray(poly) ? poly.join("⁄") : poly}(4⁄6-6⁄6)`; }, cond: (x, y) => (  x  ) >= 4   && y == 6, },
	{ get_name(poly) { return     `7${Array.isArray(poly) ? poly.join("⁄") : poly}(4⁄7-7⁄7)`; }, cond: (x, y) => (  x  ) >= 4   && y == 7, },
	{ get_name(poly) { return  `8-10${Array.isArray(poly) ? poly.join("⁄") : poly}(≥50%)`;    }, cond: (x, y) => (x / y) >= 0.5 && y >= 8  && y <= 10, },
	{ get_name(poly) { return `11-20${Array.isArray(poly) ? poly.join("⁄") : poly}(≥50%)`;    }, cond: (x, y) => (x / y) >= 0.5 && y >= 11 && y <= 20, },
	{ get_name(poly) { return   `≥21${Array.isArray(poly) ? poly.join("⁄") : poly}(≥50%)`;    }, cond: (x, y) => (x / y) >= 0.5 && y >= 21,            },
];

/**
 * @param {Record<string, string} fa
 * @param {(string[] | string)[]} poly_task ["Q", "N", ["Q", "N"]]
 * @param {{ has_subtype: boolean; has_SCD: boolean; }} options
 */
function print_polyQ_resdue(fa, poly_task = ["Q", "N", ["Q", "N"]], options = { has_subtype: true, has_SCD: false, }) {
	const seq_names = Object.keys(fa);
	
	const n_total_protein = seq_names.length;
	const results = {
		"protein #": n_total_protein, // total protein #
	};
	if (options.has_SCD) {
		results["SCD proteins #"] = 0;
		Object.values(fa).forEach((seq, seq_idx) => {
			const scd_list = find_SCD_v4(seq, MAX_DISTANCE);
			if (scd_list.length > 0) {
				results["SCD proteins #"] += 1;
			}
		});
		results["SCD proteins %"] = results["SCD proteins #"] / n_total_protein;
	}

	poly_task.map(aa => {
		const aa_name = Array.isArray(aa) ? aa.join("⁄") : aa;
		results[`total Poly${aa_name} proteins #`] = 0;
		results[`total Poly${aa_name} proteins %`] = 0;

		if (options.has_subtype) {
			polyAA_subTypes.forEach(asd => {
				const name = asd.get_name(aa);
				results[`${name} proteins #`] = 0;
				results[`${name} proteins %`] = 0;
			});
		}
	});

	Object.values(fa).forEach((seq, seq_idx) => {
		poly_task.map(aa => {
			const aa_name = Array.isArray(aa) ? aa.join("⁄") : aa;
			const _poly_type = findPolyXX(seq, 4, aa).map(v => {
				const {
					num_AA: num_AA,
					length: polyQ_region_length,
					} = v;
				const poly_type = polyAA_subTypes.find(d => d.cond(num_AA, polyQ_region_length));
				if (poly_type) {
					return {
						num_AA,
						polyQ_region_length,
						subtype: poly_type.get_name(aa),
					};
				}
			});
			const poly_type_list = new Set(_poly_type.filter(a => a).map(a => a.subtype));//.sort((a, b) => b.polyQ_region_length - a.polyQ_region_length)[0];
			if (poly_type_list.size) {
				results[`total Poly${aa_name} proteins #`] += 1;
				if (options.has_subtype) {
					[...poly_type_list].forEach(poly_type => {
						results[`${poly_type} proteins #`] += 1;
					});
				}
			}
		});
	});

	poly_task.map(aa => {
		const aa_name = Array.isArray(aa) ? aa.join("⁄") : aa;
		results[`total Poly${aa_name} proteins %`] =  results[`total Poly${aa_name} proteins #`] / n_total_protein;
		if (options.has_subtype) {
			polyAA_subTypes.forEach(asd => {
				const name = asd.get_name(aa);
				results[`${name} proteins %`] = results[`${name} proteins #`] / n_total_protein;
			});
		}
	});

	return results;
}

function cmd_main() {
	const argv = new ArgvParser(process.argv);
	const aa = argv.get(/--aa=(.+)/, (arg, args) => args[0] ? args[0] : null);
	const path_to_fa = argv.get(/--fa=(.+)/, (arg, args) => args[0] ? args[0] : null);
	if (path_to_fa && aa) {
		const fa = readFasta(path_to_fa);

		const results = print_polyQ_resdue(fa, [[...aa]], {
			has_subtype: true,
			has_SCD: false,
		});
	
		const header_cols = Object.keys(results);
		const fin = header_cols.map(k => results[k]).join("\t");
	
		console.log(header_cols.join("\t"));
		console.log(fin);
	}
	else {
		console.error("Usage:", "AS-Xcontent-7polyX --aa=Q --fa=protein.fa");
		console.error("\t--aa: amino acid symbol.");
		console.error("\t--fa: protein fasta file");
	}
}
if (__filename == fs.realpathSync(process.argv[1])) {
	cmd_main();
}

