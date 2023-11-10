// @ts-check

const fs = require("fs");
const { strictEqual, notEqual, notStrictEqual } = require("assert");

const { ArgvParser, } = require("./ArgvParser.js");
const { readFasta, } = require("./fasta_util.js");

/**
 * @param {{ [seq_name: string]: string; }} fa_colle
 * @param {number|string} geneticCode
 * @param {null | { [seq_name: string]: string; }} [protein_fa]
 */
function _fasta_codonUsage(fa_colle, geneticCode, protein_fa = null, print_log = false) {
	const codon_tab = generateCodonTable(geneticCode);
	const keys = Object.keys(codon_tab);
	const rev_codon_tab = (function () {
		/** @type {{ [aa: string]: string[]; }} */
		const m = {};
		Object.entries(codon_tab).forEach(b => {
			const d = b[0];
			const a = b[1];
			m[a] ??= [];
			m[a].push(d);
		});
		return m;
	})();
	
	function gen_codon_usage_tab() {
		const tab = Object.fromEntries(Object.entries(codon_tab).map(([dna, aa]) => {
			const obj = {
				dna: dna,
				rna: dna.replace(/T/g, "U"),
				aa: aa,
				num: 0,
				w_p: 0,
				// s_p: "",
			};
			return [dna, obj];
		}));
		return tab;
	}

	const fin_tab = gen_codon_usage_tab();
	
	const co_aa_tab = generateAATable_co();
	const sp_aa_tab = generateAATable_sp();

	const diff_list = [];
	const sp_list = [];
	const unknow_list = [];

	let num_notEqual = 0;

	let total_proteome_aa = 0;

	Object.entries(fa_colle).forEach(([ fa_name, fa_seq, ]) => {
		const gg_tab = gen_codon_usage_tab();
		const seq = fa_seq.toUpperCase();

		const protein_length = seq.length / 3;

		let is_diff = false;

		if (protein_fa && !protein_fa[fa_name]) {
			console.log("undefined protein_fa <-", fa_name);
		}

		// assert.strictEqual(seq.length % 3, 0, `seq.length % 3 == 0; // ${seq.length} % 3 -> ${seq.length % 3}`);
		if (seq.length % 3 == 0 && (protein_fa == null || Math.abs(protein_length - protein_fa[fa_name]?.length) <= 1)) {
			for (let i = 0, j = 0; i < seq.length, j < protein_length; i += 3, ++j) {
				const codon = seq.slice(i, i + 3);
				const codon_no_N = codon.replace(/N/g, ".");
				const aa = protein_fa && protein_fa[fa_name] ? protein_fa[fa_name][j] : null;

				const sp_aa = sp_aa_tab[aa];
				if (sp_aa) {
					if (!sp_aa.alias || sp_aa.alias.some(a => rev_codon_tab[a].find(ccc => ccc == codon || ccc.match(codon_no_N)))) {
						const sp = {
							seq_id: fa_name,
							codon: codon,
							aa: aa,
							alias: sp_aa.alias,
						};
						sp_list.push(sp);
					}
					else {
						push_diff_codon();
					}
				}
				else if (gg_tab[codon]) {
				// else if (gg_tab[codon] || ccc.match(codon_no_N)) {
					// if (codon[0] != "N" && codon[1] != "N" && codon[2] != "N") {
						++gg_tab[codon].num;
					// }
					// else {
					// }

					if (protein_fa != null) {
						if (codon_tab[codon] == "*" && aa == null) {
						}
						else if (co_aa_tab[aa]?.alias?.includes?.(codon_tab[codon])) {
							const comboo = {
								seq_id: fa_name,
								codon: codon,
								aa_1: gg_tab[codon],
								aa_2: aa,
							};
							console.log("else if (co_aa_tab[aa]?.alias?.includes?.(codon_tab[codon])) {", comboo);
						}
						else if (codon_tab[codon] != aa) {//asdasdasd
							push_diff_codon();
							// ++num_notEqual;
							// return;
						}
					}
				}
				else if (print_log) {
					push_unknow_codon_N();
				}

				function push_unknow_codon_N() {
					const unknow = {
						seq_id: fa_name,
						codon: codon,
						// aa: gg_tab[codon],
						aa: aa,
						"CDS": `${i}/${seq.length}`,
						"protein": `${j}/${protein_length}`,
						"is last": (i - 3 + 1) == seq.length,
					};
					console.log("if (tab[codon]) {", unknow);
					unknow_list.push(unknow);
				}

				function push_diff_codon() {
					const except = {
						seq_id: fa_name,
						"aa not equal": true,
						"seq name": fa_name,
						"CDS pos": `${i} (${i / 3})`,
						"protein pos": `${j} (${j * 3})`,
						"codon": `${codon}`,
						"translate": `${codon_tab[codon]}`,
						"amino acid": `${aa}`, // expected
					};
					if (print_log) {
						console.log(except);
					}
					diff_list.push(except);
					is_diff = true;
				}
			}
		}
		else if (print_log) {
			console.log("if (seq.length % 3 == 0) {", fa_name, protein_fa ? protein_fa[fa_name]?.length : "NA", seq.length, Math.trunc(protein_length), seq.length % 3);
			return;
		}

		if (is_diff) {
			++num_notEqual;
		}
		total_proteome_aa += protein_length;
		keys.forEach(code => fin_tab[code].num += gg_tab[code].num);
	});

	// console.log(total_proteome_aa);

	const values = Object.values(fin_tab).sort((a, b) => a.aa.localeCompare(b.aa));
	// return values.map(data => [data.aa, data.dna, data.num, (100 * data.num / total_proteome_aa).toFixed(2) + "%"]);
	
	const _total_proteome_aa = values.reduce((acc, data) => acc + data.num, 0);
	if (total_proteome_aa != _total_proteome_aa) {
		console.log({ total_proteome_aa, _total_proteome_aa, });
	}

	values.forEach(data => {
		data.w_p = data.num / _total_proteome_aa;//(100 * data.num / total_proteome_aa).toFixed(2);
	});
	return {
		num_CDS_input: Object.keys(fa_colle).length,
		num_protein_input: protein_fa ? Object.keys(protein_fa).length : null,
		num_notEqual: num_notEqual,
		codonUsage: values,
		codonUsage_dna_p_map: Object.fromEntries(values.map(v => [v.dna, v.w_p])),
		diff_list: diff_list,
		unknow_list: unknow_list,
		sp_list: sp_list,
	};
}


/**
 * @param {string} cds_seq
 * @param {number|string} geneticCode
 */
 function cds_translate(cds_seq, geneticCode) {
	strictEqual(cds_seq.length % 3, 0, "cds_seq.length % 3 == 0");

	const codon_tab = generateCodonTable(geneticCode);
	// const pep_len = cds_seq.length / 3;

	let pep_seq = "";
	for (let i = 0; i < cds_seq.length; i += 3) {
		const aa = codon_tab[cds_seq.slice(i, i + 3)];
		notStrictEqual(aa, undefined);
		pep_seq += aa;
	}

	return pep_seq;
}

/**
 * @param {number|string} geneticCode
 */
function generateCodonTable(geneticCode) {
	const codon_table = {
		TTT: "F",
		TTC: "F",
		TTA: "L",
		TTG: "L",
		
		CTT: "L",
		CTC: "L",
		CTA: "L",
		CTG: "L",
		CTN: "L",
		
		ATT: "I",
		ATC: "I",
		ATA: "I",
		ATG: "M",
		
		GTT: "V",
		GTC: "V",
		GTA: "V",
		GTG: "V",
		GTN: "V",
		
		TCT: "S",
		TCC: "S",
		TCA: "S",
		TCG: "S",
		TCN: "S",
		
		CCT: "P",
		CCC: "P",
		CCA: "P",
		CCG: "P",
		CCN: "P",
		
		ACT: "T",
		ACC: "T",
		ACA: "T",
		ACG: "T",
		ACN: "T",
		
		GCT: "A",
		GCC: "A",
		GCA: "A",
		GCG: "A",
		GCN: "A",
		
		TAT: "Y",
		TAC: "Y",
		TAA: "*",
		TAG: "*",
		
		CAT: "H",
		CAC: "H",
		CAA: "Q",
		CAG: "Q",
		
		AAT: "N",
		AAC: "N",
		AAA: "K",
		AAG: "K",
		
		GAT: "D",
		GAC: "D",
		GAA: "E",
		GAG: "E",
		
		TGT: "C",
		TGC: "C",
		TGA: "*",
		TGG: "W",
		
		CGT: "R",
		CGC: "R",
		CGA: "R",
		CGG: "R",
		CGN: "R",
		
		AGT: "S",
		AGC: "S",
		AGA: "R",
		AGG: "R",
		
		GGT: "G",
		GGC: "G",
		GGA: "G",
		GGG: "G",
		GGN: "G",
	};

	if (geneticCode == 1) {
	}
	else if (geneticCode == 6) { // NCBI code: 6
		codon_table["TAA"] = "Q";
		codon_table["TAG"] = "Q";
	}
	else if (geneticCode == 4) { // NCBI code: 4
		codon_table["TGA"] = "W"; // not *
	}
	else if (geneticCode == 10) {// NCBI code: 10
		codon_table["TGA"] = "C"; // not *
	}
	else if (geneticCode == 27) { // NCBI code: 27
		codon_table["TAA"] = "Q"; // not *
		codon_table["TAG"] = "Q"; // not *
		codon_table["TGA"] = "W"; // or *
	}
	else if (geneticCode == 28) { // NCBI code: 28
		codon_table["TAA"] = "Q"; // or *
		codon_table["TAG"] = "Q"; // or *
		codon_table["TGA"] = "W"; // or *
	}
	else if (geneticCode == 29) { // NCBI code: 29
		codon_table["TAA"] = "Y"; // not *
		codon_table["TAG"] = "Y"; // not *
	}
	else {
		throw new Error(`Unknow code: ${geneticCode}`);
	}

	return codon_table;
}

function generateAATable_sp() {
	const list = {
		"U": {
			aa: "U",
			alias: [],
		},
		"O": {
			aa: "O",
			alias: [],
		},
		"X": {
			// name: "Xaa",
			// name: "Unk",
			aa: "X",
			alias: undefined,// any, unknow
		},
	};
	return list;
}
function generateAATable_co() {
	const list = {
		"B": {
			// name: "Asx",
			aa: "B",
			alias: [
				"D", "N",
			],
		},
		"Z": {
			// name: "Glx",
			aa: "Z",
			alias: [
				"E", "Q",
			],
		},
		"J": {
			// name: "Xle",
			aa: "J",
			alias: [
				"L", "I",
			],
		},
	};
	return list;
}

/**
 * @param {string} path_to_cds_fa
 * @param {number|string} geneticCode
 * @param {null|string} [path_to_protein_fa]
 */
function fasta_codonUsage(path_to_cds_fa, geneticCode, path_to_protein_fa = null, print_log = false) {
	const cds = readFasta(path_to_cds_fa);
	const protein = path_to_protein_fa ? readFasta(path_to_protein_fa) : null;
	return _fasta_codonUsage(cds, geneticCode, protein, print_log);
}

function cmd_main() {
	const argv = new ArgvParser(process.argv);
	const path_to_cds_fa = argv.get(/--cds=(.+)/, (arg, args) => args[0] ? args[0] : null);
	const geneticCode = argv.get(/--genetic=(.+)/, (arg, args) => args[0] ? args[0] : null);
	const path_to_protein_fa = argv.get(/--protein=(.+)/, (arg, args) => args[0] ? args[0] : null);
	const verbose = !!argv.get(/--verbose/, (arg, args) => true);
	
	if (path_to_cds_fa != null && geneticCode != null) {
		const result = fasta_codonUsage(path_to_cds_fa, geneticCode, path_to_protein_fa, verbose);
		console.log(Object.entries(result.codonUsage_dna_p_map).map(kv => kv.join("\t")).join("\n"));
	}
	else {
		console.error("Usage:", "AS-codon-usage --genetic=1 --cds=test.cds.fa");
		console.error("\t--genetic: NCBI genetic code");
		console.error("\t--cds: cds fasta file");
	}
}
if (__filename == fs.realpathSync(process.argv[1])) {
	cmd_main();
}

module.exports.generateAATable_sp = generateAATable_sp;
module.exports.generateAATable_co = generateAATable_co;
module.exports.generateCodonTable = generateCodonTable;
module.exports.cds_translate = cds_translate;

module.exports.fasta_codonUsage = fasta_codonUsage;

