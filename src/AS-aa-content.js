// @ts-check

const fs = require("fs");

const { ArgvParser, } = require("./ArgvParser.js");
const { readFasta, } = require("./fasta_util.js");

function cmd_main() {
	const argv = new ArgvParser(process.argv);
	const path_to_fa = argv.get(/--fa=(.+)/, (arg, args) => args[0] ? args[0] : null);
	if (path_to_fa) {
		const fa = readFasta(path_to_fa);

		/** @type {Record<string, number>} */
		const content = {
			"A": 0, "R": 0, "N": 0, "D": 0, "C": 0,
			"Q": 0, "E": 0, "G": 0, "H": 0, "I": 0,
			"L": 0, "K": 0, "M": 0, "F": 0, "P": 0,
			"S": 0, "T": 0, "W": 0, "Y": 0, "V": 0,
			"*": 0,
			"U": 0, "O": 0,
		};
		let total_sum = 0;
		Object.values(fa).forEach(seq => {
			for (let i = 0; i < seq.length; ++i) {
				const aa = seq[i];
				content[aa] = (content[aa] ?? 0) + 1;
				total_sum += 1;
			}
		});
		const tsv = Object.keys(content).map(aa => [aa, content[aa] / total_sum].join("\t")).join("\n");
		console.log(tsv);
	}
	else {
		console.error(`AS-aa-content --fa=test.protein.fa`);
	}
}

if (__filename == fs.realpathSync(process.argv[1])) {
	cmd_main();
}

