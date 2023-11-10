// @ts-check

const fs = require("fs");

const { groupBy, } = require("./util.js");


/**
 * @param {string[]} cols
 */
function InterProScan_fromColumns(cols) {
	const obj = {};
	let str_seq_len;
	let str_start;
	let str_stop;
	let str_score;
	let str_date;
	let GO_annotations;

	obj.seq_name =              cols[ 0]; // Protein accession (e.g. P51587)
	obj.md5 =                   cols[ 1]; // Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
	str_seq_len =               cols[ 2]; // Sequence length (e.g. 3418)
	obj.program =               cols[ 3]; // Analysis (e.g. Pfam / PRINTS / Gene3D)
	obj.signature_accession =   cols[ 4]; // Signature accession (e.g. PF09103 / G3DSA:2.40.50.140)
	obj.signature_description = cols[ 5]; // Signature description (e.g. BRCA2 repeat profile)
	str_start =                 cols[ 6]; // Start location
	str_stop =                  cols[ 7]; // Stop location
	str_score =                 cols[ 8]; // Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
	obj.status =                cols[ 9]; // Status - is the status of the match (T: true)
	str_date =                  cols[10]; // Date - is the date of the run
	obj.InterPro_accession =    cols[11]; // InterPro annotations - accession (e.g. IPR002093)
	obj.InterPro_description =  cols[12]; // InterPro annotations - description (e.g. BRCA2 repeat)
	GO_annotations =            cols[13]; // (GO annotations (e.g. GO:0005515) - optional column; only displayed if –goterms option is switched on)
	obj.Pathways_annotations =  cols[14]; // (Pathways annotations (e.g. REACT_71) - optional column; only displayed if –pathways option is switched on)

	obj.seq_len = Number(str_seq_len) | 0;
	obj.start = Number(str_start) | 0;
	obj.stop = Number(str_stop) | 0;
	obj.str_score = Number(str_score) | 0;

	obj.GO_annotations = (!GO_annotations || GO_annotations == "-") ? [] : GO_annotations.split("|");

	return obj;
}

/**
 * @param {string} path_to_tsv
 */
function parse_InterPro_tsv(path_to_tsv) {
	const text = fs.readFileSync(path_to_tsv).toString();
	const tab = text.split("\n").filter(a => a.trim()).map(a => a.split("\t").map(b => b.trim()));
	// const tab = tsv_parse(text); // has some bug
	const list = tab.map(cols => InterProScan_fromColumns(cols));
	return groupBy(list, "seq_name");
}

module.exports.parse_InterPro_tsv = parse_InterPro_tsv;
