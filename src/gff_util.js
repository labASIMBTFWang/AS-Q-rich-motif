// @ts-check

const fs = require("fs");

const { readFasta, } = require("./fasta_util.js");

const GFF_header = [
	"seqid", "source", "type", "start", "end", "score", "strand", "phase", "_attributes"
];

class GFF_ROW {
	/** @type {string} */
	seqid;
	/** @type {string} */
	type;
	/** @type {number} - raw pos */
	start;
	/** @type {number} - raw pos */
	end;
	/** @type {number} - raw pos */
	length;
	/** @type {number} */
	score;
	/** @type {number} */
	strand;
	/** @type {number} */
	phase;
	/** @type {string} */
	_attributes;
	/** @type {{ [attribute: string]: any; }} */
	$attributes;

	/**
	 * @param {Partial<GFF_ROW>} gff
	 */
	constructor(gff) {
		if (gff.seqid != null) this.seqid = gff.seqid;
		// if (gff.source != null) this.source = gff.source;
		if (gff.type != null) this.type = gff.type;
		if (gff.start != null) this.start = Number(gff.start);
		if (gff.end != null) this.end = Number(gff.end);
		this.length = this.end - this.start + 1;
		if (gff.score) this.score = gff.score;
		if (gff.strand) this.strand = gff.strand;
		if (gff.phase) this.phase = gff.phase;
		//this.attributes = gff.attributes;
		if (gff._attributes) this._attributes = gff._attributes;
		// this.$attributes = null;
	}

	get geneID() {
		const gene_sub_ID = this.attributes.ID;
		return gene_sub_ID.match(/^.*_\d+/)[0];
	}

	/**
	 * = *._\d+-T\d+
	 * gene isoform ID
	 */
	get Geneid() {
		const gene_sub_ID = this.attributes.ID;
		return gene_sub_ID.match(/^.*_\d+-T\d+/)[0];
	}

	// /**
	//  * @type {number}
	//  */
	// get end() {
	// 	return this.start + this.length - 1;
	// }

	/**
	 * lazy load
	 * @type {{[attribute:string]:any}}
	 */
	get attributes() {
		if (this.$attributes) {
			return this.$attributes;
		}
		else {
			const attributes = this._parse_gff_attributes(this._attributes);
			this.$attributes = attributes;
			return attributes;
		}
	}

	_parse_gff_attributes(text_attributes) {
		let attributes = {};
		//let regexp = /ID=(?<ID>[^;]+);/;
		// attributes.ID =            text_attributes.match(/ID=[^;]+/);
		// attributes.Name =          text_attributes.match(/Name=[^;]+/);
		// attributes.Alias =         text_attributes.match(/Alias=[^;]+/);
		// attributes.Parent =        text_attributes.match(/Parent=[^;]+/);
		// attributes.Target =        text_attributes.match(/Target=[^;]+/);
		// attributes.Gap =           text_attributes.match(/Gap=[^;]+/);
		// attributes.Derives_from =  text_attributes.match(/Derives_from=[^;]+/);
		// attributes.Note =          text_attributes.match(/Note=[^;]+/);
		// attributes.Dbxref =        text_attributes.match(/Dbxref=[^;]+/);
		// attributes.Ontology_term = text_attributes.match(/Ontology_term=[^;]+/);
	
		text_attributes.split(";").map(col => {
			let [key, value] = col.split("=");
			attributes[key] = value;
		});
	
		return attributes;
	}

	static get GFF_header() {
		return GFF_header;
	}
}

/**
 * @template T
 * @param {string} text_gff
 * @param {{ [seqid: string]: (T extends (number|string) ? T : never); }} [seqid_to_chr]
 * @returns {{ [chr: string]: GFF_ROW[]; }}
 */
function parseGFF(text_gff, seqid_to_chr) {
	let ann_rows = _parseTable(text_gff, GFF_header, 1);// skip first line, that startsWith #

	/** @type {{ [chrName:string]: GFF_ROW[] }} */
	let group_by_seq = {};

	ann_rows.forEach(row => {
		if (row.seqid.startsWith("#")) {
			// skip #
		}
		else {
			const gff_row = new GFF_ROW(row);
			
			gff_row.start = Number(row.start);
			gff_row.end = Number(row.end);
			gff_row.length = gff_row.end - gff_row.start + 1;
			gff_row.strand = row.strand == "-" ? -1 : (row.strand == "+" ? 1 : 0);
			//gff_row.attributes = parse_gff_attributes(row._attributes);

			const chr = seqid_to_chr ? seqid_to_chr[gff_row.seqid] : gff_row.seqid;

			if (!group_by_seq[chr]) {
				group_by_seq[chr] = [];
			}
			group_by_seq[chr].push(gff_row);
		}
	});

	return group_by_seq;
}

/**
 * @param {string} text
 * @param {string[]} prop_name_list
 * @param {number} [start_line]
 * @returns {{ [prop: string]: string }[]}
 */
function _parseTable(text, prop_name_list, start_line = 0) {
	const _a = text.trim().split("\n").slice(start_line).map(line => {
		/** @type {{ [prop:string]: string }} */
		const obj = {};
		line = line.trim();
		if (line) {
			return line.split("\t").reduce((obj, col, j) => {
				obj[prop_name_list[j]] = col.trim();
				return obj;
			}, obj);
		}
		else {
			return null;
		}
	}).filter(a => a);

	return /** @type {{ [prop: string]: string }[]} */(_a);
}

module.exports.parseGFF = parseGFF;
module.exports.GFF_ROW = GFF_ROW;
