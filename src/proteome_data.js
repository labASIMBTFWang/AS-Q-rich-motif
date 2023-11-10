// @ts-check

const fs = require("fs");
const zlib = require("zlib");
const path = require("path");
const https = require("https");

const bigJSON = require("big-json");

const { execAsync, spawnAsync, } = require("./util.js");
const { parseGFF, GFF_ROW } = require("./gff_util.js");
const { readFasta, saveFastaAsync, } = require("./fasta_util.js");
const { parse_InterPro_tsv } = require("./InterProScan_util.js");


/**
 * @abstract
 */
 class ProteinInfoTemplate {
	/**
	 * @protected
	 * @type {{ name: string, ordered_locus: string, orf: string, }[]}
	 */
	__gene_list;
	/**
	 * @protected
	 * @type {string}
	 */
	_sequence;
	/**
	 * @protected
	 * @type {{ id: string, term: string, }[]}
	 */
	_go = [];
	/**
	 * @protected
	 * @type {string}
	 */
	_protein_name;
}

/**
 * @typedef ProtCollection
 * @type {{ [entry_name: string]: ProtInfo; }}
 */


/**
 * @class
 * @example ProtInfo.loadFile("a.fa", "a.gff")
 */
 class ProtInfo extends ProteinInfoTemplate {
	/**
	 * @param {string} id
	 */
	constructor(id) {
		super();

		/**
		 * entry id
		 * @readonly
		 */
		this.id = id;

		/** @readonly */
		this._description = null;

		/** @protected */
		this._protein_name = null;
	}

	// /**
	//  * @param {fasta_util.FastaCollection} fa
	//  * @param {GFF_ROW[]} gff
	//  */
	// static load(fa, gff) {
	// }

	/**
	 * @param {...string} array_path_to_fa
	 */
	static loadFastaFile(...array_path_to_fa) {
		/** @type {ProtCollection} */
		const colle = {};

		for (let path_to_fa of array_path_to_fa) {
			const fa = readFasta(path_to_fa);

			Object.entries(fa).forEach(([id, seq]) => {
				const prot = new ProtInfo(id);
				prot.__gene_list = [
					{
						name: null,
						ordered_locus: id,
						orf: id
					},
				];
				prot._sequence = seq.seq;
				prot._go = [];
				if (colle[id] != null) {
					console.log(path_to_fa, id);
				}
				colle[id] = prot;
			});
			if (array_path_to_fa.length > 1) {
				console.error(path_to_fa, Object.keys(fa).length);
			}
		}

		if (array_path_to_fa.length > 1) {
			console.error(array_path_to_fa.join(","), Object.keys(colle).length);
		}

		return colle;
	}

	/**
	 * @TODO validation faa and interproscan has same seq name
	 * @param {string} path_to_fa
	 * @param {string} path_to_interproscan
	 */
	static loadInterProScan(path_to_fa, path_to_interproscan) {
		/** @type {ProtCollection} */
		const colle = {};
		
		const fa = readFasta(path_to_fa);
		
		for (let [id, vals] of parse_InterPro_tsv(path_to_interproscan)) {
			const prot = new ProtInfo(id);

			if (fa[id] == null) {
				console.log("not found seq:", `'${id}'`, vals, `'${path_to_fa}'`);
			}
			prot._sequence = fa[id].seq;
			prot.__gene_list = [
				{
					name: id,
					ordered_locus: "",
					orf: "",
				},
			];
			prot._go = Array.from(new Set(vals.map(v => v.GO_annotations).flat(1))).map(GO_annotation => {
				return {
					id: GO_annotation,
					term: null,
				};
			});
			prot._protein_name = null;

			colle[id] = prot;
		}

		return colle;
	}

	/**
	 * @param {string} path_to_fa
	 * @param {string} path_to_gff3
	 */
	static loadGffFile(path_to_fa, path_to_gff3) {
		/** @type {ProtCollection} */
		const colle = {};

		const fa = readFasta(path_to_fa);

		const gff = parseGFF(fs.readFileSync(path_to_gff3).toString());
		
		// const gene_list = Object.values(gff).flat().filter(a => a.type == "gene");// no chr
		const gene_list = Object.values(gff).flat().filter(a => a.type == "mRNA");// no chr // gff3

		// const label = {
		// 	biological_process: "P",
		// 	cellular_component: "C",
		// 	molecular_function: "F",
		// };

		/** @type {{ [seqId: string]: GFF_ROW; }} */
		const gff_gene_map = Object.fromEntries(gene_list.map(a => {
			return [a.attributes.ID, a];
		}));

		Object.entries(fa).forEach(([id, seq]) => {
			const prot = new ProtInfo(id);
			const gff_row = gff_gene_map[id];

			// @ts-ignore
			prot._raw = gff_row;

			if (gff_row == null) {
				console.log("Error no ID:", id, path_to_fa, path_to_gff3);
			}
			if (gff_row.attributes == null) {
				console.log("Error no attributes:", gff_row, path_to_fa, path_to_gff3);
			}
			if (id != gff_row.attributes.ID) {
				console.error(id, gff_row.attributes.ID, gff_row);
			}

			prot.__gene_list = [
				{
					name: null,
					ordered_locus: id,
					orf: id
				},
			];
			prot._sequence = seq.seq;

			prot._go = gff_row.attributes.Ontology_term ? gff_row.attributes.Ontology_term.match(/GO:\d+/g).map(go_id => {
				return {
					id: go_id,
					term: null,
				};
			}) : [];
			prot._go.forEach(go_id => {
				if (!go_id.id.startsWith("GO:")) {
					console.error(id, gff_row);
					debugger;
				}
			});
			
			prot._protein_name = gff_row.attributes.product;

			colle[id] = prot;
		});

		return colle;
	}

	/** @param {string} v */
	set sequence(v) {
		/** @protected */
		this._sequence = v;
	}
	get sequence() {
		return this._sequence;
	}

	/**
	 * @readonly
	 */
	get protein_name() {
		return this._protein_name;
	}

	/**
	 * @param {{ id: string; term: string; }[]} gg
	 */
	set GO(gg) {
		/** @protected */
		this._go = gg;
	}
	get GO() {
		return this._go;
	}

	/** @param {{ name: string, ordered_locus: string; orf: string }[]} v*/
	set gene_list(v) {
		/** @protected */
		this.__gene_list = v;
	}
	get gene_list() {
		return this.__gene_list;
	}
}

class UniprotInfo extends ProtInfo {
	/**
	 * @param {string} id
	 * @param {any} entry
	 */
	constructor(id, entry) {
		super(id);

		/** @readonly */
		this.entry = entry;
	}

	static symbol_accession_map = Symbol("symbol_accession_map");
	static symbol_EntryCollection_info = Symbol("UniprotEntryCollection");

	/**
	 * @param {UniprotEntryCollection} colle
	 * @returns {UniprotEntryCollection}
	 */
	static getAccessionMap(colle) {
		return /** @type {any} */(colle[UniprotInfo.symbol_accession_map]);
	}

	/** @type {string} */
	get sequence() {
		return this.entry.sequence["#text"];
	}

	get GO() {
		if (this._GO == null) {
			/**
			 * @protected
			 */
			this._GO = this._get_GO();
		}
		return this._GO;
	}
	/**
	 * @protected
	 */
	_get_GO() {
		const entry = this.entry;
		/** @type {any[]} */
		const dbReference = Array.isArray(entry.dbReference) ? entry.dbReference : [entry.dbReference];

		if (dbReference) {
			return dbReference.filter(a => a["@type"] == "GO").map(gg => {
				return {
					id: gg["@id"],
					term: gg.property.find(a => a["@type"] == "term")["@value"],// ["P:", "C:", "F:"]
				};
			});
		}
		else {
			return [];
		}
	}

	getAlphaFold() {
		const entry = this.entry;
		/** @type {any[]} */
		const dbReference = Array.isArray(entry.dbReference) ? entry.dbReference : [entry.dbReference];
		if (dbReference) {
			return dbReference.filter(a => a["@type"] == "AlphaFoldDB").map(gg => {
				return gg["@id"];
			});
		}
		else {
			return [];
		}
	}

	/** @type {string} */
	get protein_name() {
		const recommendedName = this.entry.protein.recommendedName;
		const submittedName = this.entry.protein.submittedName;
		return (recommendedName ?? submittedName)?.fullName?.["#text"] ?? "";
	}
	
	get gene_list() {
		if (this._gene_list == null) {
			/**
			 * @protected
			 */
			this._gene_list = this._get_ordered_locus_name_ID();
		}
		return this._gene_list;
	}
	/**
	 * @protected
	 */
	_get_ordered_locus_name_ID() {
		const entry = this.entry;
		/** @type {any[]} */
		const gene_list = entry.gene ? (Array.isArray(entry.gene) ? entry.gene : [entry.gene]).filter(a => a) : [];
		if (gene_list.length > 0) {
			try {
				return gene_list.map(gene => {
					const name_list = Array.isArray(gene.name) ? gene.name : [gene.name];
					const ol = name_list.find(v => v["@type"] == "ordered locus");
					const primary = name_list.find(v => v["@type"] == "primary");
					const orf = name_list.find(v => v["@type"] == "ORF");
					if (primary == null && ol == null && orf == null) {
						console.log(this);
						console.log(gene_list);
						debugger;
					}
					const primary_name = primary?.["#text"];
					const ordered_locus = ol?.["#text"];
					const orf_name = orf?.["#text"];

					if (primary_name || ordered_locus) {
						if (orf_name == null) {
							// console.error("not found ORF ID:", primary_name, ordered_locus);
							debugger;
						}
					}

					return {
						/** gene name */
						name: primary_name,

						/** locus */
						ordered_locus: ordered_locus ?? orf_name ?? primary_name,
						
						orf: orf_name ?? ordered_locus ?? primary_name,
					};
				});
			}
			catch (ex) {
				console.log(this);
				console.log(gene_list);
				debugger
				return [
					{
						ordered_locus: entry.name,
						name: entry.name,
						orf: entry.name,
					},
				];
			}
		}
		else if (entry.protein && entry.protein.recommendedName) {
			return [
				{
					ordered_locus: entry.name,
					name: entry.protein.recommendedName.fullName ?? entry.protein.recommendedName.shortName,
					orf: entry.protein.recommendedName.fullName ?? entry.protein.recommendedName.shortName,
				},
			];
		}
		else {
			return [
				{
					ordered_locus: entry.name,
					name: entry.name,
					orf: entry.name,
				},
			];
		}
	}
	
	/**
	 * @param {string} uniprot_proteome_id
	 */
	static async download_proteome_xml(uniprot_proteome_id) {
		const fileURL = new URL("/uniprotkb/stream", "https://rest.uniprot.org");
		// fileURL.searchParams.set("compressed", "true");           // default: true
		fileURL.searchParams.set("download", "true");                // ???
		// fileURL.searchParams.set("format", "xml");                // defaultt: JSON
		fileURL.searchParams.set("query", `(proteome:${uniprot_proteome_id})`); //
		
		const xml_filename = path.resolve(`uniprot-proteome_${uniprot_proteome_id}.xml`);
		const xml_json_filename = path.resolve(`uniprot-proteome_${uniprot_proteome_id}.json`);
	
		if (!fs.existsSync(xml_filename) || fs.statSync(xml_filename).size <= 0) {
			console.log("begin download proteome:", xml_filename, fileURL);
			
			const promise = new Promise(function (resolve, reject) {
				const wstream = fs.createWriteStream(xml_filename);
				const gun = zlib.createGunzip();
				const gun2 = zlib.createGunzip();
				https.get(fileURL, {
					headers: {
						"Accept": "application/xml",
						// URL query param "compressed" default is true
						"Accept-Encoding": "gzip", // ?compressed=true
					},
				}, function (rsp) {
					console.log(fileURL, rsp.headers);

					rsp.pipe(gun).pipe(gun2).pipe(wstream);
					wstream.on("close", function onClose() {
						resolve(xml_filename);
					});
					rsp.on("error", function onRecvResponseError(err) {
						reject(err);
					});
				}).on("error", function onSendRequestError(err) {
					reject(err);
				});
			});
			try {
				await promise;
				console.log("downloaded proteome:", path.resolve(xml_filename), fileURL);
			}
			catch (ex) {
				console.log("error download proteome:", xml_filename, fileURL, ex);
			}
		}
		return await info_XML_to_JSON(xml_filename);
	}
}

/**
 * @typedef UniprotEntryCollection
 * @type {{ [entry_name: string]: UniprotInfo; }}
 */

/**
 * @param {string} uniprot_proteome_id
 */
async function load_proteome_xml_json(uniprot_proteome_id) {
	const path_to_json = await UniprotInfo.download_proteome_xml(uniprot_proteome_id);

	const stat = fs.statSync(path_to_json);
	const uniprot_objet = stat.size > 200_000_000 ?
		await JSONParseAsync(path_to_json)
		:
		JSON.parse(fs.readFileSync(path_to_json).toString())
	;

	/** @type {Record<string, any>[]} */
	const list = uniprot_objet?.uniprot?.entry;

	const entries = list.map(entry => {
		return [
			entry.name,
			new UniprotInfo(entry.name, entry),
		];
	});
	
	const path_to_fa = `${uniprot_proteome_id}.fa`;
	if (!fs.existsSync(path_to_fa)) {
		const fa = Object.fromEntries(entries.map(([k, pp]) => {
			return [k, pp.sequence];
		}));
		const promise = await saveFastaAsync(path_to_fa, fa);
	}

	return Object.fromEntries(entries);
}

/**
 * @param {string} proteins_fa
 * @param {string} gff3
 */
function load_gff_as_proteome(proteins_fa, gff3) {
	return ProtInfo.loadGffFile(proteins_fa, gff3);
}

/**
 * @param {string} proteins_fa
 * @param {string} interproscan_tsv
 */
function load_InterProScan_as_proteome(proteins_fa, interproscan_tsv) {
	return ProtInfo.loadInterProScan(proteins_fa, interproscan_tsv);
}

/**
 * @param {string} path_to_file
 */
function JSONParseAsync(path_to_file) {
	return new Promise((resolve, reject) => {
		const parseStream = bigJSON.createParseStream();
		
		const readStream = fs.createReadStream(path_to_file);
		
		parseStream.on("error", reject);

		parseStream.on("data", function (pojo) {
			resolve(pojo);
		});
		 
		// @ts-ignore
		readStream.pipe(parseStream);
	});
}

/**
 * pip install yq
 * wget https://github.com/jqlang/jq/releases/download/jq-1.7/jq-linux-amd64 -O ~/bin/jq
 * chmod +x ~/bin/jq
 * @param {string} path_to_xml
 */
async function info_XML_to_JSON(path_to_xml) {
	const path_to_json = path_to_xml.replace(/\.xml$/, ".json");
	if (!fs.existsSync(path_to_json) || fs.statSync(path_to_json).size <= 0) {
		const cmd = `cat ${path_to_xml} | xq . > ${path_to_json}`;
		console.error(cmd);
		await execAsync(cmd).promise;
	}
	else {
		console.error("found JSON:", path_to_json);
	}
	return path_to_json;
}

module.exports.ProteinInfoTemplate = ProteinInfoTemplate;
module.exports.ProtInfo = ProtInfo;
module.exports.UniprotInfo = UniprotInfo;

module.exports.load_proteome_xml_json = load_proteome_xml_json;
module.exports.load_gff_as_proteome = load_gff_as_proteome;
module.exports.load_InterProScan_as_proteome = load_InterProScan_as_proteome;

