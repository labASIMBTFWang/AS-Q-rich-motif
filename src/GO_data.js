// @ts-check
// version: 2023_11_10

const fs = require("fs");

// https://www.npmjs.com/package/fast-xml-parser
const { XMLParser, XMLBuilder, XMLValidator} = require("fast-xml-parser");


class GO_data {
	/**
	 * @protected
	 * @type {{ [GO: string]: GO_data; }}
	 */
	_children;
	
	constructor(go) {
		/** @protected */
		this.go = go;
	}
	
	/** @type {string} */
	get id() { return this.go["oboInOwl:id"]; };
	/** @type {string} */
	get domain() { return this.go["oboInOwl:hasOBONamespace"]; };
	/** @type {string} */
	get name() { return this.go["rdfs:label"]?.["#text"] ?? this.go["rdfs:label"]; };

	/**
	 * @readonly
	 * @type {GO_data[] | null}
	 */
	get children() {
		if (this._children) {
			return Object.values(this._children);
		}
		else {
			return null;
		}
	}

	get parent() {
		if (this.go["rdfs:subClassOf"] != null) {
			if (this._parent == null) {
				const that = this;
		
				/** @type {any[]} */
				const p_list = Array.isArray(this.go["rdfs:subClassOf"]) ? this.go["rdfs:subClassOf"] : [this.go["rdfs:subClassOf"]];
		
				/** @protected */
				this._parent = p_list.map(p => {
					const res = p["@rdf:resource"];
					if (res) {
						const parent_id = "GO:" + res.match(/\/GO_(\d+)/)[1];
						const parent = GO_data.map[parent_id];
						
						parent._children ??= {};
						parent._children[that.id] = that;
						return parent;
					}
				}).filter(a => a != null);
			}
			return this._parent;
		}
		else {
			return null;
		}
	}

	get level() {
		if (this._level == null) {
			if (this.parent) {
				const lv_list = this.parent.map(a => a.level + 1);
				this._level = Math.min(...lv_list);
			}
			else {
				return GO_data.LEVEL_ROOT;
			}
		}

		return this._level;
	}

	/**
	 * @returns {this|GO_data[]}
	 */
	get root() {
		if (this.parent == null) {
			return this;
		}
		else {
			return this.parent.map(a => a.root).flat(1).filter(a => a.isRoot());
		}
	}

	isRoot() {
		return this.root == this;
	}

	toJSON() {
		return {
			id: this.id,
			domain: this.domain,
			name: this.name,
			leve: this.level,
			parent: Array.isArray(this.parent) ? Object.values(this.parent) : [],
			children: Array.isArray(this.children) ? Object.values(this.children) : [],
		}
	}

	/** @readonly */
	static LEVEL_ROOT = 1;

	/** @type {{ [GO_id: string]: GO_data; }} */
	static map;
	/** @type {{ [GO_term: string]: GO_data; }} */
	static go_term_map;
}

/**
 * @returns {{ [GO_id: string]: GO_data; }}
 * @see {@link http://purl.obolibrary.org/obo/go.owl}
 */
function get_GO_data_map() {
	const options = {
		attributeNamePrefix : "@",
		ignoreAttributes: false,
		tagValueProcessor: (tagName, tagValue, jPath, hasAttributes, isLeafNode) => {
			if (isLeafNode) {
				return tagValue;
			}
			return "";
		},
	};
	const parser = new XMLParser(options);
	const jObj = parser.parse(fs.readFileSync("./go.owl"));
	
	/** @type {{ [GO_id: string]: GO_data; }} */
	const go_map = {
	};
	/** @type {{ [GO_term: string]: GO_data; }} */
	const go_term_map = {
	};

	jObj["rdf:RDF"]["owl:Class"].forEach(go => {
		const go_id = go["oboInOwl:id"];

		const gdata = new GO_data(go);

		go_map[go_id] = gdata;

		"oboInOwl:hasExactSynonym"
		"oboInOwl:hasBroadSynonym"
		if (go["oboInOwl:hasNarrowSynonym"]) {
			const hasNarrowSynonym_list = Array.isArray(go["oboInOwl:hasNarrowSynonym"]) ? go["oboInOwl:hasNarrowSynonym"] : [go["oboInOwl:hasNarrowSynonym"]];
			go_term_map[gdata.name] = gdata;

			hasNarrowSynonym_list.forEach((alt_go_term, i) => {
				go_term_map[alt_go_term["#text"] ?? alt_go_term] = gdata;
			});
		}

		const alternative_list = Array.isArray(go["oboInOwl:hasAlternativeId"]) ? go["oboInOwl:hasAlternativeId"] : [go["oboInOwl:hasAlternativeId"]];
		alternative_list.forEach(alt_go_id => {
			go_map[alt_go_id] = gdata;
		});
	});

	GO_data.map = go_map;
	GO_data.go_term_map = go_term_map;

	return go_map;
}


module.exports.GO_data = GO_data;
module.exports.get_GO_data_map = get_GO_data_map;
