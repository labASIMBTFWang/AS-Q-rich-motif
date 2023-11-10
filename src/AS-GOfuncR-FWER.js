// @ts-check

const fs = require("fs");
const path = require("path");

const { ArgvParser, } = require("./ArgvParser.js");
const { spawnAsync, } = require("./util.js");
const { readFasta, } = require("./fasta_util.js");
const { find_SCD_v4, MAX_DISTANCE, } = require("./find_SCD.js");
const { findPolyXX, } = require("./find_polyX.js");
const {
	load_proteome_xml_json,
	load_gff_as_proteome,
	load_InterProScan_as_proteome,
} = require("./proteome_data.js");

/**
 * ["Q", "N"].join("⁄") => Q⁄N
 */
const all_poly_task = [
	"Q", ["Q", "N"], "N",
	"Y", "W", "C", "P", "A",
	"S", "T",
	"D", "E", "G", "F", "H", "I", "L", "V", "M", "R", "K",
];

/**
 * 
 * @param {string} root_dir
 * @param {string} species_tag
 * @param {(string | string[])[]} poly_task
 * @param {string} path_to_fa
 */
async function GOfuncR(root_dir, species_tag, poly_task, path_to_fa) {
	const subtypes = [
		"SCD",
		...poly_task.map(a => Array.isArray(a) ? a.join("⁄") : a).map(v => `poly${v}`),
	];

	const fa = readFasta(path_to_fa);
	if (!fa) {
		throw new Error("required fasta");
	}

	{
		const scd_genes = Object.keys(fa).filter(key => {
			const seq = fa[key];
			const scd_list = find_SCD_v4(seq, MAX_DISTANCE);
			return scd_list.length > 0;
		});
		// console.log("scd_genes:", info.species_tag, scd_genes.length);
		out_tsv("SCD", scd_genes);
	}

	poly_task.map(aa => {
		const genes = Object.keys(fa).filter(key => {
			const seq = fa[key];
			const regions = findPolyXX(seq, 4, aa);
			return regions.length;
		});
		out_tsv(`poly${(Array.isArray(aa) ? aa.join("⁄") : aa)}`, genes);
	});

	/**
	 * @param {string} subtype
	 * @param {string[]} genes
	 */
	function out_tsv(subtype, genes) {
		const header_cols = ["gene_ids", "is_candidate"];
		const ffdd = path.join("./", "GO_Enrichment", species_tag);
		const ffpp = path.join(ffdd, subtype + ".tsv");
		if (!fs.existsSync(ffdd)) {
			fs.mkdirSync(ffdd, { recursive: true, });
		}
		fs.writeFileSync(ffpp, header_cols.join("\t") + "\n" + genes.map(v => [v, 1].join("\t")).join("\n"));
		console.log("*ffpp:", genes.length, ffpp);
	}

	const src = get_src_GO_Enrich("./", subtypes, species_tag);
	const stats_r = path.resolve(root_dir, species_tag, "run_GOfuncR.R");
	fs.writeFileSync(stats_r, src);
	
	const proc = spawnAsync("R", [
		"--quiet",
		"--no-save",
	]);
	console.log("[t:" + species_tag + "]", "run");

	if (proc.stdin) {
		proc.stdin.write(src);
		proc.stdin.end(); // close
		const ffdd = path.join("./", "GO_Enrichment", species_tag);
		const log_out = path.resolve(path.join(ffdd, "run.stdout.txt"));
		const log_err = path.resolve(path.join(ffdd, "run.stderr.txt"));
		proc.stdoutToFile(log_out);
		proc.stderrToFile(log_err);
		await proc.promise;
		console.log("[t:" + species_tag + "] exit:", proc.exitCode, proc.signalCode, log_err, log_out);
	}
	else {
		console.log("[t:" + species_tag + "] no run R");
	}
}


/**
 * 
 * @param {string} root_dir
 * @param {string[]} subtypes
 * @param {string} species_tag
 */
function get_src_GO_Enrich(root_dir, subtypes, species_tag) {
return String.raw`
library(GOfuncR)
library(openxlsx)
library(dplyr)

root_dir <- "${root_dir}"

poly_list <- c(${subtypes.map(v => `"${v}"`).join(",")})

make_polyAA_GO_enrich <- function (species, bSave = FALSE) {
	path_to_bg <- paste0(root_dir, "annotations/", species, "_annotations.tsv")
	bg <- read.csv(path_to_bg, header = T, sep = "\t")
	cat("path_to_bg:", path_to_bg, "\t", nrow(bg), "\n")
	
	wb = createWorkbook()

	for (polyType in poly_list) {
		path_to_fg <- paste0(root_dir, species, "/", polyType, ".tsv")
		fg <- read.csv(path_to_fg, header = T, sep = "\t")
		cat("path_to_fg:", path_to_fg, "\t", nrow(fg), "\n")

		done = FALSE
		path_to_stats_tsv <- paste0(root_dir, species, "/", polyType, "_stats.tsv")

		if (bSave & nrow(fg) > 0 & nrow(bg) > 0) {
			bg_fg_common_genes <- inner_join(data.frame(gene=fg$gene_ids, is_candidate=fg$is_candidate), bg)
			if (nrow(bg_fg_common_genes) > 0) {
				res_hyper_anno <- go_enrich(fg, annotations=bg)
				stats <- res_hyper_anno[[1]]

				colnames(stats) <- c(
					"ontology",
					"GO ID",
					"GO term",
					"raw p-value under-representation",
					"raw p-value over-representation",
					"FWER under-representation",
					"FWER over-representation"
				)
				
				count_polyQ_group_by_GO <- length(intersect(fg[[1]], bg[[1]]))

				write.table(stats, path_to_stats_tsv, sep = "\t", row.names=FALSE, quote=F)
				cat("path_to_stats_tsv:", path_to_stats_tsv, "\n")

				addWorksheet(wb, polyType)
				writeData(wb, sheet=polyType, x=stats, rowNames=FALSE)

				done = TRUE
			}
		}
		if (bSave & !done) {
			# make dummy
			stats <- data.frame(ontology=c(""), GO.ID=c(""), GO.term=c(""), raw.p.u=c(""), raw.p.o=c(""), FWER.u=c(""), FWER.o=c(""))
			colnames(stats) <- c(
				"ontology",
				"GO ID",
				"GO term",
				"raw p-value under-representation",
				"raw p-value over-representation",
				"FWER under-representation",
				"FWER over-representation"
			)
			addWorksheet(wb, polyType)
			writeData(wb, sheet=polyType, x=stats, rowNames=FALSE)
			
			write.table(stats, path_to_stats_tsv, sep = "\t", row.names=FALSE, quote=F)
			cat("path_to_stats_tsv <empty>:", path_to_stats_tsv, "\n")
		}
	}
	if (bSave) {
		path_to_stats_xlsx <- paste0(root_dir, species, "_stats.xlsx")
		saveWorkbook(wb, path_to_stats_xlsx)
		cat("path_to_stats_xlsx:", path_to_stats_xlsx, "\n")
	}
	return(wb)
}
make_polyAA_GO_enrich("${species_tag}", TRUE)
`;
}

// UP000230249
async function cmd_main() {
	const argv = new ArgvParser(process.argv);
	const prefix = argv.get(/--prefix=(.+)/, (arg, args) => args[0] ? args[0] : null);
	const uniprot_proteome_id = argv.get(/--uniprot-proteome-id=(.+)/, (arg, args) => args[0] ? args[0] : null);
	const path_to_gff = argv.get(/--gff=(.+)/, (arg, args) => args[0] ? args[0] : null);
	const path_to_interproscan = argv.get(/--interproscan=(.+)/, (arg, args) => args[0] ? args[0] : null);
	let path_to_fa = argv.get(/--fa=(.+)/, (arg, args) => args[0] ? args[0] : null);
	
	if (!prefix) {
		console.error("required argument: --prefix=<output prefix>");
		return;
	}

	if (uniprot_proteome_id) {
		path_to_fa = `${uniprot_proteome_id}.fa`;
	}

	const proteome = await (async function () {
		if (uniprot_proteome_id) {
			return await load_proteome_xml_json(uniprot_proteome_id);
		}
		else if (
			path_to_fa && fs.existsSync(path_to_fa) &&
			path_to_gff && fs.existsSync(path_to_gff)
		) {
			return await load_gff_as_proteome(path_to_fa, path_to_gff);
		}
		else if (
			path_to_fa && fs.existsSync(path_to_fa) &&
			path_to_interproscan && fs.existsSync(path_to_interproscan)
		) {
			return await load_InterProScan_as_proteome(path_to_fa, path_to_interproscan);
		}
		else {
			console.error(`usage: AS-GOfuncR-FWER --uniprot-proteome=UP000230249`);
			console.error(`usage: AS-GOfuncR-FWER --gff=annotation.gff`);
			console.error(`usage: AS-GOfuncR-FWER --interproscan=interproscan.output.tsv`);
			console.error(`	--uniprot-proteome=<ID|json>`);
			console.error(`	--gff=annotation.gff`);
			console.error(`	--interproscan=interproscan.output.tsv`);
		}
	})();

	if (proteome != null && path_to_fa != null) {
		const path_to_anno_tsv = `${prefix}_annotations.tsv`;

		await fs.promises.writeFile(path_to_anno_tsv, [
			["gene", "go_id"].join("\t"),
			Object.keys(proteome).map(entId => {
				const entry = proteome[entId];
				return entry.GO.map(v => [entry.id, v.id].join("\t"));
			}).flat(1).join("\n"),
		].join("\n")); // annotations/

		await GOfuncR("./", prefix, all_poly_task, path_to_fa);

		[
			"SCD",
			...all_poly_task.map(a => Array.isArray(a) ? a.join("⁄") : a).map(v => `poly${v}`),
		].forEach(polyType => {
			const path_to_stats_tsv = path.resolve("./", prefix, "/", polyType, "_stats.tsv");
			console.log("output:", path_to_stats_tsv);
		})
	}
}
if (__filename == fs.realpathSync(process.argv[1])) {
	cmd_main();
}

