// @ts-check

const child_process = require("child_process");
const fs = require("fs");


class ProcessUtiil {
	/**
	 * exitPromise
	 * @type {Promise<child_process.ChildProcess>}
	 */
	promise;
	
	/**
	 * do not await proc.promise // process no write stdio if exit
	 */
	async stdoutToString() { return ""; }
	/**
	 * do not await proc.promise // process no write stdio if exit
	 */
	async stderrToString() { return ""; }
	
	/**
	 * @param {string} filename
	 * @returns {fs.WriteStream|null}
	 */
	stdoutToFile(filename) { return null; }
	/**
	 * @param {string} filename
	 * @returns {fs.WriteStream|null}
	 */
	stderrToFile(filename) { return null; }
	/**
	 * @param {string} out_filename
	 * @param {string} err_filename
	 * @returns {(fs.WriteStream|null)[]}
	 */
	stdioToFile(out_filename, err_filename) { return [null, null, null]; }

	/**
	 * @param {import("child_process").ChildProcess} proc
	 * @param {Promise<child_process.ChildProcess>} exit_promise
	 */
	static assign(proc, exit_promise) {
		const stdoutToString = async () => proc.stdout ? await streamToString(proc.stdout) : ""; // do not await proc.promise // process write stdio when live
		const stderrToString = async () => proc.stderr ? await streamToString(proc.stderr) : ""; // do not await proc.promise // process write stdio when live
		/** @param {string} filename */
		const stdoutToFile = (filename) => proc.stdout ? proc.stdout.pipe(fs.createWriteStream(filename)) : null;
		/** @param {string} filename */
		const stderrToFile = (filename) => proc.stderr ? proc.stderr.pipe(fs.createWriteStream(filename)) : null;
		/**
		 * @param {string} out_filename
		 * @param {string} err_filename
		 */
		const stdioToFile = (out_filename, err_filename) => {
			return [
				null,
				stdoutToFile(out_filename),
				stderrToFile(err_filename),
			];
		};

		/**
		 * @type {ProcessUtiil}
		 */
		const uu = {
			/** exit_promise */
			promise: exit_promise,
			stdoutToString: stdoutToString,
			stderrToString: stderrToString,
			stdoutToFile: stdoutToFile,
			stderrToFile: stderrToFile,
			stdioToFile: stdioToFile,
			// then: promise.then,
			// catch: promise.catch,
		};
		
		/** @type {child_process.ChildProcess & ProcessUtiil} */
		const v = Object.assign(proc, uu)
		return v;
	}
}

execAsync.log = "log";
/**
 * @param {string} cmd
 * @param {child_process.ExecOptions} [options]
 * returns {child_process.ChildProcess & { promise: Promise<child_process.ChildProcess>; then: typeof Promise.prototype.then; catch: typeof Promise.prototype.catch; stdoutToString: () => Promise<string>; stderrToString: () => Promise<string>; }}
 * @returns {child_process.ChildProcess & ProcessUtiil}
 */
function execAsync(cmd, options) {
	console[execAsync.log]("exec:", cmd);

	const proc = options ? child_process.exec(cmd, options) : child_process.exec(cmd);

	const forward_err = new Error(cmd);// forward error
	
	const promise = new Promise((resolve, reject) => {
		proc.on("error", err => {
			// @ts-ignore
			forward_err.src = err.stack;
			reject(forward_err);
		});
		proc.on("exit", (code, signal) => resolve(proc));
	});

	return ProcessUtiil.assign(proc, promise);
}

spawnAsync.log = "log";
/**
 * @param {string} cmd
 * @param {string[]} args
 * @param {child_process.SpawnOptions} [options]
 * @returns {child_process.ChildProcess & ProcessUtiil}
 */
function spawnAsync(cmd, args, options) {
	try {
		console[spawnAsync.log]("spawnAsync:", [cmd, ...args.map(v => `'${v}'`)].join(" "));

		const proc = options ? child_process.spawn(cmd, args, options) : child_process.spawn(cmd, args);

		const forward_err = new Error([cmd, ...args].join(" "));// forward error

		const exit_promise =  new Promise((resolve, reject) => {
			proc.on("error", err => {
				// @ts-ignore
				forward_err.src = err.stack;
				reject(forward_err);
			});
			proc.on("exit", resolve);
		});

		return ProcessUtiil.assign(proc, exit_promise);
	}
	catch (ex) {
		console[spawnAsync.log](ex.message, cmd, args.join(" "));
		console[spawnAsync.log](ex.stack);
		throw ex;
	}
}

/**
 * @param {import("stream").Readable} stream
 * @returns {Promise<Buffer>}
 */
async function streamToBuffer(stream) {
	if (!stream) {
		throw new TypeError();
	}
	return new Promise((resolve, reject) => {
		const chunks = [];
		stream.on("error", (err) => reject(err));
		stream.on("data", (chunk) => chunks.push(Buffer.from(chunk)));
		stream.on("end", () => resolve(Buffer.concat(chunks)));
	});
}

/**
 * @param {import("stream").Readable} stream
 * @param {any} [encoding]
 * @returns {Promise<string>}
 */
async function streamToString(stream, encoding) {
	return (await streamToBuffer(stream)).toString(encoding);
}

/**
 * @template T
 * @template K
 * @param {T[]} arr
 * @param {K extends keyof T ? K : never} attr
 * @returns {Map<T[K extends keyof T ? K : never], T[]>}
 */
function groupBy(arr, attr) {
	/** @type {Map<T[K extends keyof T ? K : never], T[]>} */
	const group = new Map();

	arr.forEach(obj => {
		const k = obj[attr];
		let gs = group.get(k);
		if (!gs) {
			gs = [];
			group.set(k, gs);
		}
		gs.push(obj);
	});

	return group;
}

module.exports.groupBy = groupBy;

module.exports.execAsync = execAsync;
module.exports.spawnAsync = spawnAsync;
module.exports.streamToBuffer = streamToBuffer;
module.exports.streamToString = streamToString;
