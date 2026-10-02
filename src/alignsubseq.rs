use std::fmt::Write as FmtWrite;
use std::fs;
use std::path::Path;
use std::process;

struct Record {
    header: String,
    seq: String,
}

/*
Gaurav Sablok
gsablok@proton.me
*/

#[tokio::main]
pub async fn alignsubseq(inputpath: &str, start: &str, end: &str, prefix: &str) {
    let input_path = inputpath;
    let start: usize = start.parse::<usize>().unwrap();
    let end: usize = end.parse::<usize>().unwrap();
    if start > end {
        eprintln!(
            "Error: <start> ({}) cannot be greater than <end> ({}).",
            start, end
        );
        process::exit(1);
    }

    let stem = Path::new(input_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("region");

    let prefix = if prefix.is_empty() {
        format!("{}_{}-{}", stem, start, end)
    } else {
        prefix.to_string()
    };

    let raw = fs::read_to_string(input_path).unwrap_or_else(|e| {
        eprintln!("Error reading '{}': {}", input_path, e);
        process::exit(1);
    });

    let records = parse_fasta(&raw).unwrap_or_else(|e| {
        eprintln!("Error parsing FASTA: {}", e);
        process::exit(1);
    });

    if records.is_empty() {
        eprintln!("Error: no sequences found in '{}'.", input_path);
        process::exit(1);
    }

    // Verify this is actually an alignment: every sequence must be the same length.
    let aln_len = records[0].seq.len();
    for r in &records {
        if r.seq.len() != aln_len {
            eprintln!(
                "Error: sequences are not aligned (expected length {}, '{}' has length {}). \
                 Align your sequences first (e.g. with MAFFT/MUSCLE) before using this tool.",
                aln_len,
                r.header,
                r.seq.len()
            );
            process::exit(1);
        }
    }

    if end > aln_len {
        eprintln!(
            "Error: <end> ({}) is beyond the alignment length ({}).",
            end, aln_len
        );
        process::exit(1);
    }

    // Convert to 0-based half-open slice range.
    let (s0, e0) = (start - 1, end);

    let sliced: Vec<Record> = records
        .iter()
        .map(|r| Record {
            header: r.header.clone(),
            seq: r.seq[s0..e0].to_string(),
        })
        .collect();

    let fasta_out = format!("{}.fasta", prefix);
    let html_out = format!("{}.html", prefix);

    write_fasta(&sliced, &fasta_out).unwrap_or_else(|e| {
        eprintln!("Error writing '{}': {}", fasta_out, e);
        process::exit(1);
    });

    let html = build_html(&sliced, start, end, aln_len, input_path);
    fs::write(&html_out, html).unwrap_or_else(|e| {
        eprintln!("Error writing '{}': {}", html_out, e);
        process::exit(1);
    });

    println!(
        "Extracted columns {}-{} ({} bp) from {} sequences.",
        start,
        end,
        end - start + 1,
        sliced.len()
    );
    println!("  FASTA -> {}", fasta_out);
    println!("  HTML  -> {}", html_out);
}

/// Parse a (possibly multi-line) FASTA string into records, preserving order.
fn parse_fasta(raw: &str) -> Result<Vec<Record>, String> {
    let mut records = Vec::new();
    let mut header: Option<String> = None;
    let mut seq = String::new();

    for (i, line) in raw.lines().enumerate() {
        let line = line.trim_end_matches(['\r', '\n']);
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(h) = header.take() {
                records.push(Record {
                    header: h,
                    seq: std::mem::take(&mut seq),
                });
            }
            header = Some(rest.trim().to_string());
        } else {
            if header.is_none() {
                return Err(format!(
                    "line {}: sequence data found before any '>' header",
                    i + 1
                ));
            }
            seq.push_str(line.trim());
        }
    }
    if let Some(h) = header.take() {
        records.push(Record { header: h, seq });
    }
    Ok(records)
}

/// Write records as wrapped FASTA (60 chars per line).
fn write_fasta(records: &[Record], path: &str) -> std::io::Result<()> {
    let mut out = String::new();
    for r in records {
        writeln!(out, ">{}", r.header).unwrap();
        for chunk in r.seq.as_bytes().chunks(60) {
            out.push_str(std::str::from_utf8(chunk).unwrap());
            out.push('\n');
        }
    }
    fs::write(path, out)
}

/// Map a base to a CSS class for coloring.
fn base_class(c: char) -> &'static str {
    match c.to_ascii_uppercase() {
        'A' => "b-a",
        'T' => "b-t",
        'G' => "b-g",
        'C' => "b-c",
        'N' => "b-n",
        '-' | '.' => "b-gap",
        _ => "b-other",
    }
}

fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// Build a self-contained HTML page showing the extracted region, colored by base,
/// with a position ruler and a summary of where the region sits in the full alignment.
fn build_html(
    records: &[Record],
    start: usize,
    end: usize,
    aln_len: usize,
    source: &str,
) -> String {
    let region_len = end - start + 1;

    let mut rows = String::new();
    for r in records {
        let mut bases = String::new();
        for c in r.seq.chars() {
            let _ = write!(
                bases,
                "<span class=\"{}\">{}</span>",
                base_class(c),
                html_escape(&c.to_string())
            );
        }
        let _ = write!(
            rows,
            "<tr><th class=\"name\">{}</th><td class=\"seqcell\">{}</td></tr>\n",
            html_escape(&r.header),
            bases
        );
    }

    // Ruler: tick marks every 10 columns, labeled with alignment position.
    let mut ruler = String::new();
    for i in 0..region_len {
        let pos = start + i;
        if pos % 10 == 0 || i == 0 {
            let _ = write!(
                ruler,
                "<span class=\"tick\" style=\"left:{}ch\">{}</span>",
                i, pos
            );
        }
    }

    format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Alignment region {start}-{end}</title>
<style>
  :root {{
    --a: #2e7d32; --t: #c62828; --g: #37474f; --c: #1565c0;
    --n: #9e9e9e; --gap: #cfd8dc; --other: #8e24aa;
    --bg: #fafafa; --panel: #ffffff; --border: #e0e0e0; --text: #263238;
  }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    background: var(--bg); color: var(--text); margin: 0; padding: 24px;
  }}
  h1 {{ font-size: 18px; margin: 0 0 4px 0; }}
  .meta {{ color: #607d8b; font-size: 13px; margin-bottom: 18px; }}
  .meta code {{ background: #eceff1; padding: 1px 5px; border-radius: 4px; }}
  .panel {{
    background: var(--panel); border: 1px solid var(--border); border-radius: 10px;
    padding: 16px; overflow-x: auto;
  }}
  table {{ border-collapse: collapse; font-family: "SF Mono", Menlo, Consolas, monospace; font-size: 13px; }}
  th.name {{
    text-align: right; padding: 2px 12px 2px 0; white-space: nowrap;
    font-weight: 600; color: #37474f; position: sticky; left: 0; background: var(--panel);
  }}
  td.seqcell {{ white-space: pre; letter-spacing: 0; line-height: 1.5; }}
  td.seqcell span {{ display: inline-block; width: 1ch; text-align: center; }}
  .b-a {{ color: #fff; background: var(--a); }}
  .b-t {{ color: #fff; background: var(--t); }}
  .b-g {{ color: #fff; background: var(--g); }}
  .b-c {{ color: #fff; background: var(--c); }}
  .b-n {{ color: #fff; background: var(--n); }}
  .b-gap {{ color: #90a4ae; background: var(--gap); }}
  .b-other {{ color: #fff; background: var(--other); }}
  .ruler-row {{ position: relative; height: 18px; margin-bottom: 4px; margin-left: 0; }}
  .ruler {{ position: relative; height: 16px; font-size: 11px; color: #78909c; }}
  .tick {{ position: absolute; transform: translateX(-2px); }}
  .legend {{ margin-top: 16px; font-size: 12px; display: flex; gap: 14px; flex-wrap: wrap; color: #455a64; }}
  .legend span.swatch {{ display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 4px; vertical-align: middle; }}
</style>
</head>
<body>
  <h1>Aligned region {start}&ndash;{end}</h1>
  <div class="meta">
    Source: <code>{source}</code> &middot;
    Region length: <code>{region_len} bp</code> &middot;
    Full alignment length: <code>{aln_len} bp</code> &middot;
    Sequences shown: <code>{n_seqs}</code>
  </div>
  <div class="panel">
    <table>
      <tr>
        <th class="name"></th>
        <td class="seqcell"><div class="ruler-row"><div class="ruler">{ruler}</div></div></td>
      </tr>
      {rows}
    </table>
  </div>
  <div class="legend">
    <div><span class="swatch" style="background:var(--a)"></span>A</div>
    <div><span class="swatch" style="background:var(--t)"></span>T</div>
    <div><span class="swatch" style="background:var(--g)"></span>G</div>
    <div><span class="swatch" style="background:var(--c)"></span>C</div>
    <div><span class="swatch" style="background:var(--n)"></span>N</div>
    <div><span class="swatch" style="background:var(--gap)"></span>gap (-)</div>
  </div>
</body>
</html>
"#,
        start = start,
        end = end,
        source = html_escape(source),
        region_len = region_len,
        aln_len = aln_len,
        n_seqs = records.len(),
        ruler = ruler,
        rows = rows,
    )
}
