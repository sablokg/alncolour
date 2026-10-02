use std::collections::HashMap;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::process;

/*
Gaurav Sablok
gsablok@proton.me
 */

#[derive(Debug, Clone)]
pub struct Record {
    pub id: String,
    pub seq: Vec<char>,
}

#[tokio::main]
pub async fn alignmenthtml(inputpath: &str, outputpath: &str) {
    let input_path = inputpath.to_string();
    let output_path = outputpath.to_string();

    let raw = match fs::read_to_string(&input_path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Error reading '{}': {}", input_path, e);
            process::exit(1);
        }
    };

    let records = match parse_fasta(&raw) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error parsing FASTA: {}", e);
            process::exit(1);
        }
    };

    if records.is_empty() {
        eprintln!("No sequences found in '{}'.", input_path);
        process::exit(1);
    }

    let aln_len = records[0].seq.len();
    for r in &records {
        if r.seq.len() != aln_len {
            eprintln!(
                "Error: sequences are not aligned (same length). '{}' has length {} but expected {}.",
                r.id,
                r.seq.len(),
                aln_len
            );
            process::exit(1);
        }
    }

    let html = build_html(&records, aln_len);

    let file = match File::create(&output_path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error creating '{}': {}", output_path, e);
            process::exit(1);
        }
    };
    let mut writer = BufWriter::new(file);
    writer
        .write_all(html.as_bytes())
        .expect("failed writing HTML");

    println!(
        "Wrote color-coded alignment ({} sequences x {} columns) to '{}'",
        records.len(),
        aln_len,
        output_path
    );
}

pub fn parse_fasta(raw: &str) -> Result<Vec<Record>, String> {
    let mut records: Vec<Record> = Vec::new();
    let mut current_id: Option<String> = None;
    let mut current_seq: Vec<char> = Vec::new();

    for (lineno, line) in raw.lines().enumerate() {
        let line = line.trim_end();
        if line.trim().is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                records.push(Record {
                    id,
                    seq: current_seq.clone(),
                });
                current_seq.clear();
            }
            let id = rest.trim().to_string();
            if id.is_empty() {
                return Err(format!("empty sequence id on line {}", lineno + 1));
            }
            current_id = Some(id);
        } else {
            if current_id.is_none() {
                return Err(format!(
                    "sequence data on line {} before any '>' header",
                    lineno + 1
                ));
            }
            for c in line.chars() {
                if !c.is_whitespace() {
                    current_seq.push(c.to_ascii_uppercase());
                }
            }
        }
    }
    if let Some(id) = current_id.take() {
        records.push(Record {
            id,
            seq: current_seq,
        });
    }
    Ok(records)
}

pub fn color_for(base: char) -> &'static str {
    match base {
        'A' => "#5fbf5f", // green
        'C' => "#5f9be0", // blue
        'G' => "#e0c25f", // gold
        'T' => "#e07f7f", // red/salmon
        'U' => "#e07f7f", // treat RNA U like T
        'N' => "#c9c9c9", // grey - ambiguous
        '-' => "#f2f2f2", // gap - near white
        _ => "#d9a6e0",   // anything else (ambiguity codes etc.) - purple
    }
}

/// Text color chosen for contrast against the background.
pub fn text_color_for(base: char) -> &'static str {
    match base {
        '-' => "#999999",
        _ => "#1a1a1a",
    }
}

pub fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

pub fn consensus(records: &[Record], aln_len: usize) -> (Vec<char>, Vec<f64>) {
    let mut cons = Vec::with_capacity(aln_len);
    let mut ident = Vec::with_capacity(aln_len);

    for col in 0..aln_len {
        let mut counts: HashMap<char, usize> = HashMap::new();
        for r in records {
            *counts.entry(r.seq[col]).or_insert(0) += 1;
        }
        // pick the most common non-gap symbol if possible, else gap
        let best_non_gap = counts
            .iter()
            .filter(|(c, _)| **c != '-')
            .max_by_key(|(_, count)| **count)
            .map(|(c, count)| (*c, *count));

        let (best_char, best_count) = match best_non_gap {
            Some((c, count)) => (c, count),
            None => ('-', *counts.get(&'-').unwrap_or(&0)),
        };

        cons.push(best_char);
        ident.push(best_count as f64 / records.len() as f64);
    }
    (cons, ident)
}

pub fn build_html(records: &[Record], aln_len: usize) -> String {
    let (cons, ident) = consensus(records, aln_len);

    let max_id_len = records
        .iter()
        .map(|r| r.id.chars().count())
        .max()
        .unwrap_or(4)
        .max(10);

    let mut out = String::new();
    out.push_str("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"UTF-8\">\n");
    out.push_str("<title>DNA Alignment</title>\n<style>\n");
    out.push_str(&format!(
        r#"
:root {{
  --bg: #ffffff;
  --panel: #f7f7f9;
  --border: #dcdce2;
  --text: #1a1a1a;
  --muted: #6b6b76;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  padding: 24px;
  background: var(--bg);
  color: var(--text);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
}}
h1 {{
  font-size: 20px;
  margin: 0 0 4px 0;
}}
.meta {{
  color: var(--muted);
  font-size: 13px;
  margin-bottom: 18px;
}}
.legend {{
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  margin-bottom: 16px;
  font-size: 12px;
}}
.legend span.swatch {{
  display: inline-block;
  width: 12px;
  height: 12px;
  border-radius: 3px;
  margin-right: 4px;
  vertical-align: middle;
  border: 1px solid rgba(0,0,0,0.15);
}}
.legend .item {{
  display: inline-flex;
  align-items: center;
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 999px;
  padding: 3px 10px;
}}
.aln-wrap {{
  overflow-x: auto;
  border: 1px solid var(--border);
  border-radius: 8px;
  background: var(--panel);
  padding: 12px;
}}
.aln {{
  border-collapse: collapse;
  font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
  font-size: 13px;
  white-space: nowrap;
}}
.aln td, .aln th {{
  padding: 0;
  text-align: left;
}}
.seq-id {{
  padding: 1px 10px 1px 4px;
  font-weight: 600;
  color: var(--text);
  position: sticky;
  left: 0;
  background: var(--panel);
  min-width: {max_id_len}ch;
}}
.ruler-id {{
  padding: 1px 10px 1px 4px;
  color: var(--muted);
  position: sticky;
  left: 0;
  background: var(--panel);
  min-width: {max_id_len}ch;
}}
.base {{
  display: inline-block;
  width: 18px;
  text-align: center;
  line-height: 20px;
}}
.ruler-mark {{
  display: inline-block;
  width: 18px;
  text-align: center;
  color: var(--muted);
  font-size: 10px;
}}
.cons-row td {{
  border-top: 2px solid var(--border);
  padding-top: 4px;
}}
.ident-bar {{
  display: inline-block;
  width: 18px;
  height: 6px;
  vertical-align: bottom;
}}
tr:hover .base {{
  filter: brightness(0.92);
}}
footer {{
  margin-top: 14px;
  font-size: 11px;
  color: var(--muted);
}}
"#,
        max_id_len = max_id_len
    ));
    out.push_str("</style>\n</head>\n<body>\n");

    out.push_str("<h1>DNA Alignment</h1>\n");
    out.push_str(&format!(
        "<div class=\"meta\">{} sequences &middot; {} columns</div>\n",
        records.len(),
        aln_len
    ));

    // legend
    out.push_str("<div class=\"legend\">\n");
    for (label, base) in [
        ("A", 'A'),
        ("C", 'C'),
        ("G", 'G'),
        ("T / U", 'T'),
        ("N (ambiguous)", 'N'),
        ("- (gap)", '-'),
    ] {
        out.push_str(&format!(
            "<span class=\"item\"><span class=\"swatch\" style=\"background:{};\"></span>{}</span>\n",
            color_for(base),
            html_escape(label)
        ));
    }
    out.push_str("</div>\n");

    out.push_str("<div class=\"aln-wrap\">\n<table class=\"aln\">\n");

    // Ruler row: mark every 10th column with its 1-based position
    out.push_str("<tr>\n<td class=\"ruler-id\"></td>\n<td>\n");
    for col in 0..aln_len {
        let pos = col + 1;
        if pos % 10 == 0 {
            out.push_str(&format!(
                "<span class=\"ruler-mark\" title=\"position {p}\">{p}</span>",
                p = pos
            ));
        } else {
            out.push_str("<span class=\"ruler-mark\">&middot;</span>");
        }
    }
    out.push_str("</td>\n</tr>\n");

    // Sequence rows
    for r in records {
        out.push_str("<tr>\n");
        out.push_str(&format!(
            "<td class=\"seq-id\">{}</td>\n",
            html_escape(&r.id)
        ));
        out.push_str("<td>");
        for &c in &r.seq {
            out.push_str(&format!(
                "<span class=\"base\" style=\"background:{};color:{};\" title=\"{}\">{}</span>",
                color_for(c),
                text_color_for(c),
                html_escape(&r.id),
                c
            ));
        }
        out.push_str("</td>\n</tr>\n");
    }
    out.push_str("<tr class=\"cons-row\">\n<td class=\"seq-id\">Consensus</td>\n<td>");
    for &c in &cons {
        out.push_str(&format!(
            "<span class=\"base\" style=\"background:{};color:{};\">{}</span>",
            color_for(c),
            text_color_for(c),
            c
        ));
    }
    out.push_str("</td>\n</tr>\n");

    out.push_str("<tr>\n<td class=\"seq-id\">Conservation</td>\n<td>");
    for &frac in &ident {
        let shade = 210 - (frac * 170.0) as i32; // 210 (light) .. 40 (dark)
        let shade = shade.clamp(40, 220);
        out.push_str(&format!(
            "<span class=\"ident-bar\" style=\"background:rgb({s},{s},{s});\" title=\"{pct:.0}% identity\"></span>",
            s = shade,
            pct = frac * 100.0
        ));
    }
    out.push_str("</td>\n</tr>\n");

    out.push_str("</table>\n</div>\n");
    out.push_str("<footer>Generated by dna_align_html (Rust)</footer>\n");
    out.push_str("</body>\n</html>\n");

    out
}
